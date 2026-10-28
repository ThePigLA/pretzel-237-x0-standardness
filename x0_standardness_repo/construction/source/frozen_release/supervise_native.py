#!/usr/bin/env python3
"""External cumulative CPU / process-group RSS / wall-clock supervisor."""
from pathlib import Path
import argparse,json,os,resource,signal,subprocess,sys,time
import psutil
ROOT=Path(__file__).resolve().parent
CAP_CPU=4*3600
CAP_RSS=8*1024**3

def atomic(p,x):
    t=p.with_suffix('.tmp');t.write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n');t.replace(p)

def main():
    a=argparse.ArgumentParser();a.add_argument('--wall',type=float,default=600);a.add_argument('--cpu',type=float);a.add_argument('--ledger-name',default='native_budget_ledger.json');a.add_argument('--worker-gib',type=float,default=8);a.add_argument('args',nargs=argparse.REMAINDER);x=a.parse_args()
    if not x.args:raise ValueError('worker arguments required')
    folder=ROOT/'control';folder.mkdir(exist_ok=True)
    if Path(x.ledger_name).name!=x.ledger_name:raise ValueError('ledger must be a filename in control/')
    lp=folder/x.ledger_name
    ledger=json.loads(lp.read_text()) if lp.exists() else {'schema':'externally-supervised-native-budget-v1','cpu_cap_seconds':CAP_CPU,'rss_cap_bytes':CAP_RSS,'max_distinct_seeds':16,'runs':[]}
    used=sum(r.get('cpu_seconds',0) for r in ledger['runs']);remaining=CAP_CPU-used
    if x.cpu is not None:remaining=min(remaining,x.cpu)
    rss_limit=min(CAP_RSS,int(x.worker_gib*1024**3))
    if remaining<=0:raise RuntimeError('cumulative CPU cap already consumed')
    if x.args[0]=='search':
        seed=int(x.args[x.args.index('--seed')+1]) if '--seed' in x.args else 0
        seeds={r.get('seed') for r in ledger['runs'] if r.get('mode')=='search'}|{seed}
        if len(seeds)>16:raise RuntimeError('16-seed cap')
    else:seed=None
    run={'mode':x.args[0],'arguments':x.args,'seed':seed,'wall_cap_seconds':x.wall,'cpu_remaining_before':remaining,'status':'RUNNING'}
    ledger['runs'].append(run);atomic(lp,ledger)
    logprefix='native_run' if x.ledger_name=='native_budget_ledger.json' else Path(x.ledger_name).stem+'_run'
    log=folder/f'{logprefix}_{len(ledger["runs"]):02d}.log';start=time.monotonic();ru0=resource.getrusage(resource.RUSAGE_CHILDREN)
    env=dict(os.environ);env['PYTHONUNBUFFERED']='1'
    def enforce_limits():
        resource.setrlimit(resource.RLIMIT_AS,(rss_limit,rss_limit))
        resource.setrlimit(resource.RLIMIT_CPU,(max(1,int(remaining)),max(2,int(remaining)+1)))
    with log.open('w') as fh:
        child=subprocess.Popen([sys.executable,str(ROOT/'native_control.py'),*x.args],stdout=fh,stderr=subprocess.STDOUT,env=env,start_new_session=True,preexec_fn=enforce_limits)
        seen={};peak=0;reason=None;host_pid=None
        while child.poll() is None:
            rss=0
            try:
                # The executor uses a PID namespace with a host /proc mount.
                # Popen.pid is namespace-local, whereas psutil enumerates host
                # pids. Resolve the exact executable/argument vector first.
                if host_pid is None:
                    matches=[]
                    expected=[sys.executable,str(ROOT/'native_control.py'),*x.args]
                    for candidate in psutil.process_iter(['pid','cmdline','create_time']):
                        if candidate.info['cmdline']==expected and candidate.info['create_time']>=time.time()-(time.monotonic()-start)-2:
                            matches.append(candidate.info['pid'])
                    if len(matches)==1:host_pid=matches[0]
                if host_pid is None:
                    processes=[]
                else:
                    parent=psutil.Process(host_pid);processes=[parent,*parent.children(recursive=True)]
                for p in processes:
                    try:
                        c=p.cpu_times();seen[p.pid]=max(seen.get(p.pid,0),c.user+c.system);rss+=p.memory_info().rss
                    except (psutil.NoSuchProcess,psutil.AccessDenied):pass
                peak=max(peak,rss)
            except psutil.NoSuchProcess:pass
            elapsed=time.monotonic()-start;cpu=sum(seen.values())
            if rss>rss_limit:reason='MEMORY_CAP'
            elif cpu>remaining:reason='CUMULATIVE_CPU_CAP'
            elif elapsed>x.wall:reason='WALL_CAP'
            if reason:
                os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=3)
                except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL)
                break
            time.sleep(.25)
        child.wait()
    run.update({'status':reason or ('COMPLETED' if child.returncode==0 else 'WORKER_FAILED'),'exit_code':child.returncode,
                'wall_seconds':time.monotonic()-start,'cpu_seconds':sum(seen.values()),'peak_process_group_rss_bytes':peak,'log':log.name})
    ru1=resource.getrusage(resource.RUSAGE_CHILDREN)
    run['cpu_seconds']=max(run['cpu_seconds'],ru1.ru_utime+ru1.ru_stime-ru0.ru_utime-ru0.ru_stime)
    run['peak_waited_worker_rss_bytes']=int(ru1.ru_maxrss*1024)
    run['peak_process_group_rss_bytes']=max(run['peak_process_group_rss_bytes'],run['peak_waited_worker_rss_bytes'])
    run['hard_kernel_limits']={'address_space_bytes':rss_limit,'worker_cpu_seconds':max(1,int(remaining))}
    run['monitor_pid_resolution']='exact worker executable and argument vector in host /proc; Popen pid may be namespace-local'
    run['cpu_accounting']='max of process-tree sampled user+system CPU and waited-child getrusage user+system CPU'
    atomic(lp,ledger);print(log.read_text(),end='');print(json.dumps(run),flush=True)
    if run['status']!='COMPLETED':sys.exit(2)

if __name__=='__main__':main()

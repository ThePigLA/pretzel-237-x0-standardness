"""Reuse exactly equal X2 labelled vertex links; produce new bistellar certificates for changed X0 links."""
import argparse,gzip,hashlib,importlib.util,json,math,os,random,resource,signal,sys,time
from pathlib import Path
sys.dont_write_bytecode=True

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return json.loads(gzip.decompress(p.read_bytes()))
def atomic_gz(p,obj):
    temp=p.with_name(p.name+'.tmp')
    temp.write_bytes(gzip.compress(json.dumps(obj,separators=(',',':')).encode(),compresslevel=1,mtime=0))
    os.replace(temp,p)
def atomic_json(p,obj):
    temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(obj,indent=2)+'\n');os.replace(temp,p)
def links_of(d):
    F=[tuple(sorted(f)) for f in d['facets']]
    if len(set(F))!=len(F) or not all(len(f)==5 and len(set(f))==5 and all(type(v)==int for v in f) for f in F):raise ValueError('bad simplicial master')
    L={}
    for f in F:
        for v in f:L.setdefault(v,[]).append(tuple(x for x in f if x!=v))
    return {v:sorted(ts) for v,ts in L.items()}
class BudgetExceeded(Exception):pass
def alarm(signum,frame):raise BudgetExceeded('per-link or global CPU limit')
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--x0-master',type=Path,required=True);p.add_argument('--x2-master',type=Path,required=True)
p.add_argument('--x2-certificate',type=Path,required=True);p.add_argument('--checker',type=Path,required=True)
p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
p.add_argument('--cpu-seconds',type=float,default=120);p.add_argument('--memory-gib',type=float,default=4)
a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
if sys.flags.optimize:raise RuntimeError('Run without Python optimization; replay assertions are proof checks.')
start=time.process_time();wall=time.monotonic();used=math.ceil(start+a.cpu_seconds)
resource.setrlimit(resource.RLIMIT_AS,(int(a.memory_gib*1024**3),int(a.memory_gib*1024**3)))
resource.setrlimit(resource.RLIMIT_CPU,(used,used+2))
signal.signal(signal.SIGPROF,alarm);signal.signal(signal.SIGXCPU,alarm)
spec=importlib.util.spec_from_file_location('reviewer_link_checker',a.checker);checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)
D0=rd(a.x0_master);D2=rd(a.x2_master);L0=links_of(D0);L2=links_of(D2);C2=rd(a.x2_certificate)
if set(map(str,L2))!=set(C2['links']):raise ValueError('X2 certificate vertices mismatch')
reuse=[];changed=[];result={}
for v,L in sorted(L0.items()):
    if v in L2 and L==L2[v]:reuse.append(v);result[str(v)]=C2['links'][str(v)]
    else:changed.append(v)
binding={'x0_master_file_sha256':sha(a.x0_master),'x0_master_facets_sha256':D0['facets_sha256'],'x2_master_file_sha256':sha(a.x2_master),'x2_certificate_file_sha256':sha(a.x2_certificate),'checker_file_sha256':sha(a.checker)}
new=[];failed=[];attempts=[];status='PARTIAL';reason=None
certpath=out/'x0_vertex_links.json.gz'
def save():
    cert={'schema':'x0-vertex-links-bistellar-v1','vertices':len(L0),'links':result,'flips_total':sum(map(len,result.values())),'largest_link_tets':max(map(len,L0.values())),'complete':len(result)==len(L0),'input_bindings':binding,'reused_vertices':reuse,'new_vertices':new,'failed_vertices':failed,'reuse_rule':'exact equality of sorted labelled tetrahedron lists'}
    atomic_gz(certpath,cert)
    receipt={'status':'PRODUCED_AWAITING_FULL_REPLAY' if cert['complete'] else 'PARTIAL','input_bindings':binding,'vertices':len(L0),'source_facets':len(D0['facets']),'exactly_reused_links':len(reuse),'changed_links':len(changed),'new_links_produced_and_locally_replayed':len(new),'links_completed':len(result),'failed_vertices':failed,'new_vertex_records':attempts,'new_flips':sum(len(result[str(v)]) for v in new),'total_flips':cert['flips_total'],'cpu_seconds':time.process_time()-start,'wall_seconds':time.monotonic()-wall,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'cpu_ceiling_seconds':a.cpu_seconds,'address_space_ceiling_bytes':int(a.memory_gib*1024**3),'termination_reason':reason,'certificate_file_sha256':sha(certpath),'certificate_bytes':certpath.stat().st_size,'full_replay_passed':False}
    atomic_json(out/'production_receipt.json',receipt)
    return receipt
print(json.dumps({'event':'LINK_COMPARISON','vertices':len(L0),'unchanged':len(reuse),'changed':len(changed),'changed_tetrahedra':sum(len(L0[v]) for v in changed),'largest_changed_link':max(len(L0[v]) for v in changed)}),flush=True)
save()
# Do small links first so any difficult cap-center link cannot erase completed work.
for v in sorted(changed,key=lambda v:(len(L0[v]),v)):
    remain=a.cpu_seconds-5-(time.process_time()-start)
    if remain<=0:reason='CPU budget reserve reached';break
    link_start=time.process_time();log=None;failure=None
    perlink=min(remain,30 if len(L0[v])>1000 else 2)
    signal.setitimer(signal.ITIMER_PROF,perlink)
    try:
        rng=random.Random(v)
        for attempt in range(3):
            log=checker.reduce3(L0[v],rng,max_steps=5000)
            if log is not None:break
        if log is not None and not checker.verify3(L0[v],log):raise ValueError('new trace failed immediate replay')
    except BudgetExceeded as e:failure=str(e);log=None
    finally:signal.setitimer(signal.ITIMER_PROF,0)
    if log is None:failed.append(v);failure=failure or 'bounded reducer did not reach boundary of 4-simplex'
    else:result[str(v)]=log;new.append(v)
    attempts.append({'vertex':v,'tetrahedra':len(L0[v]),'flips':None if log is None else len(log),'passed':log is not None,'cpu_seconds':time.process_time()-link_start,'failure':failure})
    if (len(new)+len(failed))%25==0 or len(L0[v])>1000:
        save();print(json.dumps({'event':'PROGRESS','new':len(new),'failed':len(failed),'remaining':len(changed)-len(new)-len(failed),'cpu_seconds':time.process_time()-start}),flush=True)
if len(result)==len(L0):reason='All vertex-link sequences produced; full replay remains separate.'
elif reason is None:reason='Some bounded link reductions did not complete.'
r=save();print(json.dumps({k:v for k,v in r.items() if k!='new_vertex_records'},indent=2),flush=True)

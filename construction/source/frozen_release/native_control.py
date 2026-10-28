#!/usr/bin/env python3
"""Native admission and explicitly replayed simplification of a marked control.

Runs only on saved numbered gluing tables. Marks remain in the immutable master.
Native simplify internals are an explicit software trust boundary. This does not
pretend a numbered macro receipt is a portable elementary Pachner transcript.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,random,resource,sys,threading,time,warnings
import regina
warnings.filterwarnings('ignore',category=DeprecationWarning)
ROOT=Path(__file__).resolve().parent
from torus_adaptor import load,save
WORK_SOURCE=False
MACRO_POLICY='full'
INITIAL_LOCAL=False
RESUME=False

def filehash(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def native(d):return regina.Triangulation4.fromGluings(d['simplex_count'],[(i,k,j,regina.Perm5(*p)) for i,k,j,p in d['gluings']])
def numbered(T):
    rows=[[(p.adjacentPentachoron(f).index(),str(p.adjacentGluing(f))) if p.adjacentPentachoron(f) else None for f in range(5)] for p in T.pentachora()]
    return hashlib.sha256(json.dumps(rows).encode()).hexdigest()
def simple(T,policy=None):
    policy=policy or MACRO_POLICY
    if policy=='local':
        regina.RandomEngine.reseedWithDefault();T.simplifyToLocalMinimum();return
    regina.RandomEngine.reseedWithDefault();tracker=regina.ProgressTrackerObjective(T.size());done=threading.Event()
    def observe():
        while not done.wait(30):
            print('NATIVE_PROGRESS',tracker.objective(),tracker.description(),flush=True)
    observer=threading.Thread(target=observe,daemon=True);observer.start()
    try:T.simplify(tracker)
    finally:done.set()
def atomic(p,d):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n');tmp.replace(p)
def source_data(which='cover'):
    if WORK_SOURCE and which=='cover':
        cfg=load(ROOT/'control/native_work_source.json')
        p=ROOT/cfg['relative_path']
        if filehash(p)!=cfg['source_file_sha256']:raise ValueError('certified work source changed')
        for name,h in cfg['prefix_file_hashes'].items():
            if filehash(ROOT/name)!=h:raise ValueError('unmarked PL prefix data changed')
    else:p=ROOT/'control'/f'{which}_regina.json.gz'
    d=load(p);return p,d
def admit(which):
    st=time.monotonic();cpu=time.process_time();p,d=source_data(which);T=native(d)
    r={'schema':'native-marked-control-admission-v1','source':p.name,'source_sha256':filehash(p),'regina':regina.versionString(),
       'simplex_count':T.size(),'valid':T.isValid(),'orientable':T.isOrientable(),'connected':T.isConnected(),
       'closed':T.isClosed(),'ideal':T.isIdeal(),'boundary_components':T.countBoundaryComponents(),
       'f_vector':[T.countVertices(),T.countEdges(),T.countTriangles(),T.countTetrahedra(),T.size()],
       'chi':T.eulerCharTri(),'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu,
       'vertex_link_scope':'Regina isValid and isClosed/isIdeal test closed vertex links and recognize sphere links; no homology substitution for link recognition.'}
    r['passed']=r['valid'] and r['orientable'] and r['connected'] and r['closed'] and not r['ideal'] and r['chi']==2
    if WORK_SOURCE:r['unmarked_PL_prefix']=load(ROOT/'control/native_work_source.json')
    atomic(ROOT/'control'/f'{which}{"_work" if WORK_SOURCE else ""}_native_admission.json',r);print(json.dumps(r),flush=True)
    if not r['passed']:raise ValueError('native admission failed')

def invariants(which):
    st=time.monotonic();cpu=time.process_time();p,d=source_data(which);T=native(d);print('INVARIANT_INPUT',p.name,T.size(),flush=True);simple(T)
    print('INVARIANT_SIMPLIFIED',T.size(),round(time.monotonic()-st,2),flush=True)
    h=numbered(T);g=T.group();before=str(g);g.simplify()
    r={'schema':'native-unmarked-invariant-crosscheck-v1','source':p.name,'source_sha256':filehash(p),'regina':regina.versionString(),
       'macro':'RandomEngine.reseedWithDefault(); '+('simplifyToLocalMinimum()' if MACRO_POLICY=='local' else 'simplify()'),'macro_policy':MACRO_POLICY,'endpoint_numbered_sha256':h,
       'endpoint_pentachora':T.size(),'endpoint_isoSig':T.isoSig(),'valid':T.isValid(),'closed':T.isClosed(),
       'H1':str(T.homology()),'H2':str(T.homology(2)),'chi':T.eulerCharTri(),
       'group_before_simplification':before,'group_after_simplification':str(g),'group_generators':g.countGenerators(),
       'trust_boundary':'Native PL simplify and native group simplification. These crosschecks alone do not recognize a sphere.'}
    T2=native(d);simple(T2);r['macro_replay_match']=numbered(T2)==h
    r['wall_seconds']=time.monotonic()-st;r['cpu_seconds']=time.process_time()-cpu
    r['passed']=r['macro_replay_match'] and r['valid'] and r['closed'] and r['H1']=='0' and r['H2']=='0' and r['group_generators']==0
    if WORK_SOURCE:r['unmarked_PL_prefix']=load(ROOT/'control/native_work_source.json')
    atomic(ROOT/'control'/f'{which}_invariants.json',r);print(json.dumps(r),flush=True)
    if not r['passed']:raise ValueError('control invariants/replay did not pass')

def search(seed,seconds,max_trials):
    st=time.monotonic();cpu=time.process_time();p,d=source_data();rng=random.Random(seed)
    out=ROOT/'control'/f'standardization_seed{seed}.json'
    # Every source pentachoron is an embedded simplex of the saved abstract
    # simplicial complex. Deleting its interior creates one standard S3 boundary.
    source=load(ROOT/'control/marked_control.json.gz')
    work=load(ROOT/'control/native_work_source.json') if WORK_SOURCE else None
    working_complex=load(ROOT/work['facets_relative_path']) if WORK_SOURCE else source
    if RESUME:
        cert=load(out)
        if cert['source_sha256']!=filehash(p) or cert['source_facets_sha256']!=source['facets_sha256'] or cert['unmarked_PL_prefix']!=work:raise ValueError('resume source/prefix mismatch')
        puncture=cert['puncture'];pent=working_complex['facets'][puncture]
        if cert['source_puncture_vertices']!=pent or len(set(pent))!=5:raise ValueError('resume puncture mismatch')
        T=regina.Triangulation4.tightDecoding(cert['endpoint_tight_encoding'])
        if not cert['blocks'] or numbered(T)!=cert['blocks'][-1]['hash'] or T.size()!=cert['blocks'][-1]['pentachora']:raise ValueError('resume endpoint differs from last certified block')
        # A bare decoded 4D triangulation loses Regina's known-manifold
        # properties. Its next skeleton query can re-run difficult 3-sphere
        # recognition. Restore through the actual legal history instead;
        # this also replays the preceding proof before extending it.
        T=native(d);T.removePentachoron(T.pentachoron(puncture))
        for bi,b in enumerate(cert['blocks']):
            for i in b['kicks']:
                if not T.pachner(T.tetrahedron(i),True,True):raise ValueError('illegal preceding resume move')
            simple(T,b.get('macro_policy',cert.get('macro_policy','full')))
            if numbered(T)!=b['hash'] or T.size()!=b['pentachora']:raise ValueError('resume history replay mismatch')
            print('RESUME_PREFIX_REPLAYED',seed,bi+1,T.size(),flush=True)
        phases=cert.setdefault('resume_history',[]);rng_seed=seed+1000003*(len(phases)+1);rng=random.Random(rng_seed)
        phases.append({'starting_blocks':len(cert['blocks']),'starting_pentachora':T.size(),'starting_numbered_hash':numbered(T),
                       'previous_phase_wall_seconds':cert['wall_seconds'],'previous_phase_cpu_seconds':cert['cpu_seconds'],
                       'rng_seed':rng_seed,'requested_phase_seconds':seconds,'restore_method':'full preceding native history replayed from original work source',
                       'reason':'continue an atomic completed-macro checkpoint; incomplete trial is discarded'})
        print('RESUME',seed,len(cert['blocks']),T.size(),flush=True)
    else:
        T=native(d);puncture=rng.randrange(T.size());pent=working_complex['facets'][puncture]
        if len(set(pent))!=5:raise ValueError('source puncture not embedded')
        T.removePentachoron(T.pentachoron(puncture))
        cert={'schema':'source-bound-native-macro-standardization-v1','source':p.name,'source_sha256':filehash(p),
          'source_facets_sha256':source['facets_sha256'],'puncture':puncture,'source_puncture_vertices':pent,'unmarked_PL_prefix':work,
          'regina_version':regina.versionString(),'seed':seed,'blocks':[],'macro_policy':MACRO_POLICY,
          'limits':{'wall_seconds_this_seed':seconds,'trials':max_trials,'max_pentachora':250000},
          'trust_boundary':'Guarded native 2-4 moves and deterministic native simplify macros; every accepted macro endpoint is hashed and replayed.'}
    def checkpoint():
        cert.update({'final_size':T.size(),'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu})
        cert['endpoint_tight_encoding']=T.tightEncoding()
        atomic(out,cert)
    if not RESUME:
        first_policy='local' if INITIAL_LOCAL else MACRO_POLICY
        simple(T,first_policy);cert['blocks'].append({'kicks':[],'hash':numbered(T),'pentachora':T.size(),'macro_policy':first_policy});checkpoint()
        print('INITIAL',seed,puncture,T.size(),round(time.monotonic()-st,2),flush=True)
    else:
        checkpoint()
        prefix={'schema':'native-resume-prefix-replay-v1','source_sha256':filehash(p),'certificate_sha256':filehash(out),
                'recorded_blocks_sha256':hashlib.sha256(json.dumps(cert['blocks'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                'blocks':len(cert['blocks']),'guarded_kicks':sum(len(b['kicks']) for b in cert['blocks']),
                'endpoint_pentachora':T.size(),'endpoint_numbered_sha256':numbered(T),'replay_passed':True,
                'standardization_passed':False,'ball':False,'endpoint_isoSig':None,
                'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu,
                'scope':'Replay of the complete saved prefix at resume time; subsequent accepted blocks require their own final replay.'}
        atomic(ROOT/'control'/f'standardization_seed{seed}_resume_prefix_replay.json',prefix)
    best=T.size();stall=0;trials=0
    while T.size()>1 and time.monotonic()-st<seconds and trials<max_trials:
        U=regina.Triangulation4(T);kicks=[];trials+=1
        print('TRIAL',seed,trials,'source',T.size(),flush=True)
        for _ in range(rng.randint(1,10 if stall<60 else 20)):
            ix=[i for i in range(U.countTetrahedra()) if not U.tetrahedron(i).isBoundary()];rng.shuffle(ix)
            for i in ix:
                if U.pachner(U.tetrahedron(i),True,True):kicks.append(i);break
        if U.size()>250000:raise RuntimeError('pentachoron cap during native search')
        print('MACRO_START',seed,trials,len(kicks),U.size(),flush=True)
        simple(U);sideways=stall>=60 and U.size()<=T.size()+4
        if U.size()<T.size() or (sideways and rng.random()<0.5):
            T=U;cert['blocks'].append({'kicks':kicks,'hash':numbered(T),'pentachora':T.size(),'macro_policy':MACRO_POLICY})
            if T.size()<best:
                best=T.size();stall=0;print('REDUCED',seed,best,len(cert['blocks']),round(time.monotonic()-st,2),flush=True)
            else:stall=max(0,stall-30)
            checkpoint()
        else:stall+=1
    cert['trials']=trials;cert['ball']=T.size()==1 and T.isIsomorphicTo(regina.Example4.ball()) is not None
    cert['final_isoSig']=T.isoSig() if T.size()<=250 else None;checkpoint()
    print('SEARCH_DONE',seed,T.size(),cert['ball'],cert['final_isoSig'],flush=True)

def replay(seed):
    global MACRO_POLICY
    st=time.monotonic();cpu=time.process_time();out=ROOT/'control'/f'standardization_seed{seed}.json';c=json.loads(out.read_text());p,d=source_data()
    MACRO_POLICY=c.get('macro_policy','full')
    if c['source_sha256']!=filehash(p):raise ValueError('source changed')
    if c['regina_version']!=regina.versionString():raise ValueError('Regina API version differs from saved certificate')
    source=load(ROOT/'control/marked_control.json.gz')
    work=load(ROOT/'control/native_work_source.json') if WORK_SOURCE else None
    if c['source_facets_sha256']!=source['facets_sha256'] or c['unmarked_PL_prefix']!=work:raise ValueError('replay master/prefix mismatch')
    working_complex=load(ROOT/work['facets_relative_path']) if WORK_SOURCE else source
    puncture=c['puncture']
    if not isinstance(puncture,int) or not 0<=puncture<len(working_complex['facets']):raise ValueError('puncture index out of range')
    pent=working_complex['facets'][puncture]
    if c['source_puncture_vertices']!=pent or len(set(pent))!=5:raise ValueError('replay puncture is not the saved embedded simplex')
    if not c['blocks']:raise ValueError('empty native replay path')
    T=native(d);T.removePentachoron(T.pentachoron(c['puncture']))
    for bi,b in enumerate(c['blocks']):
        for i in b['kicks']:
            if not T.pachner(T.tetrahedron(i),True,True):raise ValueError('illegal recorded move')
        simple(T,b.get('macro_policy',MACRO_POLICY))
        if numbered(T)!=b['hash'] or T.size()!=b['pentachora']:raise ValueError('macro endpoint mismatch')
        print('REPLAY_BLOCK',seed,bi+1,T.size(),flush=True)
    if T.size()!=c['final_size'] or T.tightEncoding()!=c['endpoint_tight_encoding']:raise ValueError('saved final checkpoint differs from replay')
    ball=T.size()==1 and T.isIsomorphicTo(regina.Example4.ball()) is not None
    r={'schema':'native-control-standardization-replay-v1','source_sha256':filehash(p),'certificate_sha256':filehash(out),
       'source_facets_sha256':source['facets_sha256'],'source_puncture_vertices':pent,'puncture':puncture,
       'native_program_sha256':filehash(Path(__file__)),'regina_version':regina.versionString(),
       'endpoint_numbered_sha256':numbered(T),'endpoint_tight_encoding_match':True,
       'blocks':len(c['blocks']),'guarded_kicks':sum(len(b['kicks']) for b in c['blocks']),
       'endpoint_pentachora':T.size(),'endpoint_isoSig':T.isoSig() if T.size()<=250 else None,'ball':ball,
       'replay_passed':True,'standardization_passed':ball,'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu}
    atomic(ROOT/'control'/f'standardization_seed{seed}_replay.json',r);print(json.dumps(r),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['admit','invariants','search','replay']);a.add_argument('--which',choices=['cover','quotient'],default='cover');a.add_argument('--seed',type=int,default=0);a.add_argument('--seconds',type=float,default=120);a.add_argument('--trials',type=int,default=1000);a.add_argument('--work-source',action='store_true');a.add_argument('--local',action='store_true');a.add_argument('--initial-local',action='store_true');a.add_argument('--resume',action='store_true');x=a.parse_args();WORK_SOURCE=x.work_source;MACRO_POLICY='local' if x.local else 'full';INITIAL_LOCAL=x.initial_local;RESUME=x.resume
    if x.mode=='admit':admit(x.which)
    elif x.mode=='invariants':invariants(x.which)
    elif x.mode=='search':search(x.seed,x.seconds,x.trials)
    else:replay(x.seed)

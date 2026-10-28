"""Fresh complete X0 vertex-link replay, using only the Python standard library."""
import argparse,gzip,hashlib,importlib.util,json,os,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode=True

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def read(p):return json.loads(gzip.decompress(p.read_bytes()))
base=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--master',type=Path,default=base.parent/'candidate/marked_candidate.json.gz')
p.add_argument('--certificate',type=Path,default=base/'x0_vertex_links.json.gz')
p.add_argument('--checker',type=Path,default=base/'check_links.py')
p.add_argument('--receipt',type=Path,default=base/'replay_receipt.json')
a=p.parse_args()
require(not sys.flags.optimize,'Run without Python optimization: assertions are proof checks.')
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_CPU,(120,122))
t=time.process_time();wall=time.monotonic()
C=read(a.certificate);D=read(a.master)
require(C.get('complete') is True,'certificate is partial')
require(sha(a.master)==C['input_bindings']['x0_master_file_sha256'],'master bytes differ')
expected_checker_sha='aa6a86a25e6a79975cc3d90076e2861e2e407ec216a903b26331ac9e0bb2638c'
require(sha(a.checker)==C['input_bindings']['checker_file_sha256']==expected_checker_sha,'checker bytes differ from the audited implementation')
F=[tuple(f) for f in D['facets']]
require(len(set(F))==len(F) and all(len(f)==5 and len(set(f))==5 and all(type(v)==int for v in f) for f in F),'invalid master facets')
V=set().union(*map(set,F))
require(set(C['links'])==set(map(str,V)) and C['vertices']==len(V),'vertex coverage differs')
flips=0
for vertex,log in C['links'].items():
    require(type(log)==list,'invalid trace')
    for step in log:
        require(type(step)==list and len(step)==2,'invalid flip record')
        S,T=step
        require(type(S)==list and type(T)==list and S and T,'empty flip simplex')
        require(all(type(v)==int for v in S+T),'invalid vertex labels')
        require(len(set(S))==len(S) and len(set(T))==len(T) and not set(S)&set(T) and len(S)+len(T)==5,'invalid flip schema')
        flips+=1
require(flips==C['flips_total'],'total move count differs')
# This independently reruns all actual local preconditions from the frozen X0 master.
spec=importlib.util.spec_from_file_location('reviewer_link_checker',a.checker);checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)
checker.verify_saved(str(a.master),str(a.certificate))
r={'status':'PASS','master_file_sha256':sha(a.master),'master_facets_sha256':D['facets_sha256'],'certificate_file_sha256':sha(a.certificate),'checker_file_sha256':sha(a.checker),'master_facets':len(F),'vertices':len(V),'links_replayed':len(V),'bistellar_flips_replayed':flips,'every_vertex_link_reaches_boundary_of_4_simplex':True,'all_move_schemas_checked':True,'complete_source_vertex_coverage':True,'assertions_enabled':True,'cpu_seconds':time.process_time()-t,'wall_seconds':time.monotonic()-wall,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'address_space_ceiling_bytes':4*1024**3,'cpu_ceiling_seconds':120,'native_engine_used':False,'conclusion':'The stored X0 cover master is a closed combinatorial 4-manifold. This is a manifold certificate, not a standardization of X0.'}
a.receipt.parent.mkdir(parents=True,exist_ok=True);temp=a.receipt.with_name(a.receipt.name+'.tmp');temp.write_text(json.dumps(r,indent=2)+'\n');os.replace(temp,a.receipt)
print(json.dumps(r,indent=2),flush=True)

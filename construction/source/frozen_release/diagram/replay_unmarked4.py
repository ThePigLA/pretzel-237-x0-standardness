#!/usr/bin/env python3
"""Independent exact replay of closed PL4 full-link edge contractions.

Does not import the producer or reuse its link cache.  Every step checks all
simplices of both links, oriented replacement boundaries, and actual full-face
intersection with the retained complex.  The PL-manifold premise comes from
the source-bound native admission and is inherited through Nevo's theorem.
"""
from collections import Counter,defaultdict,deque
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,gzip,json,resource,time
ROOT=Path(__file__).resolve().parents[1]

def req(x,msg):
    if not x:raise ValueError(msg)
def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def filehash(p):return sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    b=Path(p).read_bytes();return json.loads(gzip.decompress(b) if str(p).endswith('.gz') else b)
def save(p,x):
    b=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(p).write_bytes(gzip.compress(b,mtime=0) if str(p).endswith('.gz') else b)
def full_faces(fs):
    return {s for f in fs for k in range(1,len(f)+1) for s in combinations(f,k)}
def face_counts(fs):
    return Counter(s for f in fs for k in range(1,len(f)+1) for s in combinations(f,k))
def boundary_chain(orient):
    c=Counter()
    for f,sign in orient.items():
        for i in range(len(f)):c[f[:i]+f[i+1:]]+=sign*(-1)**i
    return {f:s for f,s in c.items() if s}
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def replay(folder):
    started=time.monotonic();cpu=time.process_time()
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    source=load(folder/'source_facets.json.gz');cert=load(folder/'certificate.json.gz');endpoint=load(folder/'endpoint.json.gz')
    req(cert['unmarked_only'] and cert['schema']=='closed-PL4-full-link-edge-contractions-v1','Certificate scope')
    req(cert['source_snapshot_sha256']==digest(source),'Source snapshot hash')
    master_path=ROOT/source['source_file'];master=load(master_path)
    req(filehash(master_path)==source['source_file_sha256'],'Immutable source file changed')
    req(master['facets']==source['facets'] and master['orientation']==source['orientation'],'Snapshot differs from marked source')
    req(digest(source['facets'])==cert['source_facets_sha256']==source['source_facets_sha256'],'Source facets hash')
    native_path=ROOT/'control/cover_regina.json.gz';admission_path=ROOT/'control/cover_native_admission.json'
    req(filehash(native_path)==source['source_native_file_sha256'],'Native source changed')
    req(filehash(admission_path)==source['native_admission_file_sha256'],'Native admission changed')
    admission=load(admission_path);native=load(native_path)
    req(admission['source_sha256']==filehash(native_path),'Admission source binding')
    req(all(admission[k] for k in ('passed','valid','closed','connected','orientable')) and not admission['ideal'],'Closed PL source premise absent')
    req(native['source_facets_sha256']==source['source_facets_sha256'],'Admitted geometry differs')
    live=dict(zip(map(tuple,source['facets']),source['orientation']))
    req(len(live)==len(source['facets']) and all(tuple(sorted(set(f)))==f and len(f)==5 for f in live),'Initial abstract simplicial complex invalid')
    req(set(live.values())<={-1,1} and not boundary_chain(live),'Initial orientation chain not closed')
    counts=face_counts(live);stars=defaultdict(set)
    for f in live:
        for v in f:stars[v].add(f)
    req(all(n==2 for f,n in counts.items() if len(f)==4),'Initial codimension-one incidence')
    stats=Counter();last_notice=time.monotonic()
    for index,step in enumerate(cert['steps']):
        u,v=step['u'],step['v'];req(u!=v and stars[u] and stars[v],'Contraction endpoint missing')
        req(len(live)==step['before'],'Recorded before count')
        old=set(stars[u]);req(any(v in f for f in old),'Not an edge')
        Ltop={tuple(w for w in f if w!=u) for f in old}
        L=full_faces(Ltop)
        V=full_faces({tuple(w for w in f if w!=v) for f in stars[v]})
        Etop={tuple(w for w in f if w not in (u,v)) for f in old if v in f}
        E=full_faces(Etop)
        req(L&V==E,'Complete simplicial link condition fails')
        req([sum(len(f)==k for f in E) for k in (1,2,3)]==step['edge_link_f_vector'],'Recorded edge-link f-vector')
        antistar={f for f in Ltop if v not in f}
        req(antistar,'Empty link antistar')
        new={}
        for f in old:
            if v in f:continue
            q=tuple(sorted((v,*(w for w in f if w!=u))))
            req(q not in new,'Two old simplices produce one new simplex')
            new[q]=live[f]*(-1)**(f.index(u)+q.index(v))
        req(set(new)=={tuple(sorted((v,*f))) for f in antistar},'Not the stated cone on the antistar')
        req(not(set(new)&(set(live)-old)),'New top cell already in retained complex')
        old_orientation={f:live[f] for f in old}
        old_boundary=boundary_chain(old_orientation);new_boundary=boundary_chain(new)
        req(set(old_boundary)==Ltop and set(old_boundary.values())<={-1,1},'Old star boundary not the complete vertex link')
        req(new_boundary==old_boundary,'Oriented replacement boundary differs')
        old_counts=face_counts(old);new_counts=face_counts(new)
        new_top_faces=Counter(q for f in new for q in combinations(f,4))
        req(all(c in (1,2) for c in new_top_faces.values()),'New ball has bad facet incidence')
        req({q for q,c in new_top_faces.items() if c==1}==Ltop,'New ball has wrong boundary triangulation')
        for face in new_counts:
            outside=counts[face]-old_counts[face]
            req(outside>=0,'Negative retained face multiplicity')
            req(not outside or face in L,'New ball meets retained complex away from common boundary')
        req(len(old)==step['deleted'] and len(new)==step['inserted'],'Replacement counts differ')
        req(digest(sorted(old))==step['old_star_sha256'] and digest(sorted(new))==step['new_star_sha256'],'Local star digest differs')
        req(len(new)<len(old),'Pentachora did not strictly decrease')
        for face,c in old_counts.items():
            counts[face]-=c
            if not counts[face]:del counts[face]
        counts.update(new_counts)
        for f in old:
            del live[f]
            for w in f:stars[w].remove(f)
        for f,sign in new.items():
            live[f]=sign
            for w in f:stars[w].add(f)
        req(len(live)==step['after'],'Recorded after count')
        stats['full_link_conditions_checked']+=1;stats['relative_4ball_replacements']+=1
        stats['full_face_overlap_checks']+=len(new_counts)
        if time.monotonic()-last_notice>=10:
            print(f'replayed {index+1}/{len(cert["steps"])} contractions; {len(live)} pentachora',flush=True)
            last_notice=time.monotonic()
    fs=sorted(live);signs=[live[f] for f in fs]
    req(fs==list(map(tuple,endpoint['facets'])) and signs==endpoint['orientation'],'Endpoint differs from replay')
    req(digest(fs)==endpoint['facets_sha256']==cert['endpoint_facets_sha256'],'Endpoint facet hash')
    req(endpoint['vertices']==sorted({v for f in fs for v in f}),'Endpoint vertex labels')
    req(not boundary_chain(live),'Endpoint orientation is not closed')
    req(all(c==2 for f,c in counts.items() if len(f)==4),'Endpoint not closed degree-two facets')
    fvec=[sum(len(f)==k for f in counts) for k in range(1,6)]
    req(sum((-1)**i*n for i,n in enumerate(fvec))==2,'Endpoint Euler characteristic changed')
    pending={};gluings=[];adj=defaultdict(set)
    for i,f in enumerate(fs):
        for k in range(5):
            face=f[:k]+f[k+1:]
            if face not in pending:pending[face]=(i,k)
            else:
                j,l=pending.pop(face)
                p=[fs[j].index(w) if a!=k else l for a,w in enumerate(f)]
                req(signs[j]==-signs[i]*parity(p),'Endpoint native orientation mismatch')
                gluings.append([i,k,j,p]);adj[i].add(j);adj[j].add(i)
    req(not pending,'Endpoint native gluing incomplete')
    seen={0};todo=deque([0])
    while todo:
        i=todo.popleft()
        for j in adj[i]:
            if j not in seen:seen.add(j);todo.append(j)
    req(len(seen)==len(fs),'Endpoint disconnected')
    exported={'schema':'source-bound-unmarked-PL4-contraction-native-endpoint-v1','simplex_count':len(fs),'gluings':gluings,
              'source_facets_sha256':endpoint['facets_sha256'],'original_marked_facets_sha256':source['source_facets_sha256'],
              'contraction_certificate_file_sha256':filehash(folder/'certificate.json.gz')}
    save(folder/'endpoint_regina.json.gz',exported)
    result={'schema':'independent-closed-PL4-contraction-replay-v1','passed':True,
            'source_facets_sha256':source['source_facets_sha256'],'endpoint_facets_sha256':endpoint['facets_sha256'],
            'source_snapshot_file_sha256':filehash(folder/'source_facets.json.gz'),
            'certificate_file_sha256':filehash(folder/'certificate.json.gz'),
            'endpoint_file_sha256':filehash(folder/'endpoint.json.gz'),
            'endpoint_native_file_sha256':filehash(folder/'endpoint_regina.json.gz'),
            'source_pentachora':len(source['facets']),'endpoint_pentachora':len(fs),'endpoint_f_vector':fvec,
            'contractions':len(cert['steps']),'checks':dict(stats),
            'PL_premise':'source-bound native recognition of all closed vertex links as PL3-spheres, inherited through the full-link theorem',
            'theorem':'Nevo Theorem1.4, https://arxiv.org/pdf/math/0602531; Section5.1 gives identical-boundary PL ball replacement',
            'orientation_preserved':True,'marking_scope':'unmarked descendant; no deck or marked surface transfer is claimed',
            'wall_seconds':time.monotonic()-started,'cpu_seconds':time.process_time()-cpu,
            'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    save(folder/'independent_replay.json',result);print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--folder',type=Path,default=ROOT/'unmarked_contractions_v1')
    replay(ap.parse_args().folder)

#!/usr/bin/env python3
"""Bounded, unmarked PL4 contractions using the complete link condition.

The marked master is read only. No deck, branch sphere, twin or chart is
asserted on the descendant. Independent replay is required before use.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,gzip,heapq,json,resource,time
ROOT=Path(__file__).resolve().parents[1]

def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def filehash(p):return sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    b=Path(p).read_bytes();return json.loads(gzip.decompress(b) if str(p).endswith('.gz') else b)
def save(p,x):
    b=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(p).write_bytes(gzip.compress(b,mtime=0) if str(p).endswith('.gz') else b)

class Complex:
    def __init__(self,facets,signs):
        self.live=dict(zip(map(tuple,facets),signs));self.stars=defaultdict(set)
        self.revision=defaultdict(int);self.cache={};self.attempts=0
        for f in self.live:
            for v in f:self.stars[v].add(f)
    def link(self,v):
        if v not in self.cache:
            top={tuple(w for w in f if w!=v) for f in self.stars[v]}
            self.cache[v]=[{w for f in top for w in f},
                           {e for f in top for e in combinations(f,2)},
                           {t for f in top for t in combinations(f,3)},top]
        return self.cache[v]
    def proposal(self,u,v):
        self.attempts+=1
        if u==v or not self.stars[u] or not self.stars[v]:return None
        U,V=self.link(u),self.link(v)
        if v not in U[0]:return None
        edge_top={tuple(w for w in f if w!=u and w!=v) for f in self.stars[u] if v in f}
        E=[{w for f in edge_top for w in f},{e for f in edge_top for e in combinations(f,2)},edge_top,set()]
        if any(U[k]&V[k]!=E[k] for k in range(4)):return None
        old=set(self.stars[u]);new={}
        for f in old:
            if v in f:continue
            q=tuple(sorted((v,*(w for w in f if w!=u))))
            if q in new or (q in self.live and q not in old):return None
            new[q]=self.live[f]*(-1)**(f.index(u)+q.index(v))
        if not new or len(new)>=len(old):return None
        return old,new,[len(z) for z in E[:3]]
    def apply(self,u,v,proposal):
        old,new,ef=proposal;before=len(self.live)
        affected={w for f in old|set(new) for w in f}
        record={'u':u,'v':v,'before':before,'after':before-len(old)+len(new),
                'deleted':len(old),'inserted':len(new),'edge_link_f_vector':ef,
                'old_star_sha256':digest(sorted(old)),'new_star_sha256':digest(sorted(new))}
        for f in old:
            del self.live[f]
            for w in f:self.stars[w].remove(f)
        for f,sign in new.items():
            self.live[f]=sign
            for w in f:self.stars[w].add(f)
        for w in affected:self.revision[w]+=1;self.cache.pop(w,None)
        return affected,record
    def run(self,deadline):
        queue=[];records=[]
        def push(v):
            if self.stars[v]:heapq.heappush(queue,(len(self.stars[v]),v,self.revision[v]))
        for v in sorted(self.stars):push(v)
        reason='local_priority_queue_exhausted';last_notice=time.monotonic()
        while queue:
            if time.monotonic()>=deadline:reason='declared_time_cap';break
            _,u,rev=heapq.heappop(queue)
            if self.revision[u]!=rev or not self.stars[u]:continue
            ns=sorted(self.link(u)[0],key=lambda v:(-len(self.stars[u]&self.stars[v]),-len(self.stars[v]),v))
            for v in ns:
                if time.monotonic()>=deadline:reason='declared_time_cap';break
                p=self.proposal(u,v)
                if p is None:continue
                affected,record=self.apply(u,v,p);records.append(record)
                for w in sorted(affected):push(w)
                break
            if time.monotonic()-last_notice>=10:
                print(f'{len(records)} contractions; {len(self.live)} pentachora; {sum(bool(s) for s in self.stars.values())} vertices',flush=True)
                last_notice=time.monotonic()
        return records,reason

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--scan-only',action='store_true')
    ap.add_argument('--seconds',type=float,default=55)
    ap.add_argument('--output',type=Path,default=ROOT/'unmarked_contractions_v1');args=ap.parse_args()
    if not 0<args.seconds<=55:raise ValueError('Initial search/preparation cap is 55s; reserve output time within60s')
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    began=time.monotonic();cpu=time.process_time();source_path=ROOT/'control/marked_control.json.gz'
    source=load(source_path);admission_path=ROOT/'control/cover_native_admission.json';adm=load(admission_path)
    if not all(adm[k] for k in ('passed','valid','closed','connected','orientable')) or adm['ideal']:raise ValueError('No closed PL-manifold admission')
    if adm['source_sha256']!=filehash(ROOT/'control/cover_regina.json.gz'):raise ValueError('Admission source changed')
    native=load(ROOT/'control/cover_regina.json.gz')
    if native['source_facets_sha256']!=source['facets_sha256'] or digest(source['facets'])!=source['facets_sha256']:raise ValueError('Admitted facets differ')
    worker=Complex(source['facets'],source['orientation'])
    if args.scan_only:
        sample=[]
        for u in sorted(worker.stars,key=lambda v:(len(worker.stars[v]),v))[:24]:
            accepted=[]
            for v in sorted(worker.link(u)[0]):
                p=worker.proposal(u,v)
                if p is not None:accepted.append({'v':v,'pentachora_removed':len(p[0])-len(p[1])})
            sample.append({'u':u,'star_size':len(worker.stars[u]),'admissible_neighbors':accepted})
        print(json.dumps({'source_pentachora':len(worker.live),'source_vertices':len(worker.stars),
                          'sampled_vertices':len(sample),'admissible_edges_in_sample':sum(len(z['admissible_neighbors']) for z in sample),
                          'sample':sample,'seconds':time.monotonic()-began},indent=2));return
    args.output.mkdir(exist_ok=False,parents=True)
    snapshot={'facets':source['facets'],'orientation':source['orientation'],
              'source_file':'control/marked_control.json.gz','source_file_sha256':filehash(source_path),
              'source_facets_sha256':source['facets_sha256'],'source_native_file_sha256':filehash(ROOT/'control/cover_regina.json.gz'),
              'native_admission_file_sha256':filehash(admission_path)}
    save(args.output/'source_facets.json.gz',snapshot)
    moves,reason=worker.run(began+args.seconds);facets=sorted(worker.live)
    endpoint={'schema':'unmarked-closed-4D-edge-contraction-endpoint-v1','facets':facets,
              'orientation':[worker.live[f] for f in facets],'facets_sha256':digest(facets),
              'vertices':sorted({v for f in facets for v in f}),
              'scope':'Unmarked oriented PL descendant only; no inherited deck, branch sphere, twin or source chart.'}
    cert={'schema':'closed-PL4-full-link-edge-contractions-v1','source_snapshot_sha256':digest(snapshot),
          'source_facets_sha256':source['facets_sha256'],'endpoint_facets_sha256':endpoint['facets_sha256'],'steps':moves,
          'theorem':{'author':'Eran Nevo','title':"Higher Minors and Van Kampen's Obstruction",
                     'reference':'Theorem1.4; proof Section5.1','url':'https://arxiv.org/pdf/math/0602531',
                     'premise':'compact closed PL manifold; full link condition in all dimensions',
                     'local_rule':'replace u*lk(u) by v*ast(v,lk(u)) along identical boundary lk(u)'},
          'unmarked_only':True,'stop_reason':reason,'attempts':worker.attempts,
          'budget':{'search_and_preparation_seconds':args.seconds,'declared_total_seconds':60,
                    'address_space_bytes':8*1024**3,'pentachora_never_increase':True}}
    save(args.output/'endpoint.json.gz',endpoint);save(args.output/'certificate.json.gz',cert)
    result={'source_pentachora':len(source['facets']),'endpoint_pentachora':len(facets),'source_vertices':len(source['vertices']),
            'endpoint_vertices':len(endpoint['vertices']),'contractions':len(moves),'attempts':worker.attempts,'stop_reason':reason,
            'wall_seconds':time.monotonic()-began,'cpu_seconds':time.process_time()-cpu,
            'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            'source_facets_sha256':source['facets_sha256'],'endpoint_facets_sha256':endpoint['facets_sha256'],
            'scope':'PRODUCER_ONLY_PENDING_INDEPENDENT_REPLAY'}
    save(args.output/'producer_summary.json',result);print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()

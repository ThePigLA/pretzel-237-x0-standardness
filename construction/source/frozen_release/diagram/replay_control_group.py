#!/usr/bin/env python3
"""Independent standard-library replay of the immutable master's pi_1 proof.

No import of either producer, group_cells, Regina, or a group algebra package.
The full source skeleton, spanning tree presentation, and every exact Tietze
operation are reconstructed.  Native admission is used only for the separate
closed-manifold premise in the integral homology deduction.
"""
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip,json,resource,time
ROOT=Path(__file__).resolve().parents[1]

def req(condition,message):
    if not condition:raise ValueError(message)
def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def filehash(p):return sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    b=Path(p).read_bytes();return json.loads(gzip.decompress(b) if str(p).endswith('.gz') else b)
def inverse(w):return tuple(-a for a in reversed(w))
def reduce_word(w):
    stack=[]
    for a in w:
        if stack and stack[-1]==-a:stack.pop()
        else:stack.append(a)
    first=0;last=len(stack)
    while last-first>1 and stack[first]==-stack[last-1]:first+=1;last-=1
    return tuple(stack[first:last])

class Presentation:
    def __init__(self,generators,relators):
        self.active=set(generators);self.relations=set();self.occurrences=defaultdict(set)
        self.eliminations=0;self.shortenings=0;self.rewritten_relators=0
        for w in relators:self.insert(tuple(w))
    def insert(self,w):
        req(w==reduce_word(w),'Relator not reduced as declared')
        req(all(abs(a) in self.active for a in w),'Inactive generator in relator')
        if not w or w in self.relations:return
        self.relations.add(w)
        for g in set(map(abs,w)):self.occurrences[g].add(w)
    def remove(self,w):
        req(w in self.relations,'Unavailable relator')
        self.relations.remove(w)
        for g in set(map(abs,w)):self.occurrences[g].remove(w)
    def eliminate(self,event):
        g=event['generator'];r=tuple(event['relator']);s=tuple(event['replacement'])
        req(g in self.active and r in self.relations,'Defining relation/generator unavailable')
        positions=[j for j,a in enumerate(r) if abs(a)==g]
        req(len(positions)==1,'Generator does not occur exactly once')
        j=positions[0];tail=r[j+1:]+r[:j]
        expected=inverse(tail) if r[j]==g else tail
        req(s==expected and all(abs(a) in self.active-{g} for a in s),'Incorrect defining-word solution')
        affected=tuple(self.occurrences[g]);req(r in affected,'Missing pivot in occurrence index')
        si=inverse(s);new=[]
        for w in affected:
            if w!=r:new.append(reduce_word(a for z in w for a in (s if z==g else si if z==-g else (z,))))
        for w in affected:self.remove(w)
        req(not self.occurrences[g],'Eliminated generator has an unprocessed relation')
        self.active.remove(g)
        for w in new:self.insert(w)
        self.eliminations+=1;self.rewritten_relators+=len(new)
    def shorten(self,event):
        r=tuple(event['rule_relator']);w=tuple(event['target'])
        req(r in self.relations and w in self.relations and r!=w,'Rule must be retained while another relator changes')
        sign=event['rule_sign'];j=event['rule_start'];p=event['match_start'];n=event['match_length']
        req(sign in (-1,1) and 0<=j<len(r) and 0<=p<len(w),'Rule orientation or rotation')
        q=r if sign==1 else inverse(r);q=q[j:]+q[:j]
        rotated=w[p:]+w[:p]
        req(len(q)//2<n<=len(q) and n<=len(w),'Rule is not a strict shortening')
        req(rotated[:n]==q[:n],'Rule subword does not match')
        # From q=U V=1, substitute U by inverse(V), keeping q's defining
        # relation r. This is invertible modulo that retained relation.
        result=reduce_word(inverse(q[n:])+rotated[n:])
        req(result==tuple(event['result']) and len(result)<len(w),'Incorrect shortened relation')
        self.remove(w);self.insert(result)
        req(r in self.relations,'Defining relation was lost')
        self.shortenings+=1;self.rewritten_relators+=1

def main():
    started=time.monotonic();cpu=time.process_time()
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    folder=ROOT/'group_audit';binding=load(folder/'source_binding.json')
    master_path=ROOT/binding['source_file'];master=load(master_path)
    req(filehash(master_path)==binding['source_file_sha256'],'Immutable master file changed')
    facets=list(map(tuple,master['facets']))
    req(len(set(facets))==len(facets) and all(len(f)==5 and tuple(sorted(set(f)))==f for f in facets),'Source is not the stated abstract simplicial complex')
    req(digest(master['facets'])==master['facets_sha256']==binding['source_facets_sha256'],'Source facet binding')
    vertices=sorted({v for f in facets for v in f})
    edges=sorted({e for f in facets for e in combinations(f,2)})
    triangles=sorted({t for f in facets for t in combinations(f,3)})
    sk=load(folder/'skeleton.json.gz')
    req(sk['vertices']==vertices and list(map(tuple,sk['edges']))==edges and list(map(tuple,sk['triangles']))==triangles,'Stored skeleton omits or changes a source cell')
    req(digest(sk)==binding['skeleton_sha256'] and sk['source_facets_sha256']==master['facets_sha256'],'Skeleton hash binding')
    initial=load(folder/'initial_presentation.json.gz')
    req(digest(initial)==binding['initial_presentation_sha256'],'Initial presentation binding')
    tree=initial['tree_edges'];req(len(tree)==len(set(tree))==len(vertices)-1,'Tree edge count')
    req(all(isinstance(i,int) and 0<=i<len(edges) for i in tree),'Tree edge outside source')
    parent={v:v for v in vertices}
    def root(v):
        while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
        return v
    for i in tree:
        a,b=edges[i];a,b=root(a),root(b);req(a!=b,'Tree contains a cycle');parent[a]=b
    req(len({root(v) for v in vertices})==1 and initial['root'] in parent,'Tree does not span connected source')
    tree=set(tree);generator_edges=[i for i in range(len(edges)) if i not in tree]
    req(initial['generator_edge_ids']==generator_edges,'Generator-to-edge correspondence differs')
    req(initial['generators']==len(generator_edges)==len(edges)-len(vertices)+1,'Incorrect number of edge generators')
    generator={edges[e]:j+1 for j,e in enumerate(generator_edges)}
    source_relators=set()
    for a,b,c in triangles:
        w=reduce_word(z for z in (generator.get((a,b),0),generator.get((b,c),0),-generator.get((a,c),0)) if z)
        if w:source_relators.add(w)
    req(sorted(source_relators)==list(map(tuple,initial['relators'])),'Initial presentation is not all source triangle boundaries')
    prefix_path=folder/'tietze_certificate.json.gz';prefix=load(prefix_path)
    req(prefix['source_binding_sha256']==digest(binding) and prefix['initial_presentation_sha256']==digest(initial),'Prefix input binding')
    state=Presentation(range(1,initial['generators']+1),source_relators);notice=time.monotonic()
    for i,event in enumerate(prefix['eliminations']):
        state.eliminate(event)
        if time.monotonic()-notice>=10:
            print(f'replayed {i+1}/{len(prefix["eliminations"])} prefix eliminations; {len(state.active)} generators',flush=True)
            notice=time.monotonic()
    req(sorted(state.active)==prefix['generators'] and sorted(state.relations)==list(map(tuple,prefix['relators'])),'Prefix residual presentation differs')
    residual_path=folder/'residual_tietze_certificate.json.gz';residual=load(residual_path)
    req(residual['prefix_file_sha256']==filehash(prefix_path),'Tail not bound to exact prefix')
    req(residual['initial_generators']==prefix['generators'] and residual['initial_relators_sha256']==digest(prefix['relators']),'Tail initial state')
    for event in residual['events']:
        if event['kind']=='eliminate':state.eliminate(event)
        elif event['kind']=='relator_shortening':state.shorten(event)
        else:raise ValueError('Unknown Tietze action')
    req(sorted(state.active)==residual['generators'] and sorted(state.relations)==list(map(tuple,residual['relators'])),'Final residual presentation differs')
    req(all(abs(a) in state.active for w in state.relations for a in w),'Final inactive generator survives')
    trivial=not state.active and not state.relations
    # The following checks serve only the manifold/homology corollary.  The
    # exact source pi_1 presentation and replay above require no native tool.
    counts=Counter();oriented_boundary=Counter()
    req(len(master['orientation'])==len(facets) and set(master['orientation'])<={-1,1},'Source orientation array')
    for f,sign in zip(facets,master['orientation']):
        for j in range(5):
            face=f[:j]+f[j+1:];counts[face]+=1;oriented_boundary[face]+=sign*(-1)**j
    req(all(c==2 for c in counts.values()) and not any(oriented_boundary.values()),'Source not closed with the supplied orientation')
    fvec=[len(vertices),len(edges),len(triangles),len(counts),len(facets)]
    chi=sum((-1)**i*n for i,n in enumerate(fvec));req(chi==2,'Source Euler characteristic')
    admission_path=ROOT/'control/cover_native_admission.json';native_path=ROOT/'control/cover_regina.json.gz'
    req(filehash(admission_path)==binding['native_admission_file_sha256'] and filehash(native_path)==binding['native_file_sha256'],'Native manifold premise binding changed')
    admission=load(admission_path);native=load(native_path)
    req(admission['source_sha256']==filehash(native_path) and native['source_facets_sha256']==master['facets_sha256'],'Manifold admission geometry binding')
    req(all(admission[k] for k in ('passed','valid','closed','connected','orientable')) and not admission['ideal'],'Closed oriented PL-manifold premise unavailable')
    req(admission['f_vector']==fvec and admission['chi']==chi,'Admitted cell counts differ')
    result={'schema':'independent-immutable-control-group-replay-v1','passed':True,
            'source_file_sha256':filehash(master_path),'source_facets_sha256':master['facets_sha256'],
            'source_f_vector':fvec,'source_chi':chi,'initial_generators':initial['generators'],
            'initial_triangle_relators':len(source_relators),'one_occurrence_eliminations':state.eliminations,
            'retained_relator_shortenings':state.shortenings,'rewritten_relators_checked':state.rewritten_relators,
            'remaining_generators':len(state.active),'remaining_relators':len(state.relations),
            'fundamental_group':'trivial' if trivial else 'residual presentation certified; triviality unresolved',
            'integral_homology':{'H0':'Z','H1':'0','H2':'0','H3':'0','H4':'Z'} if trivial else None,
            'homology_derivation':'pi1=1 gives H1=0. For the admitted connected closed oriented 4-manifold, PD/UCT give H3=0 and torsion-free H2. chi=2 forces rank H2=0.',
            'manifold_premise':'The native source-bound closed PL-manifold admission is inherited here; group replay itself is pure standard-library and independent of native simplification.',
            'scope':'Exact pi1 and integral homology certificate for the immutable control; no PL standardization is asserted.',
            'file_sha256':{p.name:filehash(p) for p in [folder/'source_binding.json',folder/'skeleton.json.gz',
                           folder/'initial_presentation.json.gz',prefix_path,residual_path,admission_path]},
            'wall_seconds':time.monotonic()-started,'cpu_seconds':time.process_time()-cpu,
            'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    (folder/'independent_replay.json').write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()

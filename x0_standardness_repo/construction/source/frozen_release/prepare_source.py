"""Cut the actual fixed-circle tube, retaining local source carriers.

The level function is the sum of barycentric coordinates of fixed vertices.
N={f>=1/2}, E={f<=1/2}, T={f=1/2}. Pulling is equivariant.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict,deque
from functools import lru_cache
import json,gzip,sys,hashlib,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from group_cells import closure,digest,presentation,eliminate,free,inv,rewrite

def read(p):
    p=Path(p);b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def save(p,x):
    p=Path(p);p.parent.mkdir(exist_ok=True,parents=True)
    b=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
    p.write_bytes(gzip.compress(b,mtime=0) if p.suffix=='.gz' else b)
def fvec(fs):
    C=closure(fs);return [sum(len(t)==k for t in C) for k in range(1,5)]
def circle(edges):
    adj=defaultdict(list)
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    if not adj or any(len(v)!=2 for v in adj.values()):raise ValueError('not a circle')
    a=min(adj);path=[a,min(adj[a])]
    while path[-1]!=a:
        v=path[-1];path.append(next(w for w in adj[v] if w!=path[-2]))
        if len(path)>len(adj)+1:raise ValueError('circle traversal')
    if len(path)!=len(adj)+1:raise ValueError('disconnected circle')
    return path

def cut_source(s):
    Y=list(map(tuple,s['facets']));d=s['deck'];r=[min(i,d[i]) for i in range(len(d))]
    F={i for i,j in enumerate(d) if i==j};fs=sorted(closure(Y),key=lambda t:(len(t),t))
    if any(len({r[v] for v in t})!=4 for t in Y):raise ValueError('folded deck cell')
    if max(len(t) for t in fs if set(t)<=F)!=2:raise ValueError('fixed locus dimension')
    core=circle([f for f in fs if len(f)==2 and set(f)<=F])
    pts=sorted({('v',v) for t in Y for v in t}|{('m',*e) for e in fs if len(e)==2 and (e[0] in F)!=(e[1] in F)})
    pi={p:i for i,p in enumerate(pts)}
    D=[pi[('v',d[p[1]])] if p[0]=='v' else pi[('m',*sorted(d[v] for v in p[1:]))] for p in pts]
    def ptkey(p):
        if p[0]=='v':return (0,r[p[1]],-1)
        a,b=p[1:];a,b=(a,b) if a in F else (b,a)
        return (1,r[a],r[b])
    cells={};by={}
    for f in fs:
        a=[v for v in f if v in F];b=[v for v in f if v not in F]
        ms=[pi[('m',*sorted((v,w)))] for v in a for w in b]
        for typ,vs,dim in [('N',[pi[('v',v)] for v in a]+ms,len(f)-1),('E',[pi[('v',v)] for v in b]+ms,len(f)-1),('T',ms,len(f)-2)]:
            if not vs:continue
            vs=tuple(sorted(vs));cells[vs]=(typ,f,dim);by[typ,f]=vs
    bd={}
    for vs,(typ,f,k) in cells.items():
        bs=set()
        if k:
            for g in combinations(f,len(f)-1):
                z=by.get((typ,g))
                if z and cells[z][2]==k-1:bs.add(z)
            if typ!='T' and ('T',f) in by:bs.add(by['T',f])
        bd[vs]=sorted(bs)
    @lru_cache(None)
    def pull(vs):
        if len(vs)==1:return (vs,)
        keys=[ptkey(pts[v]) for v in vs]
        if len(set(keys))!=len(keys):raise ValueError('pulling rank collision')
        v=min(vs,key=lambda v:ptkey(pts[v]))
        return tuple(sorted(tuple(sorted((v,*q))) for f in bd[vs] if v not in f for q in pull(f)))
    out={'schema':'fixed-circle-tube-cut-v1','source_sha256':digest(s),
         'level':'sum(fixed-vertex barycentric coordinates)=1/2','points':pts,'deck':D,
         'fixed_source_vertices':sorted(F),'source_ranks':r,'oriented_fixed_circle_source_vertices':core,
         'orientation_scope':'chosen circle orientation; comparison to the signed diagram recorded separately'}
    for typ,k in [('N',3),('E',3),('T',2)]:
        out[typ]=sorted(q for vs,(side,f,dim) in cells.items() if side==typ and dim==k for q in pull(vs))
    out['top_carriers']=[{'side':side,'source_tetrahedron_index':i,'children':list(pull(by[side,t]))}
                        for i,t in enumerate(Y) for side in ('N','E') if (side,t) in by]
    out['f_vectors']={k:fvec(out[k]) for k in ('N','E','T')}
    for side in ('N','E'):
        cc=Counter(f for t in out[side] for f in combinations(t,3))
        if {f for f,n in cc.items() if n==1}!=set(out['T']) or any(n not in (1,2) for n in cc.values()):raise ValueError('cut boundary')
        if {tuple(sorted(D[v] for v in t)) for t in out[side]}!=set(out[side]):raise ValueError('cut deck')
    edge=tuple(pi['v',v] for v in core[:2])
    link={tuple(v for v in t if v not in edge) for t in out['N'] if set(edge)<=set(t)}
    eq=circle(link)
    n=len(eq)-1
    if n%2 or any(D[eq[i]]!=eq[(i+n//2)%n] for i in range(n)):raise ValueError('normal half turn')
    punct=[t for t in out['N'] if set(edge)<=set(t)]
    bdp=sorted(tuple(sorted((v,*f))) for v in edge for f in link)
    out['puncture']={'fixed_edge_cut_vertices':edge,'source_fixed_edge':core[:2],'equator_cycle':eq[:-1],
                     'removed_N_tetrahedra':punct,'boundary_triangles':bdp,
                     'model':'edge join its normal link circle; boundary is S0 join S1',
                     'normal_orientation':'minimum-labelled neighbor chosen; final sign transport is separate'}
    return out

def skeleton(fs,D,quotient=True):
    C=closure(fs)
    rep=lambda t:min(t,tuple(sorted(D[v] for v in t))) if quotient else t
    rows=[sorted({rep(t) for t in C if len(t)==k}) for k in (1,2,3)]
    ei={e:i for i,e in enumerate(rows[1])}
    def edge(a,b):
        z=rep(tuple(sorted((a,b))));sg=1 if (a,b)==z or (quotient and (D[a],D[b])==z) else -1
        return sg*(ei[z]+1)
    return {'vertices':[v[0] for v in rows[0]],'edge_representatives':rows[1],
            'edge_endpoints':[(min(a,D[a]),min(b,D[b])) if quotient else (a,b) for a,b in rows[1]],
            'triangle_representatives':rows[2],
            'triangle_boundaries':[[edge(a,b),edge(b,c),edge(c,a)] for a,b,c in rows[2]]}

def export3(facets,D=None):
    fs=sorted(map(tuple,facets));K=set(fs);D=list(D) if D is not None else None
    rep=lambda f:min(f,tuple(sorted(D[v] for v in f))) if D is not None else f
    reps=sorted({rep(f) for f in fs});idx={t:i for i,t in enumerate(reps)}
    ridges=defaultdict(list)
    for t in fs:
        for k in range(4):ridges[t[:k]+t[k+1:]].append((t,k))
    glu=[]
    for i,t in enumerate(reps):
        for k in range(4):
            h=ridges[t[:k]+t[k+1:]]
            if len(h)==1:continue
            if len(h)!=2:raise ValueError('nonmanifold ridge')
            u,l=next((u,l) for u,l in h if (u,l)!=(t,k))
            R=rep(u);j=idx[R]
            trans=(lambda v:v) if u==R else lambda v:D[v]
            target_l=R.index(trans(u[l]))
            p=[R.index(trans(v)) if a!=k else target_l for a,v in enumerate(t)]
            if (i,k)<=(j,target_l):glu.append([i,k,j,p])
    return {'schema':'source-bound-regina3-gluings-v1','simplex_count':len(reps),'gluings':glu,
            'source_facets_sha256':digest(fs),'tetrahedron_representatives':reps,'quotient_by_deck':D is not None}

def native(x):
    import regina
    t=regina.Triangulation3();t.newTetrahedra(x['simplex_count'])
    for i,k,j,p in x['gluings']:t.tetrahedron(i).join(k,t.tetrahedron(j),regina.Perm4(*p))
    return t

def run(source,out):
    st=time.monotonic();s=read(source);x=cut_source(s);out=Path(out)
    save(out/'cut_source.json.gz',x);print('CUT',x['f_vectors'],'equator',len(x['puncture']['equator_cycle']),flush=True)
    report={'source_sha256':digest(s),'cut_sha256':digest(x),'f_vectors':x['f_vectors'],'cases':{}}
    for side in ('N','E'):
        for quotient in (False,True):
            label=side+('_quotient' if quotient else '_cover')
            z=export3(x[side],x['deck'] if quotient else None);save(out/(label+'_regina.json.gz'),z)
            t=native(z)
            import regina
            small=regina.Triangulation3(t);small.simplify()
            r={'valid':t.isValid(),'orientable':t.isOrientable(),'connected':t.isConnected(),
               'tetrahedra':t.size(),'boundary_components':t.countBoundaryComponents(),
               'H1':str(small.homology()),'euler':t.eulerCharManifold(),
               'simplified_tetrahedra':small.size(),'simplified_isoSig':small.isoSig(),
               'homology_scope':'computed on a native unmarked simplification copy'}
            if not(r['valid'] and r['orientable'] and r['connected'] and r['H1']=='Z' and r['boundary_components']==1):raise ValueError((label,r))
            report['cases'][label]=r
            print(label,r,flush=True)
            if quotient and side=='E':
                t=small
                sigs=[{'stage':'simplified_compact','size':t.size(),'isoSig':t.isoSig()}]
                t.finiteToIdeal();sigs.append({'stage':'finiteToIdeal','size':t.size()})
                t.simplify();sigs.append({'stage':'simplify','size':t.size(),'isoSig':t.isoSig()})
                print('EXTERIOR ideal',t.size(),t.isoSig(),flush=True)
                save(out/'quotient_exterior_native_path.json',{'scope':'Regina-native macro ancestry, internal legal moves trusted','states':sigs})
                save(out/'quotient_exterior_ideal.json',{'isoSig':t.isoSig(),'simplex_count':t.size()})
    report['seconds']=time.monotonic()-st
    save(out/'cut_native_admission.json',report)
    return x

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=ROOT/'coarsened/source.json.gz');p.add_argument('--output',type=Path,default=ROOT/'peripheral')
    a=p.parse_args();run(a.source,a.output)

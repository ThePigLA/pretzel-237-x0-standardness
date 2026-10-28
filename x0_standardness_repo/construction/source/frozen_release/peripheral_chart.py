"""Exact peripheral cohomology and rational flat-torus charts on the source cut.

No presentation from the target knot is inserted into the extraction.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict,deque
from itertools import combinations
import math,time,sys
from prepare_source import read,save,skeleton,ROOT
from group_cells import presentation,eliminate,free,inv,rewrite,digest

def abelian_cocycle(sk):
    import sympy as S
    p,_,_=presentation(sk);red=eliminate(p);gs=red['generators'];rr=red['relators']
    rows=[[w.count(g)-w.count(-g) for g in gs] for w in rr]
    M=S.Matrix(rows) if rows else S.zeros(0,len(gs))
    ns=M.nullspace()
    if len(ns)!=1:raise ValueError(('H1 rank',len(gs),len(ns)))
    den=math.lcm(*(int(z.q) for z in ns[0]));v=[int(z*den) for z in ns[0]]
    g=math.gcd(*v);v=[z//g for z in v]
    vals=dict(zip(gs,v));residual=dict(vals)
    for e in reversed(red['eliminations']):
        vals[e['generator']]=sum((1 if z>0 else -1)*vals[abs(z)] for z in e['replacement'])
    co=[0]*len(sk['edge_endpoints'])
    for i,e in enumerate(p['generator_edge_ids']):co[e]=vals[i+1]
    if any(sum((1 if z>0 else -1)*co[abs(z)-1] for z in t) for t in sk['triangle_boundaries']):raise ValueError('not cocycle')
    return co,{'presentation':p,'reduction':red,'residual_exponent_rows':rows,'residual_abelian_values':residual,
               'rational_nullity':1,'primitive_integer_cocycle':True}

def cycle_value(sk,D,co,path):
    ei={tuple(e):i for i,e in enumerate(sk['edge_representatives'])};z=0
    for a,b in zip(path,path[1:]):
        r=min(tuple(sorted((a,b))),tuple(sorted((D[a],D[b]))))
        z+=(1 if (a,b)==r or (D[a],D[b])==r else -1)*co[ei[r]]
    return z

def solve_integer_system(L,b):
    n=len(L);a=[list(row)+list(b[i]) for i,row in enumerate(L)];old=1
    for k in range(n-1):
        pv=a[k][k]
        if not pv:raise ValueError('unexpected zero SPD pivot')
        for i in range(k+1,n):
            aik=a[i][k]
            for j in range(k+1,n+2):
                q,rem=divmod(a[i][j]*pv-aik*a[k][j],old)
                if rem:raise ValueError('nonexact Bareiss division')
                a[i][j]=q
            a[i][k]=0
        old=pv
    sol=[[Q(0),Q(0)] for _ in range(n)]
    for i in range(n-1,-1,-1):
        for k in range(2):sol[i][k]=(Q(a[i][n+k])-sum(a[i][j]*sol[j][k] for j in range(i+1,n)))/a[i][i]
    return sol

def powword(w,n):return free((tuple(w) if n>=0 else inv(w))*abs(n))
def build(cut):
    D=cut['deck'];T=skeleton(cut['T'],D);N=skeleton(cut['N'],D);E=skeleton(cut['E'],D)
    print('EXTRACTING',len(N['edge_endpoints']),len(E['edge_endpoints']),flush=True)
    cn,ncert=abelian_cocycle(N);ce,ecert=abelian_cocycle(E)
    pts=list(map(tuple,cut['points']));pi={p:i for i,p in enumerate(pts)}
    core=[pi['v',v] for v in cut['oriented_fixed_circle_source_vertices']]
    sn=cycle_value(N,D,cn,core)
    if abs(sn)!=1:raise ValueError('core not primitive')
    cn=[a*sn for a in cn];ncert['sign_normalization']=sn
    # Normal circle in link of the first source fixed edge; the quotient is a half traversal.
    edge=cut['puncture']['source_fixed_edge'];eq=cut['puncture']['equator_cycle']
    # These equator vertices are already the level-set midpoints based at one endpoint.
    normal=eq[:len(eq)//2+1]
    if D[normal[0]]!=normal[-1]:raise ValueError('normal path does not change sheet')
    if any(p not in {v for t in cut['T'] for v in t} for p in normal):raise ValueError('normal not in torus')
    se=cycle_value(E,D,ce,normal)
    if abs(se)!=1:raise ValueError(('meridian not primitive',se))
    ce=[a*se for a in ce];ecert['sign_normalization']=se
    ixN={tuple(e):i for i,e in enumerate(N['edge_representatives'])};ixE={tuple(e):i for i,e in enumerate(E['edge_representatives'])}
    u=[cn[ixN[tuple(e)]] for e in T['edge_representatives']];v=[ce[ixE[tuple(e)]] for e in T['edge_representatives']]
    vs=T['vertices'];vi={v:i for i,v in enumerate(vs)};nn=len(vs)
    L=[[0]*(nn-1) for _ in range(nn-1)];rhs=[[0,0] for _ in range(nn-1)]
    for k,(a,b) in enumerate(T['edge_endpoints']):
        i,j=vi[a],vi[b]
        for x,y,sgn in ((i,j,1),(j,i,-1)):
            if not x:continue
            L[x-1][x-1]+=1
            if y:L[x-1][y-1]-=1
            rhs[x-1][0]+=sgn*u[k];rhs[x-1][1]+=sgn*v[k]
    xy=dict(zip(vs,[[Q(0),Q(0)]]+solve_integer_system(L,rhs)))
    up=skeleton(cut['T'],D,False);di={tuple(e):i for i,e in enumerate(T['edge_representatives'])}
    edge_map=[];adj=defaultdict(list)
    for a,b in up['edge_representatives']:
        rr=min((a,b),tuple(sorted((D[a],D[b]))));i=di[rr];sg=1 if (a,b)==rr or (D[a],D[b])==rr else -1
        edge_map.append(sg*(i+1));adj[a].append((b,v[i]%2));adj[b].append((a,v[i]%2))
    root=min(adj);eps={root:0};todo=deque([root])
    while todo:
        a=todo.popleft()
        for b,k in adj[a]:
            z=(eps[a]+k)%2
            if b in eps:
                if eps[b]!=z:raise ValueError('meridional covering parity')
            else:eps[b]=z;todo.append(b)
    if not all(eps[D[a]]==1-eps[a] for a in eps):raise ValueError('deck parity')
    ux={a:[xy[min(a,D[a])][0],(xy[min(a,D[a])][1]+eps[a])/2] for a in up['vertices']}
    uc=[]
    for (a,b),ee in zip(up['edge_endpoints'],edge_map):
        i=abs(ee)-1;sg=1 if ee>0 else -1;z=sg*v[i]+eps[a]-eps[b]
        if z%2:raise ValueError('half-integral period')
        uc.append([sg*u[i],z//2])
    tp,paths,_=presentation(T)
    val=lambda w,co:sum((1 if e>0 else -1)*co[abs(e)-1] for e in w)
    cycles=[]
    for j,(a,b) in enumerate(T['edge_endpoints']):
        w=free(paths[a]+(j+1,)+inv(paths[b]));pr=(val(w,u),val(w,v))
        if pr!=(0,0):cycles.append((len(w),pr,w))
    cycles.sort();pair=None
    for x,y in combinations(cycles,2):
        a,b=x[1];c,d=y[1];det=a*d-b*c
        if abs(det)==1:pair=x,y,a,b,c,d,det;break
    if pair is None:raise ValueError('no primitive peripheral pair')
    x,y,a,b,c,d,det=pair
    longitude=free(powword(x[2],d//det)+powword(y[2],-b//det))
    meridian=free(powword(x[2],-c//det)+powword(y[2],a//det))
    ep,epaths,eword=presentation(E);eids={tuple(e):i+1 for i,e in enumerate(E['edge_representatives'])};emap=[eids[tuple(e)] for e in T['edge_representatives']]
    peripheral={}
    for name,w in [('longitude',longitude),('meridian',meridian)]:
        ew=tuple((1 if e>0 else -1)*emap[abs(e)-1] for e in w);ew=epaths[T['vertices'][0]]+ew+inv(epaths[T['vertices'][0]])
        peripheral[name]={'torus_edge_loop':list(w),'periods':[val(w,u),val(w,v)],
                          'raw_exterior_word':list(eword(ew)),'reduced_exterior_word':rewrite(eword(ew),ecert['reduction'])}
    return {'schema':'marked-pretzel-peripheral-chart-v1','cut_sha256':digest(cut),
            'N_skeleton':N,'E_skeleton':E,'torus_downstairs':T,'torus_upstairs':up,
            'N_H1_certificate':ncert,'E_H1_certificate':ecert,'N_cocycle':cn,'E_cocycle':ce,
            'oriented_core_point_cycle':core,'oriented_base_meridian_lift':normal,
            'downstairs_periods':list(map(list,zip(u,v))),
            'downstairs_vertex_coordinates':[[v,*map(str,xy[v])] for v in vs],
            'upstairs_vertex_coordinates':[[v,*map(str,ux[v])] for v in up['vertices']],
            'upstairs_periods':uc,'upstairs_edge_projection':edge_map,'lift_parity':[[a,eps[a]] for a in sorted(eps)],
            'peripheral_loops':peripheral,
            'coordinate_meaning':{'u':'zero-framed longitude; tube-degree +1, exterior-degree 0',
                                  'v_downstairs':'meridian positive by selected normal half-circle',
                                  'v_upstairs':'double meridian; projection (u,v)->(u,2v)'},
            'diagram_sign_transport':'PENDING explicit source-oriented comparison'}

if __name__=='__main__':
    st=time.monotonic();x=build(read(ROOT/'peripheral/cut_source.json.gz'))
    save(ROOT/'peripheral/peripheral_chart.json.gz',x)
    print('CHART',len(x['torus_downstairs']['vertices']),len(x['torus_downstairs']['triangle_representatives']),
          {k:v['periods'] for k,v in x['peripheral_loops'].items()},time.monotonic()-st,flush=True)

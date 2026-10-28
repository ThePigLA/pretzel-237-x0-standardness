#!/usr/bin/env python3
"""Source-independent even-twist cap meshes, with rational common refinement.

The boundary is a specified PL S^2 x S^1.  This module does not identify or
attach it to a knot source.  All angular/time charts retain universal lifts.
No inherited trefoil builder is imported or changed.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import argparse, gzip, hashlib, json, math


def det(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def area2(p):
    return sum(p[i][0]*p[(i+1)%len(p)][1]-p[i][1]*p[(i+1)%len(p)][0] for i in range(len(p)))


def clean(p):
    out=[]
    for x in p:
        if not out or x!=out[-1]: out.append(x)
    if len(out)>1 and out[0]==out[-1]:out.pop()
    return out


def clip(p,a):
    if area2(a)<0:a=list(reversed(a))
    for x,y in zip(a,a[1:]+a[:1]):
        if not p:return []
        out=[]
        for u,v in zip(p,p[1:]+p[:1]):
            s,t=det(x,y,u),det(x,y,v)
            if s>=0:out.append(u)
            if (s<0<t) or (t<0<s):
                r=s/(s-t);out.append(tuple(u[k]+r*(v[k]-u[k]) for k in (0,1)))
        p=clean(out)
    if len(p)<3 or not area2(p):return []
    if area2(p)<0:p.reverse()
    return p


def mod(p):return tuple(v-v.numerator//v.denominator for v in p)
def enc(p):return tuple(str(x) for x in p)
def E(p):return ('E',*enc(mod(p)))


def grid(n,t,shear=0):
    out=[]
    for i in range(n):
        for j in range(t):
            def p(a,b):return (Q(a,n)+shear*Q(b,t),Q(b,t))
            for s,ij in enumerate((((i,j),(i+1,j),(i+1,j+1)),((i,j),(i+1,j+1),(i,j+1)))):
                out.append({'i':i,'j':j,'s':s,'ids':tuple((a%n,b%t) for a,b in ij),'points':tuple(p(a,b) for a,b in ij)})
    return out


def common_torus(n,t):
    a,b=grid(n,t),grid(n,t,-1)
    boxes=lambda f:tuple((min(v[k] for v in f['points']),max(v[k] for v in f['points'])) for k in (0,1))
    bb=[boxes(f) for f in b];pieces=[]
    for ia,fa in enumerate(a):
        ab=boxes(fa)
        for ib,fb in enumerate(b):
            ranges=[range(math.ceil(ab[k][0]-bb[ib][k][1]),math.floor(ab[k][1]-bb[ib][k][0])+1) for k in (0,1)]
            for shift in product(*ranges):
                bp=[tuple(p[k]+shift[k] for k in (0,1)) for p in fb['points']]
                p=clip(list(fa['points']),bp)
                if p:pieces.append({'a':ia,'b':ib,'shift':shift,'polygon':p})
    # Retain every collinear boundary vertex so all polygon fans fit exactly.
    verts={mod(p) for row in pieces for p in row['polygon']}
    for row in pieces:
        old=row['polygon'];p=[]
        for u,v in zip(old,old[1:]+old[:1]):
            cand={u,v}
            for q in verts:
                for dx in range(math.floor(min(u[0],v[0]))-1,math.ceil(max(u[0],v[0]))+2):
                    for dt in range(math.floor(min(u[1],v[1]))-1,math.ceil(max(u[1],v[1]))+2):
                        z=(q[0]+dx,q[1]+dt)
                        if det(u,v,z)==0 and all(min(u[k],v[k])<=z[k]<=max(u[k],v[k]) for k in (0,1)):cand.add(z)
            k=0 if v[0]!=u[0] else 1
            order=sorted(cand,key=lambda z:(z[k]-u[k])/(v[k]-u[k]));p.extend(order[:-1])
        row['polygon']=p
        c=tuple(sum(z[k] for z in p)/len(p) for k in (0,1))
        row['centre']=c
        row['triangles']=[(c,u,v) for u,v in zip(p,p[1:]+p[:1])]
    return a,b,pieces


def stair(a,b):
    """All monotone lattice paths in Delta(a) x Delta(b)."""
    out=[]
    def rec(i,j,path):
        if i==len(a)-1 and j==len(b)-1:out.append(tuple(path));return
        if i<len(a)-1:rec(i+1,j,path+[(a[i+1],b[j])])
        if j<len(b)-1:rec(i,j+1,path+[(a[i],b[j+1])])
    rec(0,0,[(a[0],b[0])]);return out


def determinant3(a,b,c,d):
    x,y,z=[tuple(v[k]-a[k] for k in range(3)) for v in (b,c,d)]
    return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])


def shared_boundary(n,t,pieces,a):
    parts={p:[] for p in (0,1)}
    torus=[]
    for row in pieces:
        j=a[row['a']]['j']
        for tr in row['triangles']:
            et=tuple(E(x) for x in tr);torus.append(et)
            for p in (0,1):parts[p].append((('C',p,j),*et))
    for p in (0,1):
        for j in range(t):
            c=('C',p,j)
            for jj in (j,j+1):
                for i in range(n):parts[p].append((c,('P',p,jj%t),E((Q(i,n),Q(jj,t))),E((Q(i+1,n),Q(jj,t)))))
    return parts,torus


def source_children(m,pieces):
    children=defaultdict(list)
    k=m//2
    for row in pieces:
        parent=row['a'] if not k else row['b']
        shift=(0,0) if not k else row['shift']
        for tr in row['triangles']:
            pts=tuple((p[0]-shift[0]+k*(p[1]-shift[1]),p[1]-shift[1]) for p in tr)
            tags=tuple(E(p) for p in tr)
            children[parent].append((pts,tags))
    return children


def filling_torus(n,t,m,a,pieces):
    """S^1 x D^2 plus an explicit common-refinement collar."""
    hs=[];charts=[];children=source_children(m,pieces)
    for i in range(n):
        for j in range(t):
            def v(a,b):return ('H',a%n) if b==-1 else ('I',a%n,b%t)
            for tet in stair((i,i+1),(-1,j,j+1)):hs.append(tuple(v(x,y) for x,y in tet))
    for fi,f in enumerate(a):
        pp=f['points'];ids=f['ids'];co={('I',*idx):(*p,Q(0)) for idx,p in zip(ids,pp)}
        boundary=[tuple(('I',*idx) for idx in ids)]
        top_points={}
        for pts,tags in children[fi]:
            boundary.append(tags)
            for p,v in zip(pts,tags):co[v]=(*p,Q(1));top_points[p]=v
        for u,v,iu,iv in ((pp[0],pp[1],ids[0],ids[1]),(pp[1],pp[2],ids[1],ids[2]),(pp[2],pp[0],ids[2],ids[0])):
            cand=[p for p in top_points if det(u,v,p)==0 and all(min(u[k],v[k])<=p[k]<=max(u[k],v[k]) for k in (0,1))]
            k=0 if u[0]!=v[0] else 1;cand.sort(key=lambda p:(p[k]-u[k])/(v[k]-u[k]))
            if not cand or cand[0]!=u or cand[-1]!=v:raise ValueError(('edge endpoints',fi,u,v,cand))
            edge=tuple(sorted((iu,iv)));cen=('S',*edge)
            co[cen]=((u[0]+v[0])/2,(u[1]+v[1])/2,Q(1,2))
            lo_u,lo_v=('I',*iu),('I',*iv)
            edges=[(lo_u,lo_v),(lo_v,top_points[v]),(top_points[u],lo_u)]
            edges.extend((top_points[x],top_points[y]) for x,y in zip(cand,cand[1:]))
            boundary.extend((cen,*edge) for edge in edges)
        center=('F',fi);co[center]=tuple(sum(p[k] for p in pp)/3 for k in (0,1))+(Q(1,2),)
        volume=Q(0)
        for tr in boundary:
            tet=(*tr,center);coords=tuple(co[v] for v in tet);d=abs(determinant3(*coords))
            if not d:raise ValueError(('zero collar simplex',m,fi,tet))
            volume+=d;hs.append(tet);charts.append({'parent':fi,'vertices':tet,'coords':coords})
        if volume!=3*abs(det(*pp)):raise ValueError(('collar volume',m,fi,volume,3*abs(det(*pp))))
    return hs,charts,children


def deck_tag(v,n,t):
    if v[0]=='E':return E((Q(v[1])+Q(1,2),Q(v[2])))
    if v[0] in ('A','P','C'):return v
    if v[0] in ('I','H'):return (v[0],(v[1]+n//2)%n,*v[2:])
    if v[0]=='S':return ('S',*sorted((((i+n//2)%n,j) for i,j in v[1:])))
    if v[0]=='F':
        i,r=divmod(v[1],2*t);return ('F',((i+n//2)%n)*2*t+r)
    raise ValueError(v)


def faces(facets):
    return [set(x for f in facets for x in combinations(sorted(f),k+1)) for k in range(len(facets[0]))]


def fv(facets):return [len(s) for s in faces(facets)]


def orient_and_glue(facets):
    slots=defaultdict(list)
    for i,f in enumerate(facets):
        for j in range(len(f)):slots[f[:j]+f[j+1:]].append((i,j))
    if any(len(v)>2 for v in slots.values()):raise ValueError('ridge multiplicity')
    adj=defaultdict(list);glu=[];bd=[]
    for ridge,hits in sorted(slots.items()):
        if len(hits)==1:bd.append(ridge);continue
        (i,k),(j,l)=hits
        perm=[facets[j].index(v) if q!=k else l for q,v in enumerate(facets[i])]
        sign=(-1)**sum(perm[a]>perm[b] for a in range(len(perm)) for b in range(a+1,len(perm)))
        adj[i].append((j,-sign));adj[j].append((i,-sign));glu.append([i,k,j,perm])
    signs={0:1};todo=[0]
    for i in todo:
        for j,s in adj[i]:
            if j in signs:
                if signs[j]!=signs[i]*s:raise ValueError('nonorientable')
            else:signs[j]=signs[i]*s;todo.append(j)
    if len(signs)!=len(facets):raise ValueError('disconnected')
    return glu,bd,[signs[i] for i in range(len(facets))]


def encode(obj):
    if isinstance(obj,Q):return str(obj)
    if isinstance(obj,dict):return {str(k):encode(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):return [encode(x) for x in obj]
    return obj


def digest(x):return hashlib.sha256(json.dumps(encode(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def save(path,obj):
    raw=(json.dumps(encode(obj),sort_keys=True,separators=(',',':'))+'\n').encode()
    if str(path).endswith('.gz'):
        with open(path,'wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as z:z.write(raw)
    else:Path(path).write_bytes(raw)


def boundary_refinement(m,n,t,a,pieces,btags,bd_facets,coarse):
    """Barycentric subdivision map from the shared boundary to source B0.

    This is beta_m^{-1}, linear on each saved refined simplex.  Its inverse
    gives beta_m into the common boundary triangulation.
    """
    bi={v:i for i,v in enumerate(btags)}
    ct=list(map(tup_local,coarse['vertices']));ci={v:i for i,v in enumerate(ct)}
    cf=[tuple(f) for f in coarse['facets']];cfi={f:i for i,f in enumerate(cf)}
    rows={}
    for fi,children in source_children(m,pieces).items():
        f=a[fi];ps=f['points'];den=det(*ps)
        for p in (0,1):
            center=('C',p,f['j']);cone=tuple(sorted([ci[center]]+[ci[E(v)] for v in ps]));carrier=cfi[cone]
            for pts,tags in children:
                fd=tuple(sorted([bi[center]]+[bi[v] for v in tags]));weights={bi[center]:{ci[center]:Q(1)}}
                for q,v in zip(pts,tags):
                    bary=(det(q,ps[1],ps[2])/den,det(ps[0],q,ps[2])/den,det(ps[0],ps[1],q)/den)
                    weights[bi[v]]={ci[E(ps[z])]:bary[z] for z in range(3)}
                rows[fd]={'parent':carrier,'barycentric':[[weights[v].get(w,Q(0)) for w in cone] for v in fd]}
    def inverse_tag(v):
        if v[0]=='E':return E((Q(v[1])+Q(m,2)*Q(v[2]),Q(v[2])))
        return v
    for fd in bd_facets:
        if fd in rows:continue
        ss={v:ci[inverse_tag(btags[v])] for v in fd};cone=tuple(sorted(ss.values()))
        rows[fd]={'parent':cfi[cone],'barycentric':[[Q(ss[v]==w) for w in cone] for v in fd]}
    return [rows[f] for f in bd_facets]


def tup_local(x):return tuple(tup_local(v) for v in x) if isinstance(x,(list,tuple)) else x


def build(out,n=8,t=4):
    if n<4 or n%2 or t<3 or n%t:raise ValueError('angles even >=4, times >=3, and times divide angles required')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    a,b,pieces=common_torus(n,t);parts,torus=shared_boundary(n,t,pieces,a)
    btags=sorted({v for ps in parts.values() for f in ps for v in f},key=repr);bi={v:i for i,v in enumerate(btags)}
    bd_facets=sorted(tuple(sorted(bi[v] for v in f)) for ps in parts.values() for f in ps)
    if len(bd_facets)!=len(set(bd_facets)):raise ValueError('duplicate boundary tetrahedra')
    coarse_parts,_=shared_boundary(n,t,[{'a':i,'triangles':[f['points']]} for i,f in enumerate(a)],a)
    coarse_tags=sorted({v for ps in coarse_parts.values() for f in ps for v in f},key=repr);ci={v:i for i,v in enumerate(coarse_tags)}
    coarse={'vertices':coarse_tags,'facets':sorted(tuple(sorted(ci[v] for v in f)) for ps in coarse_parts.values() for f in ps),'deck':[ci[deck_tag(v,n,t)] for v in coarse_tags]}
    shared={'schema':'even-cap-shared-boundary-v1','angles':n,'times':t,'vertices':btags,'deck':[bi[deck_tag(v,n,t)] for v in btags],'facets':bd_facets,'hemispheres':{str(p):sorted(tuple(sorted(bi[v] for v in f)) for f in fs) for p,fs in parts.items()},'equator':sorted(tuple(sorted(bi[v] for v in f)) for f in torus),'coarse_source_boundary':coarse,'base_grid':a,'sheared_grid':b,'overlay':pieces,'f_vector':fv(bd_facets),'scope':'source-independent PL S2 x S1 common refinement; no knot source attached'}
    save(out/'shared_boundary.json.gz',shared)
    summary={'schema':'even-cap-generation-v1','shared_boundary_f_vector':shared['f_vector'],'shared_boundary_sha256':digest(shared),'cases':{},'scope':'CAP_TEMPLATE_ONLY; no J source or X2 constructed'}
    for m in (0,2):
        hs,charts,children=filling_torus(n,t,m,a,pieces)
        cap_parts={p:[(('A',p),*f) for f in parts[p]+hs] for p in (0,1)}
        tags=btags+sorted({v for ps in cap_parts.values() for f in ps for v in f}-set(btags),key=repr);vi={v:i for i,v in enumerate(tags)}
        fs=sorted(tuple(sorted(vi[v] for v in f)) for ps in cap_parts.values() for f in ps)
        if len(fs)!=len(set(fs)) or any(len(set(f))!=5 for f in fs):raise ValueError('cap facet multiplicity')
        D=[vi[deck_tag(v,n,t)] for v in tags];fset=set(fs)
        if any(tuple(sorted(D[v] for v in f)) not in fset for f in fs):raise ValueError('cap deck does not preserve facets')
        if any(any(v!=D[v] and D[v] in f for v in f) for f in fs):raise ValueError('deck inversion inside simplex')
        glu,bd,signs=orient_and_glue(fs)
        if bd!=bd_facets:raise ValueError('cap boundary differs from shared boundary')
        fixed={v for v,w in enumerate(D) if v==w}
        ram=sorted({tuple(v for v in f if v in fixed) for f in fs if sum(v in fixed for v in f)==3})
        if any(sum(v in fixed for v in f)>3 for f in fs):raise ValueError('fixed dimension')
        twin=sorted(tuple(sorted((vi[('A',p)],vi[('H',i)],vi[('H',(i+1)%n)]))) for p in (0,1) for i in range(n))
        case={'schema':'even-twist-cap-template-v1','twist':m,'angles':n,'times':t,'vertices':tags,'deck':D,'facets':fs,'boundary':bd,'boundary_vertex_ids':list(range(len(btags))),'boundary_identification':'identity on shared boundary vertex ids; source charts below retain beta_m inverse','boundary_to_source_refinement':boundary_refinement(m,n,t,a,pieces,btags,bd_facets,coarse),'gluings':glu,'orientation':signs,'hemisphere_cones':{str(p):sorted(tuple(sorted(vi[v] for v in f)) for f in ps) for p,ps in cap_parts.items()},'filling_torus':sorted(tuple(sorted(vi[v] for v in f)) for f in hs),'collar_charts':[{'parent':r['parent'],'vertices':[vi[v] for v in r['vertices']],'coords':r['coords']} for r in charts],'source_boundary_children':{str(f):[{'source_lifts':ps,'target_vertices':[vi[v] for v in vs]} for ps,vs in rows] for f,rows in children.items()},'ramification_triangles':ram,'twin_triangles':twin,'f_vector':fv(fs),'beta':{'formula':'(theta,t) -> (theta - m*t/2,t), coned hemisphere extension','inverse_formula':'(theta,t) -> (theta + m*t/2,t), inverse coned extension','unwrapped_displacement':[-m*n//2,t],'normalized_displacement':[-m//2,1],'downstairs_winding':-m,'bulk_deck_parity':m%2,'cap_full_turns':-m//2},'shared_boundary_sha256':digest(shared),'scope':'Source-independent S2 x D2 cap with explicit boundary attachment class. No global relative knot-source transfer or X2 admission.'}
        save(out/f'cap_m{m}.json.gz',case)
        summary['cases'][str(m)]={'f_vector':case['f_vector'],'pentachora':len(fs),'ramification_f_vector':fv(ram),'twin_f_vector':fv(twin),'sha256':digest(case)}
    save(out/'generation.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'data');p.add_argument('--angles',type=int,default=8);p.add_argument('--times',type=int,default=4);a=p.parse_args();build(a.output,a.angles,a.times)

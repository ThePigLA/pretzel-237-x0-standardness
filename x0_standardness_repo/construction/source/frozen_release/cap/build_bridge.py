#!/usr/bin/env python3
"""Exact product-boundary bridge: source suspension(C4) x C4 to cap boundary.

An exact rational tetrahedron intersection complex gives a common refinement.
Two coned prism collars connect the coarse boundaries to that refinement.
The result is a source-bound S2 x S1 x I, not a closed X2 or X0.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse,gzip,hashlib,json,math,time
from build_caps import stair,faces,fv,orient_and_glue,encode,digest,save


def load(p):
    with gzip.open(p,'rt') if str(p).endswith('.gz') else open(p) as f:return json.load(f)
def tup(x):return tuple(tup(v) for v in x) if isinstance(x,list) else x
def mod(x):return x-x.numerator//x.denominator
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def det3(a,b,c,d):return dot(sub(b,a),cross(sub(c,a),sub(d,a)))


def detmatrix(a):
    a=[list(r) for r in a];out=Q(1)
    for j in range(len(a)):
        p=next((i for i in range(j,len(a)) if a[i][j]),None)
        if p is None:return Q(0)
        if p!=j:a[j],a[p]=a[p],a[j];out=-out
        q=a[j][j];out*=q
        for k in range(j,len(a)):a[j][k]/=q
        for i in range(j+1,len(a)):
            q=a[i][j]
            for k in range(j,len(a)):a[i][k]-=q*a[j][k]
    return out


def bary(p,tet):
    den=det3(*tet)
    if not den:raise ValueError('zero parent tetrahedron')
    return tuple(det3(*(list(tet[:j])+[p]+list(tet[j+1:])))/den for j in range(4))


def planes(tet):
    out=[]
    for j in range(4):
        f=[tet[i] for i in range(4) if i!=j];n=cross(sub(f[1],f[0]),sub(f[2],f[0]));d=-dot(n,f[0])
        if dot(n,tet[j])+d<0:n=tuple(-x for x in n);d=-d
        out.append((*n,d))
    return out


def value(plane,p):return dot(plane[:3],p)+plane[3]
def inside(p,ps):return all(value(h,p)>=0 for h in ps)
def box(t):return tuple((min(p[k] for p in t),max(p[k] for p in t)) for k in range(3))
def overlap(a,b):return all(max(a[k][0],b[k][0])<min(a[k][1],b[k][1]) for k in range(3))


def intersection(a,b,pa,pb):
    verts=set(p for p in a if inside(p,pb))|set(p for p in b if inside(p,pa))
    for tet,ps in ((a,pb),(b,pa)):
        for u,v in combinations(tet,2):
            for h in ps:
                x,y=value(h,u),value(h,v)
                if x*y>=0:continue
                q=x/(x-y);p=tuple(u[k]+q*(v[k]-u[k]) for k in range(3))
                if inside(p,pa) and inside(p,pb):verts.add(p)
    if len(verts)<4:return None
    verts=sorted(verts)
    if not any(det3(verts[0],*s) for s in combinations(verts[1:],3)):return None
    fset={}
    for h in pa+pb:
        f=frozenset(p for p in verts if value(h,p)==0)
        if len(f)>=3 and any(cross(sub(x,tuple(f)[0]),sub(y,tuple(f)[0]))!=(0,0,0) for x,y in combinations(tuple(f)[1:],2)):
            fset[f]=h
    return verts,list(fset.items())


def geomtag(p,hem):
    x,y,t=p
    return (hem if abs(x)+abs(y)!=1 else -1,str(x),str(y),str(mod(t)))


def deckgeom(v):return (v[0],str(-Q(v[1])),str(-Q(v[2])),v[3])


def orbitkey(p,hem):
    tag=geomtag(p,hem);x,y=Q(tag[1]),Q(tag[2]);a=min((x,y),(-x,-y));return (tag[0],*a,Q(tag[3]))


def polygon_order(vertices,plane):
    drop=max(range(3),key=lambda k:abs(plane[k]));axes=[k for k in range(3) if k!=drop]
    def xy(p):return tuple(p[k] for k in axes)
    def turn(a,b,c):
        a,b,c=xy(a),xy(b),xy(c);return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    ps=sorted(vertices,key=xy);lo=[];hi=[]
    for p in ps:
        while len(lo)>1 and turn(lo[-2],lo[-1],p)<=0:lo.pop()
        lo.append(p)
    for p in reversed(ps):
        while len(hi)>1 and turn(hi[-2],hi[-1],p)<=0:hi.pop()
        hi.append(p)
    polygon=lo[:-1]+hi[:-1]
    if len(polygon)!=len(vertices):raise ValueError('non-extreme point in intersection face')
    return polygon


def pulling(poly,hem):
    verts,facets=poly;root=min(verts,key=lambda p:orbitkey(p,hem));out=[]
    keys=[orbitkey(p,hem) for p in verts]
    if len(keys)!=len(set(keys)):raise ValueError('deck orbit repeats in intersection cell')
    for vs,h in facets:
        if root in vs:continue
        cyc=polygon_order(vs,h);k=min(range(len(cyc)),key=lambda i:orbitkey(cyc[i],hem));cyc=cyc[k:]+cyc[:k]
        for i in range(1,len(cyc)-1):
            te=(root,cyc[0],cyc[i],cyc[i+1])
            if not det3(*te):raise ValueError('zero pulling tetrahedron')
            out.append(te)
    return out


def circlepoint(theta):
    q=mod(theta)*4;i=q.numerator//q.denominator;s=q-i;v=[(Q(1),Q(0)),(Q(0),Q(1)),(Q(-1),Q(0)),(Q(0),Q(-1))]
    return tuple((1-s)*v[i][k]+s*v[(i+1)%4][k] for k in (0,1))


def source_boundary(interface):
    poles=interface['poles'];eq=interface['equator_cycle_diagram_positive'];T=interface['time_intervals'];coords={v:circlepoint(Q(i,4)) for i,v in enumerate(eq)}
    coords.update({p:(Q(0),Q(0)) for p in poles});rows=[];tags=sorted((v,j) for v in coords for j in range(T));vi={v:i for i,v in enumerate(tags)}
    for tri in interface['ordered_triangles']:
        pp=[p for p in poles if p in tri]
        if len(pp)!=1:raise ValueError('source triangle has wrong poles')
        hem=poles.index(pp[0])
        for j in range(T):
            for path in stair(tri,(j,j+1)):
                tags4=tuple((v,k%T) for v,k in path);pts=tuple((*coords[v],Q(k,T)) for v,k in path)
                rows.append({'vertices':tuple(vi[v] for v in tags4),'points':pts,'hemisphere':hem,'slab':j})
    fs=sorted(tuple(sorted(r['vertices'])) for r in rows);idx={f:i for i,f in enumerate(fs)}
    if len(fs)!=len(set(fs)):raise ValueError('source boundary duplicate')
    for r in rows:
        r['index']=idx[tuple(sorted(r['vertices']))]
        ordered=sorted(zip(r['vertices'],r['points']));r['vertices']=tuple(x[0] for x in ordered);r['points']=tuple(x[1] for x in ordered)
    D=[vi[(interface['deck'][str(v)],j)] for v,j in tags]
    return {'vertices':tags,'facets':fs,'deck':D,'charts':rows}


def cap_charts(b):
    tags=list(map(tup,b['vertices']));T=b['times'];rows=[]
    for i,f in enumerate(b['facets']):
        cs=[tags[v] for v in f if tags[v][0]=='C']
        if len(cs)!=1:raise ValueError('cap boundary hemisphere slab')
        _,hem,j=cs[0];pts=[]
        for v in f:
            tag=tags[v]
            if tag[0]=='E':
                x,y=circlepoint(Q(tag[1]));t=Q(tag[2]);t=t+1 if j==T-1 and t==0 else t;pts.append((x,y,t))
            elif tag[0]=='P':
                t=Q(tag[2],T);t=t+1 if j==T-1 and t==0 else t;pts.append((Q(0),Q(0),t))
            elif tag[0]=='C':pts.append((Q(0),Q(0),Q(2*j+1,2*T)))
            else:raise ValueError(tag)
        rows.append({'vertices':tuple(f),'points':tuple(pts),'hemisphere':hem,'slab':j,'index':i})
    return rows


def common(source,b):
    a=source['charts'];c=cap_charts(b);groups=defaultdict(list)
    for r in a:groups[(r['hemisphere'],r['slab'])].append(r)
    cachea={r['index']:(planes(r['points']),box(r['points'])) for r in a};rows=[];counts=Counter();start=time.monotonic()
    for ix,r in enumerate(c):
        pc,bc=planes(r['points']),box(r['points']);hem=r['hemisphere']
        for s in groups[(hem,r['slab'])]:
            pa,ba=cachea[s['index']]
            if not overlap(ba,bc):continue
            p=intersection(s['points'],r['points'],pa,pc)
            if not p:continue
            for te in pulling(p,hem):
                tags=tuple(geomtag(q,hem) for q in te)
                rows.append({'vertices':tags,'points':te,'source_parent':s['index'],'cap_parent':r['index'],'source_barycentric':[bary(q,s['points']) for q in te],'cap_barycentric':[bary(q,r['points']) for q in te]})
                counts[(hem,r['slab'])]+=1
        if ix%224==223:print('overlay',ix+1,'/',len(c),'common_tets',len(rows),'seconds',round(time.monotonic()-start,2),flush=True)
    tags=sorted({v for r in rows for v in r['vertices']});vi={v:i for i,v in enumerate(tags)}
    out=[]
    for r in rows:
        order=sorted(range(4),key=lambda j:vi[r['vertices'][j]])
        out.append({'vertices':tuple(vi[r['vertices'][j]] for j in order),'points':tuple(r['points'][j] for j in order),'source_parent':r['source_parent'],'cap_parent':r['cap_parent'],'source_barycentric':[r['source_barycentric'][j] for j in order],'cap_barycentric':[r['cap_barycentric'][j] for j in order]})
    out.sort(key=lambda r:r['vertices']);fs=[r['vertices'] for r in out]
    if len(fs)!=len(set(fs)):raise ValueError('common tetrahedron duplicate')
    D=[vi[deckgeom(v)] for v in tags];fset=set(fs)
    if any(tuple(sorted(D[v] for v in f)) not in fset for f in fs):raise ValueError('common deck action')
    glu,bd,ori=orient_and_glue(fs)
    if bd:raise ValueError('common boundary not closed')
    for name,coarse in (('source',source),('cap',b)):
        vols=Counter();images={}
        for r in out:
            parent=coarse['facets'][r[name+'_parent']];mat=r[name+'_barycentric'];vols[r[name+'_parent']]+=abs(detmatrix(mat))
            for v,weights in zip(r['vertices'],mat):
                image=tuple((w,x) for w,x in zip(parent,weights) if x)
                if v in images and images[v]!=image:raise ValueError(('common image mismatch',name,v,images[v],image))
                images[v]=image
        if len(vols)!=len(coarse['facets']) or any(x!=1 for x in vols.values()):raise ValueError(('common volumes',name,len(vols),[(k,v) for k,v in vols.items() if v!=1][:5]))
    return {'vertices':tags,'facets':fs,'deck':D,'refinements':out,'orientation':ori,'f_vector':fv(fs)}


def prism_collar(side,coarse,fine):
    ff=faces(fine['facets']);cf=faces(coarse['facets']);images={};parents=side+'_parent';barys=side+'_barycentric'
    for r in fine['refinements']:
        f=coarse['facets'][r[parents]]
        for v,baryv in zip(r['vertices'],r[barys]):images[v]={w:x for w,x in zip(f,baryv) if x}
    top=[defaultdict(list) for _ in range(4)]
    for k in range(4):
        for f in ff[k]:
            carrier=tuple(sorted({w for v in f for w in images[v]}))
            if len(carrier)==k+1:top[k][carrier].append(tuple(('R',v) for v in f))
    prisms={};out=[];charts=[]
    for k in range(4):
        for f in sorted(cf[k]):
            if k==0:
                tops=top[k][f]
                if len(tops)!=1:raise ValueError(('coarse vertex preimage',side,f,tops))
                prisms[f]=[(('B',side,f[0]),tops[0][0])];continue
            # Pull from the least bottom vertex.  Only the top and the
            # opposite face prism are visible.  A global deck-orbit order
            # makes these triangulations restrict identically to shared
            # faces; no symmetry-breaking center choices are introduced.
            ranks=[min(v,coarse['deck'][v]) for v in f]
            if len(ranks)!=len(set(ranks)):raise ValueError('coarse face repeats a deck orbit')
            pivot=min(f,key=lambda v:min(v,coarse['deck'][v]));opposite=tuple(v for v in f if v!=pivot)
            boundaries=top[k][f]+prisms[opposite]
            root=('B',side,pivot);prisms[f]=[(*g,root) for g in boundaries]
            if k==3:
                parent_index={tuple(x):i for i,x in enumerate(coarse['facets'])}[f]
                def coords(v):
                    if v[0]=='B':im={v[2]:Q(1)};height=Q(0)
                    elif v[0]=='R':im=images[v[1]];height=Q(1)
                    else:im={w:Q(1,len(v[2])) for w in v[2]};height=Q(1,2)
                    return tuple(im.get(w,Q(0)) for w in f)+(height,)
                volume=Q(0)
                for pent in prisms[f]:
                    cs=[coords(v) for v in pent];mat=[tuple(cs[i][k]-cs[0][k] for k in (1,2,3,4)) for i in range(1,5)];d=abs(detmatrix(mat))
                    if not d:raise ValueError('degenerate collar pentachoron')
                    volume+=d;out.append(pent);charts.append({'parent':parent_index,'vertices':pent,'barycentric_height':cs})
                if volume!=4:raise ValueError(('collar volume',side,f,volume))
    return out,charts


def build(interface_path,shared_path,out):
    start=time.monotonic();interface=load(interface_path);b=load(shared_path);source=source_boundary(interface)
    if interface['time_intervals']!=b['times']:raise ValueError('time intervals mismatch')
    fine=common(source,b);print('COMMON',fine['f_vector'],flush=True)
    parts={};chs={}
    for side,coarse in (('source',source),('cap',b)):
        parts[side],chs[side]=prism_collar(side,coarse,fine);print('COLLAR',side,len(parts[side]),flush=True)
    if sum(map(len,parts.values()))>60000:raise RuntimeError('bridge preallocation cap: more than 60,000 pentachora')
    tags=sorted({v for fs in parts.values() for f in fs for v in f},key=repr);vi={v:i for i,v in enumerate(tags)}
    fs=sorted(tuple(sorted(vi[v] for v in f)) for ps in parts.values() for f in ps)
    def Dtag(v):
        if v[0]=='R':return ('R',fine['deck'][v[1]])
        co=source if v[1]=='source' else b
        if v[0]=='B':return ('B',v[1],co['deck'][v[2]])
        return ('M',v[1],tuple(sorted(co['deck'][i] for i in v[2])))
    D=[vi[Dtag(v)] for v in tags];fset=set(fs)
    if len(fs)!=len(fset):raise ValueError('bridge duplicate pentachora')
    if any(tuple(sorted(D[v] for v in f)) not in fset for f in fs):raise ValueError('bridge deck')
    glu,bd,orientation=orient_and_glue(fs)
    ports={side:{'vertices':[vi[('B',side,i)] for i in range(len(co['vertices']))],'facets':sorted(tuple(sorted(vi[('B',side,v)] for v in f)) for f in co['facets'])} for side,co in (('source',source),('cap',b))}
    if sorted(ports['source']['facets']+ports['cap']['facets'])!=bd:raise ValueError('bridge boundary ports')
    fixed={v for v,w in enumerate(D) if v==w};ram=sorted({tuple(v for v in f if v in fixed) for f in fs if sum(v in fixed for v in f)==3})
    if any(sum(v in fixed for v in f)>3 for f in fs):raise ValueError('bridge fixed locus dimension')
    source_cut_path=Path(interface_path).parent/'cut_source.json.gz'
    cut_binding={'file_sha256':hashlib.sha256(source_cut_path.read_bytes()).hexdigest(),'puncture_canonical_sha256':digest(load(source_cut_path)['puncture'])} if source_cut_path.exists() else None
    result={'schema':'source-cap-boundary-bridge-v1','scope':'Certified-boundary refinement collar only; no closed X2 or X0 constructed','source_interface':interface,'source_interface_file_sha256':hashlib.sha256(Path(interface_path).read_bytes()).hexdigest(),'source_cut_binding':cut_binding,'cap_shared_canonical_sha256':digest(b),'source_boundary':source,'common_boundary':fine,'vertices':tags,'facets':fs,'deck':D,'ports':ports,'gluings':glu,'orientation':orientation,'parts':{s:sorted(tuple(sorted(vi[v] for v in f)) for f in ps) for s,ps in parts.items()},'collar_charts':{s:[{'parent':r['parent'],'vertices':[vi[v] for v in r['vertices']],'barycentric_height':r['barycentric_height']} for r in rows] for s,rows in chs.items()},'ramification_triangles':ram,'f_vector':fv(fs),'elapsed_seconds':time.monotonic()-start}
    out=Path(out);out.mkdir(parents=True,exist_ok=True);save(out/'source_cap_bridge.json.gz',result)
    summary={'schema':'source-cap-bridge-generation-v1','pentachora':len(fs),'f_vector':result['f_vector'],'common_boundary_f_vector':fine['f_vector'],'source_boundary_tetrahedra':len(source['facets']),'cap_boundary_tetrahedra':len(b['facets']),'parts':{s:len(p) for s,p in parts.items()},'ramification_f_vector':fv(ram),'seconds':result['elapsed_seconds'],'canonical_sha256':digest(result)};save(out/'bridge_generation.json',summary);print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    root=Path(__file__).resolve().parent;p=argparse.ArgumentParser();p.add_argument('--interface',type=Path,default=root.parent/'peripheral/puncture_product_interface.json');p.add_argument('--shared',type=Path,default=root/'data/shared_boundary.json.gz');p.add_argument('--output',type=Path,default=root/'data');a=p.parse_args();build(a.interface,a.shared,a.output)

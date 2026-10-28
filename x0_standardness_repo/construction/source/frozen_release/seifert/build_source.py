#!/usr/bin/env python3
"""Explicit equivariant Montesinos source for pretzel parameters (-2,3,7).

Build a reflected pair of pants times a circle, then attach three reflected
solid tori by A_p=[[p,1],[1,0]]. All boundary common refinements use exact
rational coordinates. This script makes geometry; it does not identify a
manifold from a census name or substitute an expected homology calculation.
"""
from __future__ import annotations
from fractions import Fraction as Q
from itertools import product, combinations
from collections import defaultdict, Counter
from pathlib import Path
import argparse, gzip, hashlib, json, math, time

N = 4

def ordered(x):
    return tuple(sorted(x, key=repr))

def parity(seq):
    return (-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))

def order_sign(seq):
    target=ordered(seq)
    return parity([target.index(v) for v in seq])

def mod(x):
    return x - x.numerator // x.denominator

def pointkey(i, p):
    return ('T', i, str(mod(p[0])), str(mod(p[1])))

def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]

def sub(a,b):
    return (a[0]-b[0],a[1]-b[1])

def det(a,b,c):
    return cross(sub(b,a),sub(c,a))

def area(P):
    return sum(cross(a,b) for a,b in zip(P,P[1:]+P[:1]))

def clean(P):
    out=[]
    for p in P:
        if not out or p!=out[-1]:out.append(p)
    if len(out)>1 and out[0]==out[-1]:out.pop()
    return out

def clip(P,A):
    P=list(P);A=list(A)
    if area(A)<0:A.reverse()
    for a,b in zip(A,A[1:]+A[:1]):
        if not P:return []
        out=[]
        for c,d in zip(P,P[1:]+P[:1]):
            x,y=det(a,b,c),det(a,b,d)
            if x>=0:out.append(c)
            if (x<0<y) or (y<0<x):
                t=x/(x-y)
                out.append(tuple(c[k]+t*(d[k]-c[k]) for k in range(2)))
        P=clean(out)
    if len(P)<3 or area(P)==0:return []
    if area(P)<0:P.reverse()
    return P

def avg(P):
    return tuple(sum(p[k] for p in P)/len(P) for k in range(2))

def save(path,data):
    data=json.dumps(data,sort_keys=True,separators=(',',':')).encode()+b'\n'
    if path.suffix=='.gz':data=gzip.compress(data,mtime=0)
    path.write_bytes(data)

def digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()


class Block:
    def __init__(self,name,triangles,reflection,cycles):
        self.name=name
        self.triangles=[tuple(t) for t in triangles]
        self.reflection=reflection
        self.cycles=cycles  # (interface, boundary cycle, target matrix)
        self.vertices=set(reflection)
        self.edges={ordered(e) for t in triangles for e in combinations(t,2)}
        self.boundary_location={}
        self.boundary_edges={}
        for interface,cyc,matrix in cycles:
            for k,v in enumerate(cyc):
                self.boundary_location[v]=(interface,Q(k,len(cyc)),matrix)
                self.boundary_edges[ordered((v,cyc[(k+1)%len(cyc)]))]=(interface,k,cyc,matrix)

    def v(self,v,j):
        j%=N
        if v in self.boundary_location:
            i,x,A=self.boundary_location[v];y=Q(j,N)
            return pointkey(i,(A[0][0]*x+A[0][1]*y,A[1][0]*x+A[1][1]*y))
        return ('V',self.name,v,j)

    def edge(self,a,b):
        v,j=a;w,k=b
        if j==k:return ('E',self.name,'s',ordered((v,w)),j)
        if v!=w:raise ValueError('not an original product edge')
        if (j+1)%N==k:return ('E',self.name,'c',v,j)
        if (k+1)%N==j:return ('E',self.name,'c',v,k)
        raise ValueError('not a circle edge')

    def endpoints(self,e):
        if e[2]=='s':return (self.v(e[3][0],e[4]),self.v(e[3][1],e[4]))
        return (self.v(e[3],e[4]),self.v(e[3],(e[4]+1)%N))

    def face(self,t,j,kind):
        return ('F',self.name,kind,ordered(t),j)

    def faces_of_prism(self,t,j):
        return [self.face(t,j,'s'),self.face(t,(j+1)%N,'s')]+[self.face(e,j,'q') for e in combinations(t,2)]

    def face_loop(self,f):
        _,_,kind,t,j=f
        if kind=='s':return [(v,j) for v in t]
        a,b=t
        return [(a,j),(b,j),(b,(j+1)%N),(a,(j+1)%N)]

    def deck(self,k):
        if k[0]=='V':return ('V',self.name,self.reflection[k[2]],(-k[3])%N)
        if k[0]=='F':
            _,_,kind,t,j=k
            return self.face(tuple(self.reflection[v] for v in t),(-j if kind=='s' else -j-1)%N,kind)
        if k[0]=='C':
            _,_,t,j=k
            return ('C',self.name,ordered(self.reflection[v] for v in t),(-j-1)%N)
        raise ValueError(k)

    def boundary_records(self):
        records=[]
        for i,cyc,A in self.cycles:
            for k,v in enumerate(cyc):
                w=cyc[(k+1)%len(cyc)]
                for j in range(N):
                    xy=[(Q(k,len(cyc)),Q(j,N)),(Q(k+1,len(cyc)),Q(j,N)),
                        (Q(k+1,len(cyc)),Q(j+1,N)),(Q(k,len(cyc)),Q(j+1,N))]
                    ps=[(A[0][0]*x+A[0][1]*y,A[1][0]*x+A[1][1]*y) for x,y in xy]
                    loop=[(v,j),(w,j),(w,(j+1)%N),(v,(j+1)%N)]
                    records.append({'interface':i,'face':self.face((v,w),j,'q'),'points':ps,'loop':loop})
        return records


def make_blocks(params):
    # Double the hexagon on e01,e23,e45. Other edges are independently
    # subdivided in the two halves, giving four edges on every boundary circle.
    tris=[];R={str(i):str(i) for i in range(6)}
    for sign in ('+','-'):
        cen='c'+sign;R[cen]='c'+('-' if sign=='+' else '+')
        cyc=[]
        for k in range(3):
            b='b'+str(k)+sign;R[b]='b'+str(k)+('-' if sign=='+' else '+')
            cyc.extend([str(2*k),str(2*k+1),b])
        tris.extend((cen,cyc[j],cyc[(j+1)%len(cyc)]) for j in range(len(cyc)))
    cycles=[(k,[str(2*k+1),'b'+str(k)+'+',str((2*k+2)%6),'b'+str(k)+'-'],((1,0),(0,1))) for k in range(3)]
    blocks=[Block('pants',tris,R,cycles)]
    for i,p in enumerate(params):
        R={'o':'o',**{str(j):str((-j)%4) for j in range(4)}}
        tri=[('o',str(j),str((j+1)%4)) for j in range(4)]
        blocks.append(Block('fill'+str(i),tri,R,[(i,[str(j) for j in range(4)],((p,1),(1,0)))]))
    return blocks


def build(params=(-2,3,7)):
    blocks=make_blocks(params);byname={b.name:b for b in blocks}
    edge_points=defaultdict(dict);face_triangles=defaultdict(list);face_orientations=defaultdict(dict)
    overlay_records=[]

    def mark_edge(block,rec,a,b,p,key):
        pa,pb=rec['points'][a],rec['points'][b]
        d=sub(pb,pa)
        if det(pa,pb,p)!=0:return
        z=0 if d[0] else 1;t=(p[z]-pa[z])/d[z]
        if not 0<=t<=1:return
        va,vb=rec['loop'][a],rec['loop'][b];e=block.edge(va,vb)
        start,end=block.endpoints(e)
        if block.v(*va)!=start:t=1-t
        old=edge_points[e].get(t)
        if old is not None and old!=key:raise ValueError('inconsistent edge point')
        edge_points[e][t]=key

    total_polygons=0
    for i,p in enumerate(params):
        A=[r for r in blocks[0].boundary_records() if r['interface']==i]
        B=blocks[i+1].boundary_records();aa=[Q(0)]*len(A);bb=[Q(0)]*len(B)
        for ai,a in enumerate(A):
            ap=a['points'];boundsA=[(min(x[k] for x in ap),max(x[k] for x in ap)) for k in range(2)]
            for bi,b in enumerate(B):
                bp=b['points'];boundsB=[(min(x[k] for x in bp),max(x[k] for x in bp)) for k in range(2)]
                shifts=[range(math.ceil(boundsA[k][0]-boundsB[k][1]),math.floor(boundsA[k][1]-boundsB[k][0])+1) for k in range(2)]
                for shift in product(*shifts):
                    P=clip(ap,[(v[0]+shift[0],v[1]+shift[1]) for v in bp])
                    if not P:continue
                    total_polygons+=1;aa[ai]+=area(P);bb[bi]+=area(P)
                    centre=pointkey(i,avg(P));tris=[];points=[];natural=[]
                    for u,v in zip(P,P[1:]+P[:1]):
                        mid=avg((u,v));uk,vk,mk=pointkey(i,u),pointkey(i,v),pointkey(i,mid)
                        natural.extend([(uk,mk,centre),(mk,vk,centre)])
                        tris.extend([ordered((uk,mk,centre)),ordered((mk,vk,centre))]);points.extend([u,mid])
                    if any(len(set(t))!=3 for t in tris):raise ValueError('degenerate overlay triangle')
                    face_triangles[a['face']].extend(tris);face_triangles[b['face']].extend(tris)
                    for rec in (a,b):
                        canonical_sign=(1 if area(rec['points'])>0 else -1)*(1 if rec['loop'][0][0]==rec['face'][3][0] else -1)
                        for nt in natural:face_orientations[rec['face']][ordered(nt)]=canonical_sign*order_sign(nt)
                    for q in points:
                        key=pointkey(i,q)
                        for k in range(4):
                            mark_edge(blocks[0],a,k,(k+1)%4,q,key)
                            mark_edge(blocks[i+1],b,k,(k+1)%4,(q[0]-shift[0],q[1]-shift[1]),key)
                    overlay_records.append({'interface':i,'central_face':a['face'],'filling_face':b['face'],
                        'translation':shift,'polygon':[[str(z) for z in q] for q in P],
                        'triangles':tris})
        if any(x!=abs(area(r['points'])) for x,r in zip(aa,A)):raise ValueError('central overlay area')
        if any(x!=abs(area(r['points'])) for x,r in zip(bb,B)):raise ValueError('filling overlay area')

    prisms=[];allfaces={}
    for block in blocks:
        for tri in block.triangles:
            for j in range(N):
                fs=block.faces_of_prism(tri,j)
                prisms.append((block,('C',block.name,ordered(tri),j),fs))
                for f in fs:allfaces[f]=block
    for f,block in allfaces.items():
        if f in face_triangles:continue
        loop=block.face_loop(f);tris=[]
        for a,b in zip(loop,loop[1:]+loop[:1]):
            e=block.edge(a,b);start,end=block.endpoints(e)
            pts=dict(edge_points[e]);pts[Q(0)]=start;pts[Q(1)]=end
            seq=[key for _,key in sorted(pts.items())]
            if block.v(*a)!=start:seq.reverse()
            for u,v in zip(seq,seq[1:]):
                nt=(f,u,v);tri=ordered(nt);tris.append(tri);face_orientations[f][tri]=order_sign(nt)
        face_triangles[f]=tris
    tets=[];parts=defaultdict(list);tet_orientation={}
    for block,c,fs in prisms:
        # Original oriented product cell: reflected pants half is negative;
        # each other surface triangle uses its displayed vertex order.
        tr=next(t for t in block.triangles if ordered(t)==c[2]);j=c[3]
        sigma=-1 if block.name=='pants' and tr[0]=='c-' else 1
        signs={block.face(tr,j,'s'):-sigma*order_sign(tr),block.face(tr,(j+1)%N,'s'):sigma*order_sign(tr)}
        for a,b in ((tr[1],tr[2]),(tr[2],tr[0]),(tr[0],tr[1])):
            signs[block.face((a,b),j,'q')]=sigma*order_sign((a,b))
        for f in fs:
            for tri in face_triangles[f]:
                tet=ordered((c,*tri))
                if len(set(tet))!=4:raise ValueError('degenerate tetrahedron')
                tets.append(tet);parts[block.name].append(tet)
                tet_orientation[tet]=signs[f]*face_orientations[f][tri]*order_sign((c,*tri))
    if len(set(tets))!=len(tets):raise ValueError('duplicate tetrahedron')
    verts=sorted({v for t in tets for v in t},key=repr);ix={v:i for i,v in enumerate(verts)}
    def D(k):
        if k[0]=='T':return pointkey(k[1],(-Q(k[2]),-Q(k[3])))
        return byname[k[1]].deck(k)
    missing=[(v,D(v)) for v in verts if D(v) not in ix]
    if missing:raise ValueError(('deck vertex missing',missing[:5]))
    deck=[ix[D(v)] for v in verts]
    facets=sorted(tuple(sorted(ix[v] for v in t)) for t in tets);fs=set(facets)
    if any(tuple(sorted(deck[v] for v in t)) not in fs for t in fs):raise ValueError('deck tetrahedron missing')
    fixed={i for i,d in enumerate(deck) if i==d}
    if any(len({min(v,deck[v]) for v in t})<4 for t in facets):raise ValueError('folded tetrahedron')
    fixededges=sorted({e for t in facets for e in combinations(t,2) if set(e)<=fixed})
    if any(len(set(t)&fixed)>2 for t in facets):raise ValueError('higher fixed locus')
    fc=Counter(v for e in fixededges for v in e)
    if set(fc)!=fixed or any(x!=2 for x in fc.values()):raise ValueError('fixed locus not circles')
    fixed_components=[];todo=set(fixed)
    adj=defaultdict(set)
    for a,b in fixededges:adj[a].add(b);adj[b].add(a)
    while todo:
        seen={min(todo)};walk=list(seen)
        for v in walk:
            for w in adj[v]:
                if w not in seen:seen.add(w);walk.append(w)
        todo-=seen;fixed_components.append(sorted(seen))
    tri_slots=defaultdict(list)
    for i,t in enumerate(facets):
        for k in range(4):tri_slots[t[:k]+t[k+1:]].append((i,k))
    bad=[(f,h) for f,h in tri_slots.items() if len(h)!=2]
    if bad:raise ValueError(('triangle incidence',len(bad),bad[:3]))
    orientations=[tet_orientation[tuple(verts[v] for v in t)] for t in facets]
    gluings=[];neighbor={}
    for face,hits in sorted(tri_slots.items()):
        (i,k),(j,l)=hits
        perm=[facets[j].index(v) if a!=k else l for a,v in enumerate(facets[i])]
        gluings.append([i,k,j,perm]);neighbor[i,k]=(j,l)
        neighbor[j,l]=(i,k)
        if orientations[j]!=-parity(perm)*orientations[i]:raise ValueError(('source orientation mismatch',i,j))
    ti={t:i for i,t in enumerate(facets)}
    td=[ti[tuple(sorted(deck[v] for v in t))] for t in facets]
    if any(i==j for i,j in enumerate(td)):raise ValueError('fixed tetrahedron')
    for i,t in enumerate(facets):
        image=facets[td[i]];perm=[image.index(deck[v]) for v in t]
        if orientations[td[i]]!=parity(perm)*orientations[i]:raise ValueError('deck reverses the source orientation')
    reps=[i for i,j in enumerate(td) if i<j];ri={v:i for i,v in enumerate(reps)}
    qglu=[];done=set()
    for a,i in enumerate(reps):
        for k in range(4):
            if (a,k) in done:continue
            j,l=neighbor[i,k];flip=j>td[j];r=td[j] if flip else j;b=ri[r]
            shared=[deck[v] if flip else v for h,v in enumerate(facets[i]) if h!=k]
            target=facets[r];ll=next(h for h,v in enumerate(target) if v not in shared)
            perm=[]
            for h,v in enumerate(facets[i]):
                perm.append(ll if h==k else target.index(deck[v] if flip else v))
            qglu.append([a,k,b,perm]);done.add((a,k));done.add((b,ll))
    gluing_paths=[]
    for b in blocks[1:]:
        interface,cyc,A=b.cycles[0];paths={}
        for name in ('meridian','longitude'):
            path=[];lift=[]
            for k in range(N):
                va,vb=((cyc[k],0),(cyc[(k+1)%N],0)) if name=='meridian' else ((cyc[0],k),(cyc[0],(k+1)%N))
                e=b.edge(va,vb);start,end=b.endpoints(e)
                pts=dict(edge_points[e]);pts[Q(0)]=start;pts[Q(1)]=end
                points=sorted(pts.items())
                if b.v(*va)!=start:points=[(1-t,v) for t,v in reversed(points)]
                for t,v in points:
                    u,w=(Q(k,N)+t/N,Q(0)) if name=='meridian' else (Q(0),Q(k,N)+t/N)
                    q=(A[0][0]*u+A[0][1]*w,A[1][0]*u+A[1][1]*w)
                    if path and path[-1]==ix[v]:continue
                    path.append(ix[v]);lift.append([str(z) for z in q])
            paths[name]={'source_vertices':path,'unwrapped_target_coordinates':lift,
                'total_displacement':[A[0][0],A[1][0]] if name=='meridian' else [A[0][1],A[1][1]],
                'closed_modulo_Z2':path[0]==path[-1]}
        gluing_paths.append({'interface':interface,'matrix':A,'basis_columns':paths,
            'scope':'boundary of each rational-tangle filling disk and a chosen filling longitude; these are not the knot peripheral basis'})
    model={'schema':'marked-montesinos-seifert-source-v1','pretzel_parameters':list(params),
        'rational_tangle_slopes':['1/'+str(p) for p in params],
        'gluing_matrices':[[[p,1],[1,0]] for p in params],
        'coordinates':'central: x=oriented pants boundary, y=fibre; filling: u=positive meridian, v=longitude; (x,y)=A_p(u,v)',
        'orientation_convention':'A_p has determinant -1; central oriented boundary basis (q_i,h); product-orientation gluings reverse boundary orientation',
        'deck_convention':'pants reflection times fibre reflection; filling disk conjugation times longitude conjugation; torus (x,y)->(-x,-y)',
        'vertices':verts,'facets':facets,'closed_cover_facets':facets,'deck':deck,'tetrahedron_deck':td,
        'tetrahedron_orientations':orientations,'quotient_tetrahedron_orientations':[orientations[i] for i in reps],
        'orientation_source':'pants positive hexagon +, reflected hexagon -; filling disks positive; circle factors positive; exact oriented product-boundary chains and positive cone extension, not arbitrary BFS orientation',
        'vertex_ranks':[min(i,deck[i]) for i in range(len(deck))],
        'rank_scope':'deck-invariant ranks; distinct on every tetrahedron, checked by the no-fold condition',
        'ramification_edges':fixededges,'ramification_components':fixed_components,
        'parts':{s:sorted(tuple(sorted(ix[v] for v in t)) for t in ts) for s,ts in parts.items()},
        'overlay_polygons':overlay_records,'quotient_representatives':reps,'gluing_paths':gluing_paths,
        'blocks':[{'name':b.name,'surface_triangles':b.triangles,'surface_reflection':b.reflection,
                   'boundary_cycles':b.cycles,'circle_edges':N} for b in blocks],
        'facets_sha256':digest(facets),'source_interpretation':'explicit Montesinos pair construction; sign/endpoint comparison with external diagram is a separate recorded gate',
        'count_summary':{'vertices':len(verts),'tetrahedra':len(facets),'quotient_tetrahedra':len(reps),'overlay_polygons':total_polygons,'fixed_components':len(fixed_components)}}
    return model,{'simplex_count':len(facets),'gluings':gluings}, {'simplex_count':len(reps),'gluings':qglu}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=Path(__file__).parent/'data')
    args=parser.parse_args();args.output.mkdir(exist_ok=True,parents=True)
    if (args.output/'source.json.gz').exists():raise FileExistsError('Use a fresh output directory')
    t=time.monotonic();model,cover,quotient=build()
    save(args.output/'source.json.gz',model);save(args.output/'cover_regina.json.gz',cover);save(args.output/'quotient_regina.json.gz',quotient)
    result={'status':'FINITE_FACET_CONSTRUCTION_COMPLETE_NATIVE_NOT_YET_ADMITTED',**model['count_summary'],'seconds':time.monotonic()-t,'facets_sha256':model['facets_sha256']}
    save(args.output/'build_summary.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()

"""Separate finite checker for the chart-aware roll assembly.
No clipping, ear triangulation, cylinder builder, or 4D producer is imported.
The geometric interpretation uses the explicit PL straightening in THEORY.md.
"""
from __future__ import annotations
from fractions import Fraction as Q
from collections import Counter,defaultdict,deque
from itertools import combinations,product
from pathlib import Path
import gzip,json,hashlib,time,sys
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'

def load(p):
 p=Path(p);b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def dig(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def tup(x):return tuple(tup(z) for z in x) if isinstance(x,(tuple,list)) else x

def req(x,m):
 if not x:raise ValueError(m)
def mod(p):return tuple(z-z.numerator//z.denominator for z in p)
def det2(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def ar2(P):return sum(P[i][0]*P[(i+1)%len(P)][1]-P[(i+1)%len(P)][0]*P[i][1] for i in range(len(P)))
def inside(p,A):
 s=1 if ar2(A)>0 else -1
 return all(s*det2(a,b,p)>=0 for a,b in zip(A,A[1:]+A[:1]))
def hull(P):
 P=sorted(set(P));L=[];U=[]
 for x in P:
  while len(L)>1 and det2(L[-2],L[-1],x)<=0:L.pop()
  L.append(x)
 for x in reversed(P):
  while len(U)>1 and det2(U[-2],U[-1],x)<=0:U.pop()
  U.append(x)
 return L[:-1]+U[:-1]
def intersect_vertices(A,B):
 P={x for x in A if inside(x,B)}|{x for x in B if inside(x,A)}
 for a,b in zip(A,A[1:]+A[:1]):
  for c,d in zip(B,B[1:]+B[:1]):
   x=det2(a,b,c);y=det2(a,b,d)
   if x==y:continue
   f=x/(x-y)
   if 0<=f<=1:
    p=tuple(c[k]+f*(d[k]-c[k]) for k in (0,1))
    if all(min(a[k],b[k])<=p[k]<=max(a[k],b[k]) for k in (0,1)):P.add(p)
 return hull(P)

def overlay_check(A,B,O,deck_shift=None):
 ap=[[tuple(map(Q,p)) for p in f['points']] for f in A];bp=[[tuple(map(Q,p)) for p in f['points']] for f in B]
 aa=[Q(0)]*len(A);ab=[Q(0)]*len(B);va=[set() for _ in A];vb=[set() for _ in B];seen=set();rt=[]
 key=(lambda p:(p[0],p[1]%deck_shift)) if deck_shift else lambda p:p
 for piece in O['pieces']:
  i,j=piece['a'],piece['b'];shift=piece['shift'];ident=(i,j,tuple(shift));req(ident not in seen,'duplicate overlay intersection');seen.add(ident)
  BB=[tuple(v[k]+shift[k] for k in (0,1)) for v in bp[j]];P=[tuple(map(Q,p)) for p in piece['polygon']]
  req(hull(P)==intersect_vertices(ap[i],BB),'wrong exact intersection polygon')
  req(ar2(P)>0 and all(inside(p,ap[i]) and inside(p,BB) for p in P),'intersection orientation/carrier')
  AA=Q(0);tri=[]
  for t in piece['triangles']:
   t=[tuple(map(Q,p)) for p in t];req(det2(*t)>0 and set(t)<=set(P),'overlay triangle')
   for old in tri:
    H=intersect_vertices(t,old)
    req(len(H)<3 or ar2(H)==0,'overlapping overlay triangles')
   tri.append(t);AA+=det2(*t);rt.append(tuple(sorted((mod(p) for p in t),key=key)))
  req(AA==ar2(P),'incomplete polygon triangulation')
  aa[i]+=AA;ab[j]+=AA;va[i].update(P);vb[j].update(tuple(p[k]-shift[k] for k in (0,1)) for p in P)
 for i,p in enumerate(ap):req(aa[i]==abs(ar2(p)),'incomplete A overlay')
 for j,p in enumerate(bp):req(ab[j]==abs(ar2(p)),'incomplete B overlay')
 req(len(rt)==len({frozenset(t) for t in rt}),'duplicate common-refinement triangle')
 req(rt==[tuple(tuple(map(Q,p)) for p in t) for t in O['triangles']],'common-refinement triangle list')
 for keyname,vs in [('a_vertices',va),('b_vertices',vb)]:req(vs==[{tuple(map(Q,p)) for p in row} for row in O[keyname]],'overlay edge-vertex table')
 return {'a':va,'b':vb,'key':key,'refined_triangles':len(rt),'polygons':len(seen)}

def grid(L,M,n=0):
 fs=[]
 for i in range(L):
  for j in range(M):
   for path in [((i,j),(i+1,j),(i+1,j+1)),((i,j),(i,j+1),(i+1,j+1))]:
    ids=tuple((a%L,b%M) for a,b in path);points=tuple((Q(a,L)+n*Q(b,M),Q(b,M)) for a,b in path)
    fs.append({'ids':ids,'points':points,'order':ids})
 return fs

def same_coarse(got,want):
 req(len(got)==len(want),'coarse face count')
 for a,b in zip(got,want):
  req(tup(a['ids'])==tup(b['ids']) and tup(a['order'])==tup(b['order']),'coarse vertex labels/order')
  req([tuple(map(Q,p)) for p in a['points']]==list(b['points']),'coarse geometric chart')

def det3(cols):
 a,b,c=cols;return a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+c[0]*(a[1]*b[2]-a[2]*b[1])

def half_check(F,O,oc,side,charts,cv,rv,iv):
 children=defaultdict(list)
 for p in O['pieces']:
  sh=p['shift'] if side=='b' else (0,0)
  for tr in p['triangles']:children[p[side]].append([tuple(Q(v[k])-sh[k] for k in (0,1)) for v in tr])
 key=oc['key'];want={};totvolume=Q(0)
 for i,f in enumerate(F):
  pp={tup(v):tuple(map(Q,p)) for v,p in zip(f['ids'],f['points'])};order=tup(f['order']);C={cv(v):(*pp[v],Q(0)) for v in pp};bd=[tuple(cv(v) for v in order)]
  for tr in children[i]:
   tr.sort(key=lambda p:key(mod(p)));ids=tuple(rv(mod(p)) for p in tr);bd.append(ids)
   for v,p in zip(ids,tr):C[v]=(*p,Q(1))
  for i0,i1 in [(0,1),(0,2),(1,2)]:
   a,b=order[i0],order[i1];pa,pb=pp[a],pp[b];axis=0 if pa[0]!=pb[0] else 1
   edge=[p for p in oc[side][i] if det2(pa,pb,p)==0 and all(min(pa[k],pb[k])<=p[k]<=max(pa[k],pb[k]) for k in (0,1))]
   edge.sort(key=lambda p:(p[axis]-pa[axis])/(pb[axis]-pa[axis]));req(edge[0]==pa and edge[-1]==pb,'side endpoints')
   bd.append((cv(a),cv(b),rv(mod(pb))));C[rv(mod(pb))]=(*pb,Q(1))
   for x,y in zip(edge,edge[1:]):
    z=sorted([x,y],key=lambda p:key(mod(p)));bd.append((cv(a),*(rv(mod(p)) for p in z)))
    for p in z:C[rv(mod(p))]=(*p,Q(1))
  centre=iv(i);C[centre]=tuple(sum(p[k] for p in pp.values())/3 for k in (0,1))+(Q(1,2),)
  for tr in bd:
   te=(*tr,centre);co=tuple(C[v] for v in te);want[te]=(i,co)
   p=co[0];cols=[tuple(x[k]-p[k] for k in range(3)) for x in co[1:]];vol=abs(det3(cols))/6;req(vol>0,'degenerate prism tetrahedron');totvolume+=vol
 req(len(charts)==len(want),'prism tets count')
 actual=[]
 for c in charts:
  te=tup(c['vertices']);co=tuple(tuple(map(Q,p)) for p in c['coords']);req(te in want and want[te]==(c['parent'],co),'prism boundary cone mismatch');actual.append(te)
 req(len(set(actual))==len(want),'prism repetition');req(totvolume==1,'prism volumes not torus x I')
 return actual


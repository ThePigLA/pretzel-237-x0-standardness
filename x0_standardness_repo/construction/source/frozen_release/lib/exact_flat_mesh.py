"""Exact flat-torus overlay and branched prism collars. Fractions throughout."""
from fractions import Fraction as Q
from itertools import product
from collections import defaultdict
import math

def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def det(a,b,c):return cross(sub(b,a),sub(c,a))
def mod(p):return tuple(z-z.numerator//z.denominator for z in p)
def area2(P):return sum(cross(P[i],P[(i+1)%len(P)]) for i in range(len(P)))
def clean(P):
 out=[]
 for p in P:
  if not out or p!=out[-1]:out.append(p)
 if len(out)>1 and out[0]==out[-1]:out.pop()
 return out

def clip(P,A):
 if area2(A)<0:A=list(reversed(A))
 for a,b in zip(A,A[1:]+A[:1]):
  if not P:return []
  out=[]
  for c,d in zip(P,P[1:]+P[:1]):
   x,y=det(a,b,c),det(a,b,d)
   if x>=0:out.append(c)
   if (x<0<y) or (y<0<x):
    t=x/(x-y);out.append(tuple(c[k]+t*(d[k]-c[k]) for k in range(2)))
  P=clean(out)
 if len(P)<3 or not area2(P):return []
 if area2(P)<0:P.reverse()
 return P

def ear(P,key):
 """Deterministic convex-polygon triangulation retaining collinear side vertices."""
 P=list(P);out=[]
 if area2(P)<0:P.reverse()
 while len(P)>3:
  cand=[]
  for i,b in enumerate(P):
   a,c=P[i-1],P[(i+1)%len(P)]
   if det(a,b,c)<=0:continue
   if any(det(a,b,x)>=0 and det(b,c,x)>=0 and det(c,a,x)>=0 for j,x in enumerate(P) if j not in {(i-1)%len(P),i,(i+1)%len(P)}):continue
   cand.append(i)
  if not cand:raise ValueError('no nondegenerate ear')
  i=min(cand,key=lambda i:key(mod(P[i])))
  out.append((P[i-1],P[i],P[(i+1)%len(P)]));P.pop(i)
 if det(*P)<=0:raise ValueError('degenerate last ear')
 return out+[tuple(P)]

def grid_faces(L=3,M=4,shear=0):
 out=[]
 for i in range(L):
  for j in range(M):
   def pt(a,b):return (Q(a,L)+shear*Q(b,M),Q(b,M))
   for ij in [((i,j),(i+1,j),(i+1,j+1)),((i,j),(i,j+1),(i+1,j+1))]:
    ids=tuple((a%L,b%M) for a,b in ij);ps=tuple(pt(a,b) for a,b in ij)
    out.append({'ids':ids,'points':ps,'order':ids})
 return out

def overlay(A,B,deck_shift=None):
 key=(lambda p:(p[0],p[1]%deck_shift)) if deck_shift else lambda p:p
 pieces=[];rs=[];a_verts=[set() for _ in A];b_verts=[set() for _ in B];areasA=[Q(0) for _ in A];areasB=[Q(0) for _ in B]
 bbox=lambda ps:tuple((min(p[k] for p in ps),max(p[k] for p in ps)) for k in range(2))
 boxesB=[bbox(b['points']) for b in B]
 for i,a in enumerate(A):
  ap=list(a['points']);aa=bbox(ap)
  for j,b in enumerate(B):
   bb=boxesB[j]
   # Integer translations of B meeting A, including boundary contacts.
   ranges=[range(math.ceil(aa[k][0]-bb[k][1]),math.floor(aa[k][1]-bb[k][0])+1) for k in range(2)]
   for shift in product(*ranges):
    bp=[tuple(p[k]+shift[k] for k in range(2)) for p in b['points']]
    P=clip(ap,bp)
    if not P:continue
    if len({key(mod(p)) for p in P})<len(P):raise ValueError('deck orbit repeated in polygon')
    tri=ear(P,key);idx=len(pieces)
    pieces.append({'a':i,'b':j,'shift':shift,'polygon':P,'triangles':tri})
    ar=abs(area2(P));areasA[i]+=ar;areasB[j]+=ar
    a_verts[i].update(P);b_verts[j].update(tuple(p[k]-shift[k] for k in range(2)) for p in P)
    for t in tri:
     verts=tuple(sorted((mod(p) for p in t),key=key));rs.append(verts)
 for f,ar in zip(A,areasA):
  if abs(area2(f['points']))!=ar:raise ValueError('A overlay area')
 for f,ar in zip(B,areasB):
  if abs(area2(f['points']))!=ar:raise ValueError('B overlay area')
 if len(rs)!=len({frozenset(t) for t in rs}):raise ValueError('duplicate refinement triangles')
 # Parent-child incidences and rational chart points are preserved.
 return {'pieces':pieces,'triangles':rs,'a_vertices':a_verts,'b_vertices':b_verts,'key':key}

def cylinder_half(coarse,O,side,coarse_vertex,refined_vertex,centre_vertex):
 """Coarse torus at h=0, common subdivision at h=1; exact coned triangle prisms.
 Returns locally ordered tetrahedra and universal charts (x,y,h).
 """
 key=O['key'];out=[];charts=[];parent_tris=defaultdict(list)
 for p in O['pieces']:
  parent=p[side];shift=p['shift'] if side=='b' else (0,0)
  for tr in p['triangles']:
   ts=[tuple(v[k]-shift[k] for k in range(2)) for v in tr]
   parent_tris[parent].append(ts)
 for i,f in enumerate(coarse):
  pp=dict(zip(f['ids'],f['points']));order=f['order'];cs={coarse_vertex(v):(*p,Q(0)) for v,p in pp.items()}
  tri=[]
  tri.append(tuple(coarse_vertex(v) for v in order))
  for r in parent_tris[i]:
   r=sorted(r,key=lambda p:key(mod(p)));tr=tuple(refined_vertex(mod(p)) for p in r)
   tri.append(tr)
   for v,p in zip(tr,r):cs[v]=(*p,Q(1))
  for ia,ib in ((0,1),(0,2),(1,2)):
   a,b=order[ia],order[ib];pa,pb=pp[a],pp[b];d=sub(pb,pa)
   seg=[p for p in O[side+'_vertices'][i] if det(pa,pb,p)==0 and all(min(pa[k],pb[k])<=p[k]<=max(pa[k],pb[k]) for k in range(2))]
   k=0 if d[0] else 1;seg.sort(key=lambda p:(p[k]-pa[k])/d[k])
   if not seg or seg[0]!=pa or seg[-1]!=pb:raise ValueError('edge subdivision endpoints')
   ca,cb=coarse_vertex(a),coarse_vertex(b);rb=refined_vertex(mod(pb));tri.append((ca,cb,rb));cs[rb]=(*pb,Q(1))
   for p,q in zip(seg,seg[1:]):
    tt=sorted([p,q],key=lambda p:key(mod(p)));tr=(ca,*(refined_vertex(mod(x)) for x in tt));tri.append(tr)
    for z in tt:cs[refined_vertex(mod(z))]=(*z,Q(1))
  centre=centre_vertex(i);cp=tuple(sum(p[k] for p in f['points'])/3 for k in range(2))+(Q(1,2),);cs[centre]=cp
  for tr in tri:
   tet=(*tr,centre);out.append(tet);charts.append({'parent':i,'vertices':tet,'coords':tuple(cs[v] for v in tet)})
 return out,charts

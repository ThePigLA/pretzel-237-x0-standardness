"""Deterministic equivariant interior-star coarsening before the time product.
Producer only. A separate local-ball and face-incidence implementation replays it.
No native library is used and no exterior/tube is replaced by a named model.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from itertools import combinations
import hashlib,heapq,json,time
from branch_dpll import solve

def tup(x):return tuple(tup(z) for z in x) if isinstance(x,(list,tuple)) else str(x) if isinstance(x,Q) else x

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def determinant(C):
 a,b,c=C
 return a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+c[0]*(a[1]*b[2]-a[2]*b[1])

def bary(P,co):
 a=co[0];cols=[tuple(p[k]-a[k] for k in range(3)) for p in co[1:]];den=determinant(cols)
 if not den:return None
 rhs=tuple(P[k]-a[k] for k in range(3));w=[]
 for k in range(3):z=cols.copy();z[k]=rhs;w.append(determinant(z)/den)
 w=[1-sum(w),*w]
 return w if min(w)>=0 else None

def support_tags(ac,chart,deck):
 """Freeze the old actual source-quarter witnesses and old peripheral edge paths."""
 witnesses=[];protect=set()
 for n in (0,1):
  for s in (Q(0),Q(1,4),Q(3,4),Q(1)):
   r=3*s;side='N' if r<=1 else 'E';h=min(Q(1),max(Q(0),r-1));u,v=Q(1,17),Q(1,7)
   for time,x in ((1,u+n*(s-h)),(0,u+n*s)):
    found=None
    for half in ('a','b'):
     for cell in ac['charts'][side][half]:
      co=[tuple(map(Q,p)) for p in cell['coords']]
      if side=='N':co=[(p[0],p[1],p[2]/2 if half=='a' else 1-p[2]/2) for p in co]
      else:co=[(p[0],p[1],3-p[2]/2 if half=='a' else 2+p[2]/2) for p in co]
      for i in (-1,0,1):
       for j in (-1,0,1):
        P=(x+i,v+j,r)
        if any(not min(p[k] for p in co)<=P[k]<=max(p[k] for p in co) for k in range(3)):continue
        w=bary(P,co)
        if w is not None:found={'roll':n,'input_s':str(s),'time':time,'side':side,'half':half,'vertices':tup(cell['vertices']),'coordinates':[[str(a) for a in p] for p in co],'barycentric_weights':list(map(str,w)),'lattice_translation':[i,j]};break
       if found:break
      if found:break
     if found:break
    if not found:raise ValueError('unlocated frozen witness')
    witnesses.append(found);protect.update(found['vertices'])
 per=set()
 for loop in chart['peripheral_loops'].values():
  for e in loop['torus_edge_loop']:
   a,b=chart['torus_downstairs']['edge_representatives'][abs(e)-1];per|={('E','o',a),('E','o',b)}
 per.add(('E','o',chart['torus_downstairs']['vertices'][0]))
 protect|=per;protect|={deck[v] for v in list(protect)}
 return protect,witnesses,sorted(per,key=repr)

class Complex:
 def __init__(self,tets,deck,protected):
  self.K=set(tets);self.st=defaultdict(set);self.D=deck;self.protected=set(protected);self.cache={};self.rev=defaultdict(int)
  for t in self.K:
   for v in t:self.st[v].add(t)
 def lk(self,u):
  if u not in self.cache:
   tris={tuple(v for v in t if v!=u) for t in self.st[u]}
   self.cache[u]=[{v for f in tris for v in f},{e for f in tris for e in combinations(f,2)},tris]
  return self.cache[u]
 def proposal(self,u,v):
  if u in self.protected or not self.st[u] or not self.st[v]:return None
  U,V=self.lk(u),self.lk(v)
  if v not in U[0]:return None
  ed={tuple(w for w in t if w not in (u,v)) for t in self.st[u] if v in t}
  if U[0]&V[0]!={z for e in ed for z in e} or U[1]&V[1]!=ed or U[2]&V[2]:return None
  old=set(self.st[u]);new={tuple(sorted(v if w==u else w for w in t)) for t in old if v not in t}
  if len(new)!=sum(v not in t for t in old) or new&(self.K-old):return None
  if any(len({min(w,self.D[w]) for w in t})!=4 for t in new):return None
  return old,new
 def apply(self,old,new):
  affected={v for t in old|new for v in t}
  for t in old:
   self.K.remove(t)
   for v in t:self.st[v].remove(t)
  for t in new:
   self.K.add(t)
   for v in t:self.st[v].add(t)
  for v in affected:self.cache.pop(v,None);self.rev[v]+=1
  return affected
 def run(self,seconds=120):
  start=time.monotonic();pq=[];steps=[]
  def push(v):
   if v not in self.protected and self.st[v] and v<self.D[v]:heapq.heappush(pq,(len(self.st[v]),v,self.rev[v]))
  for v in sorted(self.st):push(v)
  while pq:
   if time.monotonic()-start>seconds:raise RuntimeError('COARSENING_CAP_NO_BUILD')
   _,u,ver=heapq.heappop(pq)
   if ver!=self.rev[u] or not self.st[u] or u in self.protected:continue
   ns=sorted(self.lk(u)[0],key=lambda v:(v not in self.protected,-len(self.st[u]&self.st[v]),-len(self.st[v]),v))
   for v in ns:
    U,V=self.D[u],self.D[v]
    if U==u or V==v or V==u or U in self.lk(u)[0]:continue
    p=self.proposal(u,v);p2=self.proposal(U,V)
    if p is None or p2 is None:continue
    count=len(self.K);aft=self.apply(*p);pp=self.proposal(U,V)
    if pp is None:self.apply(p[1],p[0]);continue
    aft|=self.apply(*pp)
    steps.append({'u':u,'v':v,'deck_u':U,'deck_v':V,'before':count,'after':len(self.K),'deleted_per_side':[len(p[0]),len(pp[0])],'inserted_per_side':[len(p[1]),len(pp[1])]})
    for w in sorted(aft):push(w)
    break
  return steps

def branching(K,initial,D,protected):
 edges=sorted({e for t in K for e in combinations(t,2)});tris=sorted({f for t in K for f in combinations(t,3)})
 def rep(e):return min(e,tuple(sorted(D[v] for v in e)))
 reps=sorted(set(map(rep,edges)));vi={e:i for i,e in enumerate(reps)}
 def lit(a,b):
  e=tuple(sorted((a,b)));rr=rep(e);positive=(a<b) if e==rr else (D[a]<D[b]);return (vi[rr]+1)*(1 if positive else -1)
 frozen={}
 for t in initial:
  for a,b in combinations(t,2):
   if a in protected and b in protected:frozen[min(a,b),max(a,b)]=a<b
 clauses=[]
 for a,b,c in tris:
  z=[lit(a,b),lit(b,c),lit(c,a)];clauses.extend((z,[-l for l in z]))
 for (a,b),orientation in sorted(frozen.items()):
  if (a,b) not in set(edges):raise ValueError('lost protected edge')
  clauses.append([lit(a,b)*(1 if orientation else -1)])
 sol,stats=solve(clauses,len(reps),seconds=60)
 def before(a,b):l=lit(a,b);return int((sol[abs(l)-1]>0)==(l>0))
 ordered=[]
 for t in sorted(K):
  out=tuple(sorted(t,key=lambda v:sum(before(w,v) for w in t if w!=v)))
  if not all(before(a,b) for a,b in combinations(out,2)):raise ValueError('cyclic local branching')
  ordered.append(out)
 # Variable elapsed time is excluded from deterministic proof objects.
 stats.pop('seconds',None)
 return ordered,stats

def construct(stages,deck,ac,chart):
 stages={side:[tup(t) for t in ts] for side,ts in stages.items()};deck={tup(a):tup(b) for a,b in deck.items()}
 protect,witnesses,peripheral=support_tags(ac,chart,deck)
 cert={'schema':'equivariant-relative-3D-coarsening-v1','parts':{},'protected_source_witnesses':witnesses,'protected_source_peripheral_vertices':peripheral,'topology_rule':'interior-star ball replacement with the identical boundary; disjoint deck-related interiors; protected subcomplex pointwise fixed'}
 compact={}
 for side in ('N','E'):
  tags=sorted({v for t in stages[side] for v in t},key=repr);idx={v:i for i,v in enumerate(tags)};D=[idx[deck[v]] for v in tags]
  initial=[tuple(idx[v] for v in t) for t in stages[side]];K={tuple(sorted(t)) for t in initial};faces=Counter(f for t in K for f in combinations(t,3));bd={f for f,c in faces.items() if c==1}
  frozen={v for f in bd for v in f}|{i for i,d in enumerate(D) if i==d}|{idx[v] for v in protect if v in idx};frozen|={D[v] for v in list(frozen)}
  worker=Complex(K,D,frozen);moves=worker.run();orders,bstats=branching(worker.K,initial,D,frozen)
  result={'tags':tags,'deck':D,'protected':sorted(frozen),'initial':sorted(K),'initial_ordered':initial,'initial_sha256':digest(sorted(K)),'endpoint':sorted(worker.K),'endpoint_sha256':digest(sorted(worker.K)),'ordered_endpoint':orders,'steps':moves,'branching_search':bstats}
  cert['parts'][side]=result;compact[side]=[tuple(tags[v] for v in t) for t in orders]
  print('3D COARSEN',side,len(K),'->',len(worker.K),'deck-pair operations',len(moves),flush=True)
 return compact,deck,cert

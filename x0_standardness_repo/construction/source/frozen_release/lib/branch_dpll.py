"""Tiny deterministic SAT search; output is only an orientation witness, checked separately."""
import time

def solve(clauses,n,preferred=None,max_nodes=200000,seconds=30):
 clauses=tuple(tuple(sorted(set(c),key=lambda z:(abs(z),z<0))) for c in clauses if not any(-x in c for x in c))
 clauses=tuple(sorted(set(clauses),key=lambda c:(len(c),c)))
 preferred=preferred or [1]*n
 start=time.monotonic();nodes=0
 def rec(a):
  nonlocal nodes
  nodes+=1
  if nodes>max_nodes or time.monotonic()-start>seconds:raise RuntimeError('BRANCHING_SEARCH_CAP')
  while True:
   units={};opencls=[]
   for c in clauses:
    rem=[];sat=False
    for l in c:
     v=a[abs(l)-1]
     if v and (v>0)==(l>0):sat=True;break
     if not v:rem.append(l)
    if sat:continue
    if not rem:return None
    if len(rem)==1:
     l=rem[0];i=abs(l)-1;v=1 if l>0 else -1
     if i in units and units[i]!=v:return None
     units[i]=v
    else:opencls.append(rem)
   if not units:break
   for i,v in units.items():
    if a[i] and a[i]!=v:return None
    a[i]=v
  if not opencls:
   return [v or preferred[i] for i,v in enumerate(a)]
  shortest=min(map(len,opencls));scores={}
  for c in opencls:
   if len(c)!=shortest:continue
   for l in c:scores[abs(l)-1]=scores.get(abs(l)-1,0)+1
  i=min(scores,key=lambda i:(-scores[i],i))
  for v in (preferred[i],-preferred[i]):
   b=a.copy();b[i]=v;ans=rec(b)
   if ans is not None:return ans
  return None
 ans=rec([0]*n)
 if ans is None:raise ValueError('NO_BRANCHING_FOUND')
 assert all(any((ans[abs(l)-1]>0)==(l>0) for l in c) for c in clauses)
 return ans,{'nodes':nodes,'seconds':time.monotonic()-start,'clauses':len(clauses),'variables':n}

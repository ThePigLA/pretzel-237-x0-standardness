"""Independent replayer: replace certified 3-balls along identical triangulated S2.
Does NOT import the contraction producer, its adjacency caches, or branching solver.
The geometric interpretation uses the relative PL ball-replacement lemma in THEORY.
"""
from __future__ import annotations
from collections import Counter,defaultdict,deque
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json,time

def req(x,m):
 if not x:raise ValueError(m)
def tup(x):return tuple(tup(z) for z in x) if isinstance(x,(list,tuple)) else str(x) if isinstance(x,Q) else x
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def faces(tets):return {f for t in tets for k in range(1,5) for f in combinations(sorted(t),k)}
def facecount(tets):return Counter(f for t in tets for k in range(1,5) for f in combinations(sorted(t),k))
def boundary(tets):return {f for f,c in Counter(f for t in tets for f in combinations(sorted(t),3)).items() if c==1}

def connected(vertices,edges):
 if not vertices:return False
 ad=defaultdict(set)
 for a,b in edges:ad[a].add(b);ad[b].add(a)
 seen={next(iter(vertices))};todo=list(seen)
 for v in todo:
  for w in ad[v]:
   if w not in seen:seen.add(w);todo.append(w)
 return seen==set(vertices)

def surface(tris,disk=False):
 tris=set(tris);V={v for t in tris for v in t};ec=Counter(e for t in tris for e in combinations(sorted(t),2));E=set(ec)
 req(tris and all(len(t)==3 and len(set(t))==3 for t in tris),'bad surface triangles')
 req(all(c in ((1,2) if disk else (2,)) for c in ec.values()),'surface edge incidence')
 req(connected(V,E),'disconnected surface')
 req(len(V)-len(E)+len(tris)==(1 if disk else 2),'surface Euler')
 bd={e for e,c in ec.items() if c==1};bv={v for e in bd for v in e}
 if disk:
  bc=Counter(v for e in bd for v in e);req(bd and all(c==2 for c in bc.values()) and connected(bv,bd),'disk boundary is not one circle')
 for v in V:
  le={tuple(w for w in t if w!=v) for t in tris if v in t};lv={w for e in le for w in e};dc=Counter(w for e in le for w in e)
  req(connected(lv,le),'surface vertex pinching')
  if v in bv:req(sum(c==1 for c in dc.values())==2 and all(c in (1,2) for c in dc.values()),'boundary vertex link is not a path')
  else:req(all(c==2 for c in dc.values()),'interior surface vertex link is not a circle')
 return (len(V),len(E),len(tris))

def check_protection(cert,stages,deck,ac,chart):
 tags={v for ts in stages.values() for t in ts for v in t};protected=set();wits=cert['protected_source_witnesses'];req(len(wits)==16,'witness coverage')
 seen=set()
 for w in wits:
  n=w['roll'];s=Q(w['input_s']);t=w['time'];req(n in (0,1) and s in (0,Q(1,4),Q(3,4),1) and t in (0,1),'witness input')
  k=(n,s,t);req(k not in seen,'repeated witness');seen.add(k)
  r=3*s;side='N' if r<=1 else 'E';req(w['side']==side,'witness side')
  vv=tup(w['vertices']);req(set(vv)<=tags and vv in [tup(z) for z in stages[side]],'witness not an original tetrahedron')
  # Verify a source chart, not a claimed endpoint coordinate.
  found=[c for c in ac['charts'][side][w['half']] if tup(c['vertices'])==vv];req(len(found)==1,'witness original chart')
  co=[tuple(map(Q,p)) for p in found[0]['coords']]
  if side=='N':co=[(p[0],p[1],p[2]/2 if w['half']=='a' else 1-p[2]/2) for p in co]
  else:co=[(p[0],p[1],3-p[2]/2 if w['half']=='a' else 2+p[2]/2) for p in co]
  req(co==[tuple(map(Q,p)) for p in w['coordinates']],'witness geometric chart changed')
  weights=list(map(Q,w['barycentric_weights']));req(len(weights)==4 and min(weights)>=0 and sum(weights)==1,'witness barycentric weights')
  h=min(Q(1),max(Q(0),r-1));x=Q(1,17)+n*((s-h) if t else s);i,j=w['lattice_translation'];req(isinstance(i,int) and isinstance(j,int),'witness lattice')
  req(tuple(sum(weights[z]*co[z][k] for z in range(4)) for k in range(3))==(x+i,Q(1,7)+j,r),'witness not the prescribed source point')
  protected.update(vv)
 per={('E','o',chart['torus_downstairs']['vertices'][0])}
 for loop in chart['peripheral_loops'].values():
  for e in loop['torus_edge_loop']:
   per.update(('E','o',v) for v in chart['torus_downstairs']['edge_representatives'][abs(e)-1])
 req(per==set(map(tup,cert['protected_source_peripheral_vertices'])),'source peripheral paths changed')
 protected|=per;protected|={deck[v] for v in list(protected)}
 return protected

def replay_part(p,initial_tags_tets,dt,required):
 tags=list(map(tup,p['tags']));req(tags==sorted({v for t in initial_tags_tets for v in t},key=repr),'initial vertex numbering')
 ix={v:i for i,v in enumerate(tags)};D=[ix[dt[v]] for v in tags];req(D==p['deck'] and all(D[D[i]]==i for i in range(len(D))),'3D deck')
 initial_ordered=[tuple(ix[v] for v in t) for t in initial_tags_tets];initial={tuple(sorted(t)) for t in initial_ordered};req(sorted(initial)==p['initial'] or sorted(initial)==list(map(tuple,p['initial'])),'initial tetrahedra')
 req(list(map(tuple,p['initial_ordered']))==initial_ordered,'initial branching')
 req(digest(sorted(initial))==p['initial_sha256'],'initial digest')
 bd=boundary(initial);protected=set(p['protected']);need={v for f in bd for v in f}|{i for i,d in enumerate(D) if i==d}|{ix[v] for v in required if v in ix}
 need|={D[v] for v in need};req(need<=protected and all(D[v] in protected for v in protected),'unprotected mark or boundary')
 live=set(initial);counts=facecount(live);original_fixed={f for f in counts if all(D[v]==v for v in f)};stats=Counter()
 def single(u,v):
  req(u!=v and u not in protected and D[u]!=u and D[v]!=v,'protected or nonfree contraction')
  old={t for t in live if u in t};req(old and any(v in t for t in old),'contraction support missing')
  L={tuple(w for w in t if w!=u) for t in old};surface(L)
  disk={t for t in L if v not in t};surface(disk,disk=True)
  new={tuple(sorted((v,*t))) for t in disk};req(new and boundary(new)==L,'new ball boundary differs')
  oldc=facecount(old);Lc={f for t in L for k in (1,2,3) for f in combinations(sorted(t),k)}
  req(not (new&(live-old)),'duplicate replacement tetrahedron')
  for f in faces(new):
   req(counts[f]-oldc[f]==0 or f in Lc,'new ball meets complement off its boundary')
  req(all(len({min(z,D[z]) for z in t})==4 for t in new),'new cell folds under deck quotient')
  for f,c in oldc.items():
   counts[f]-=c
   if not counts[f]:del counts[f]
  counts.update(facecount(new));live.difference_update(old);live.update(new)
  stats['ball_replacements']+=1;stats['spherical_links_checked']+=1;stats['disk_antistars_checked']+=1
  return old,new
 for i,s in enumerate(p['steps']):
  u,v,U,V=(s[z] for z in ('u','v','deck_u','deck_v'))
  req(U==D[u] and V==D[v] and u!=U and v!=V and U!=v and V!=u,'paired vertex map')
  req(not counts.get(tuple(sorted((u,U))),0),'paired stars have overlapping interiors')
  req(len(live)==s['before'],'recorded source size')
  old,new=single(u,v);old2,new2=single(U,V)
  mirror=lambda ts:{tuple(sorted(D[z] for z in t)) for t in ts}
  req(old2==mirror(old) and new2==mirror(new),'move pair is not exactly deck-conjugate')
  req([len(old),len(old2)]==s['deleted_per_side'] and [len(new),len(new2)]==s['inserted_per_side'],'recorded local counts')
  req(len(live)==s['after'],'recorded target size')
 req(sorted(live)==list(map(tuple,p['endpoint'])) and digest(sorted(live))==p['endpoint_sha256'],'coarsening endpoint')
 req(boundary(live)==bd,'boundary triangulation changed')
 req({f for f in counts if all(D[v]==v for v in f)}==original_fixed,'branch fixed locus changed')
 for t in live:req(tuple(sorted(D[v] for v in t)) in live,'final deck not simplicial')
 # Independently check every endpoint vertex link; the star-replacement proof is
 # stronger than this final manifold check but both are retained.
 endpoint_vertices={v for t in live for v in t};bdv={v for t in bd for v in t}
 for v in endpoint_vertices:surface({tuple(w for w in t if w!=v) for t in live if v in t},disk=(v in bdv))
 ordered=list(map(tuple,p['ordered_endpoint']));req(len(ordered)==len(live) and {tuple(sorted(t)) for t in ordered}==live,'endpoint branching facets')
 ed={};lookup={frozenset(t):t for t in ordered}
 for t in ordered:
  for a,b in combinations(t,2):
   key=tuple(sorted((a,b)));req(key not in ed or ed[key]==(a,b),'inconsistent product branching');ed[key]=(a,b)
  req(tuple(D[v] for v in t)==lookup[frozenset(D[v] for v in t)],'product branching is not equivariant')
 for t in initial_ordered:
  for a,b in combinations(t,2):
   if a in protected and b in protected:req(ed.get(tuple(sorted((a,b))))==(a,b),'frozen boundary/witness branching changed')
 f=[sum(len(x)==k for x in counts) for k in (1,2,3,4)]
 return [tuple(tags[v] for v in t) for t in ordered],{'initial_tetrahedra':len(initial),'endpoint_f':f,'deck_pair_operations':len(p['steps']),'individual_ball_replacements':stats['ball_replacements'],'endpoint_vertex_links_checked':len(endpoint_vertices),'boundary_triangles_preserved':len(bd),'fixed_cells_preserved':len(original_fixed),'protected_vertices':len(protected),'full_relative_ball_checks':dict(stats),'kernel':'standard-library ball/disk and full-face-incidence checks; no Regina'}

def replay(cert,stages,deck,ac,chart):
 start=time.monotonic();stages={s:[tup(t) for t in ts] for s,ts in stages.items()};deck={tup(a):tup(b) for a,b in deck.items()}
 req(cert['schema']=='equivariant-relative-3D-coarsening-v1','coarsening schema')
 required=check_protection(cert,stages,deck,ac,chart);out={};new={}
 for side in ('N','E'):new[side],out[side]=replay_part(cert['parts'][side],stages[side],deck,required)
 return new,deck,{'schema':'relative-equivariant-coarsening-replay-v1','passed':True,'parts':out,'witness_tetrahedra_preserved':len(cert['protected_source_witnesses']),'seconds':time.monotonic()-start}

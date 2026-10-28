"""Cell-identity based complement presentations for the static full markings.

The returned induced complexes are deformation retracts of complements, not
compact triangulated exteriors. No knot classification follows from a group.
"""
from __future__ import annotations
from itertools import combinations
from collections import defaultdict,deque,Counter
from pathlib import Path
import hashlib,json,gzip,heapq,math

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(p):
 p=Path(p);b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 b=json.dumps(x,sort_keys=True,separators=(',',':')).encode()
 p.write_bytes(gzip.compress(b,mtime=0) if p.suffix=='.gz' else b+b'\n')
def free(w):
 s=[]
 for a in w:
  if s and s[-1]==-a:s.pop()
  else:s.append(a)
 return tuple(s)
def cyc(w):
 s=free(w)
 while len(s)>1 and s[0]==-s[-1]:s=s[1:-1]
 return s
def inv(w):return tuple(-a for a in reversed(w))
def closure(facets):return {tuple(t) for f in facets for k in range(1,len(f)+1) for t in combinations(sorted(f),k)}

def skeleton(x,label,quotient=False):
 fs=closure(x['facets']);S=closure(x[label]);vs={v for f in S for v in f}
 if not S<=fs:raise ValueError('marked surface not a subcomplex')
 bad=sorted(f for f in fs if set(f)<=vs and f not in S)
 if bad:raise ValueError('MARKING_NOT_FULL: explicit logged refinement required')
 invmap=x['deck'];rep=lambda f:min(tuple(f),tuple(sorted(invmap[v] for v in f))) if quotient else tuple(f)
 if quotient and any(rep(f) not in S for f in S):raise ValueError('marking not invariant')
 rows=[sorted({rep(f) for f in fs if len(f)==k+1 and not (set(f)&vs)}) for k in range(3)]
 vertex=lambda v:min(v,invmap[v]) if quotient else v
 vids=[f[0] for f in rows[0]];es=rows[1];idx={e:i for i,e in enumerate(es)}
 ends=[(vertex(a),vertex(b)) for a,b in es]
 if any(a==b for a,b in ends):raise ValueError('edge has repeated quotient vertex')
 def oriented_edge(a,b):
  r=rep(tuple(sorted((a,b))));j=idx[r]
  if (a,b)==r:return j+1
  if (b,a)==r:return -j-1
  im=(invmap[a],invmap[b])
  if quotient and im==r:return j+1
  if quotient and im[::-1]==r:return -j-1
  raise ValueError('unmatched oriented edge')
 triangles=[]
 for a,b,c in rows[2]:triangles.append([oriented_edge(a,b),oriented_edge(b,c),oriented_edge(c,a)])
 if quotient:
  for k,row in enumerate(rows):
   for f in row:
    if tuple(sorted(invmap[v] for v in f))==f:raise ValueError('not a free cover on complement')
 return {'vertices':vids,'edge_representatives':es,'edge_endpoints':ends,
         'triangle_representatives':rows[2],'triangle_boundaries':triangles,
         'fullness_verified':True,'refinement_steps':0,
         'complement_f':[len({rep(f) for f in fs if len(f)==k+1 and not set(f)&vs}) for k in range(5)]}

def presentation(sk,root=None):
 V=sk['vertices'];E=sk['edge_endpoints'];root=V[0] if root is None else root
 adj=defaultdict(list)
 for i,(a,b) in enumerate(E):adj[a].append((b,i));adj[b].append((a,i))
 paths={root:()};tree=set();todo=deque([root])
 while todo:
  a=todo.popleft()
  for b,i in sorted(adj[a]):
   if b in paths:continue
   sign=1 if E[i][0]==a else -1
   paths[b]=paths[a]+(sign*(i+1),);tree.add(i);todo.append(b)
 if set(paths)!=set(V):raise ValueError('disconnected complement')
 gen={e:i+1 for i,e in enumerate(e for e in range(len(E)) if e not in tree)}
 def word(path):return free((1 if e>0 else -1)*gen[abs(e)-1] for e in path if abs(e)-1 in gen)
 rel=sorted(set(filter(None,(cyc(word(t)) for t in sk['triangle_boundaries']))))
 return {'root':root,'tree_edges':sorted(tree),'generator_edge_ids':[e for e in gen],
         'generators':len(gen),'relators':[list(w) for w in rel]},paths,word

def eliminate(p,max_word=200000):
 rel={};occ=defaultdict(set);heap=[];n=0;active=set(range(1,p['generators']+1));events=[]
 def put(w):
  nonlocal n
  w=cyc(w)
  if not w:return
  if len(w)>max_word:raise RuntimeError('WORD_CAP')
  n+=1;rel[n]=w
  for a in set(map(abs,w)):occ[a].add(n)
  heapq.heappush(heap,(len(w),n))
 def erase(i):
  for a in set(map(abs,rel.pop(i))):occ[a].discard(i)
 for w in p['relators']:put(w)
 while heap:
  _,i=heapq.heappop(heap)
  if i not in rel:continue
  w=rel[i];cs=Counter(map(abs,w));single=[a for a in cs if cs[a]==1]
  if not single:continue
  g=min(single,key=lambda a:(len(occ[a]),a));j=next(j for j,v in enumerate(w) if abs(v)==g)
  tail=w[j+1:]+w[:j];sub=inv(tail) if w[j]>0 else tail
  events.append({'generator':g,'relator':list(w),'replacement':list(sub)})
  affected=sorted(occ[g]);erase(i);active.remove(g)
  for z in affected:
   if z==i:continue
   new=[]
   for a in rel[z]:new.extend(sub if a==g else inv(sub) if a==-g else (a,))
   erase(z);put(new)
 return {'eliminations':events,'generators':sorted(active),'relators':[list(w) for w in sorted(set(rel.values()))]}

def rewrite(w,cert):
 w=tuple(w)
 for e in cert['eliminations']:
  g=e['generator'];s=tuple(e['replacement']);z=[]
  for a in w:z.extend(s if a==g else inv(s) if a==-g else (a,))
  w=free(z) # NOT cyclic: this is a based element.
 return list(w)

def meridian_data(x,sk,downsk):
 tri=min(map(tuple,x['ramification_triangles']));link=set()
 for t in x['facets']:
  if set(tri)<=set(t):link.add(tuple(v for v in t if v not in tri))
 adj=defaultdict(list)
 for a,b in link:adj[a].append(b);adj[b].append(a)
 if not adj or any(len(a)!=2 for a in adj.values()):raise ValueError('normal link not circle')
 start=min(adj);path=[start,min(adj[start])]
 while path[-1]!=start:
  a,b=path[-2:];c=next(v for v in adj[b] if v!=a);path.append(c)
  if len(path)>len(adj)+1:raise ValueError('normal link disconnected')
 if len(path)!=len(adj)+1:raise ValueError('link length')
 deck=x['deck'];n=len(adj);h=n//2
 if n%2 or any(deck[path[i]]!=path[(i+h)%n] for i in range(n)):raise ValueError('deck not half-turn')
 qi={e:i+1 for i,e in enumerate(map(tuple,downsk['edge_representatives']))}
 ui={e:i+1 for i,e in enumerate(map(tuple,sk['edge_representatives']))}
 def qedge(a,b):
  r=min(tuple(sorted((a,b))),tuple(sorted((deck[a],deck[b]))));z=qi[r]
  return z if (a,b)==r or (deck[a],deck[b])==r else -z
 upedges=[ui[tuple(sorted((a,b)))]*(1 if a<b else -1) for a,b in zip(path,path[1:])]
 projected=[qedge(a,b) for a,b in zip(path,path[1:])]
 if projected[:h]!=projected[h:]:raise ValueError('meridian not a double traversal')
 # Vertices 0..h-1 give one base meridian; it lifts to an open path.
 eps=[]
 for a,b in downsk['edge_representatives']:eps.append(int(a!=min(a,deck[a]))^int(b!=min(b,deck[b])))
 for t in downsk['triangle_boundaries']:
  if sum(eps[abs(e)-1] for e in t)%2:raise ValueError('not a sheet cocycle')
 if sum(eps[abs(e)-1] for e in projected[:h])%2!=1:raise ValueError('meridian parity not one')
 return {'ramification_triangle':list(tri),'upstairs_normal_vertex_cycle':path,
         'upstairs_normal_edge_cycle':upedges,'downstairs_normal_edge_cycle':projected[:h],
         'projected_upstairs_normal_edge_cycle':projected,'upstairs_meridian_projects_to_square':True,
         'downstairs_sheet_cocycle_on_edges':eps,'downstairs_meridian_sheet_parity':1}

def build(name,fixtures):
 x=read(fixtures/f'{name}_cover4.json.gz');result={'schema':'static-branch-exterior-v2','fixture':name,'source_sha256':digest(x),'cases':{}}
 up=skeleton(x,'ramification_triangles');dn=skeleton(x,'ramification_triangles',True)
 mer=meridian_data(x,up,dn)
 for role,sk in [('ramification_upstairs',up),('branch_downstairs',dn),('twin_upstairs',skeleton(x,'twin_triangles'))]:
  p,paths,word=presentation(sk);c=eliminate(p)
  rec={'skeleton':sk,'presentation':p,'reduction':c,'scope':'fundamental group of complement deformation retract, not a knot classification'}
  if role!='twin_upstairs':
   cycle=mer['upstairs_normal_edge_cycle' if role=='ramification_upstairs' else 'downstairs_normal_edge_cycle']
   vertex=mer['upstairs_normal_vertex_cycle'][0]
   if role=='branch_downstairs':vertex=min(vertex,x['deck'][vertex])
   based=paths[vertex]+tuple(cycle)+inv(paths[vertex]);w=word(based)
   rec['meridian_raw_word']=list(w);rec['meridian_reduced_word']=rewrite(w,c)
  result['cases'][role]=rec
 # Verify quotient 2-skeleton is genuinely the two-sheet orbit cover.
 for k,key in [(0,'vertices'),(1,'edge_representatives'),(2,'triangle_representatives')]:
  if len(up[key])!=2*len(dn[key]):raise ValueError('not two sheets')
 # Record edge-path covering map; no claim of a peripheral group isomorphism from strings.
 deck=x['deck'];edge_idx={tuple(e):i+1 for i,e in enumerate(dn['edge_representatives'])}
 emap=[]
 for a,b in up['edge_representatives']:
  r=min((a,b),tuple(sorted((deck[a],deck[b]))));sgn=1 if (a,b)==r or (deck[a],deck[b])==r else -1
  emap.append(sgn*edge_idx[r])
 mer['upstairs_to_downstairs_oriented_edge_map']=emap
 canon=lambda w:min(tuple(v[i:]+v[:i]) for v in (list(w),list(inv(w))) for i in range(3))
 base_triangles={canon(t) for t in dn['triangle_boundaries']}
 for tr in up['triangle_boundaries']:
  images=[emap[abs(e)-1]*(1 if e>0 else -1) for e in tr]
  # Each lifted 2-cell maps to one base 2-cell, up to orientation and cyclic start.
  if canon(images) not in base_triangles:raise ValueError('cover face map')
 result['covering_meridian']=mer
 return result

def canonical_relator(w):
 w=cyc(w)
 if not w:return ()
 return min(v[i:]+v[:i] for v in (w,inv(w)) for i in range(len(v)))
def substitute(w,images):
 return free(a for g in w for a in (images[g] if g>0 else inv(images[-g])))

def amalgam_nf(word,p,q):
 """Exact normal form in <x,y | x^p=y^q>, central Z plus C_p*C_q syllables."""
 central=0;syll=[]
 for z in word:
  g=abs(z);n=1 if z>0 else -1
  if g not in (1,2):raise ValueError('unknown amalgam generator')
  if syll and syll[-1][0]==g:n+=syll.pop()[1]
  c,r=divmod(n,p if g==1 else q);central+=c
  if r:syll.append((g,r))
 return central,tuple(syll)

def identify(case,fixture,role):
 red=case['reduction'];gs=red['generators'];rr=red['relators'];ren={g:[i+1] for i,g in enumerate(gs)}
 rel=[substitute(w,ren) for w in rr]
 if fixture=='unknot':
  if len(gs)!=1 or rel:raise ValueError('expected infinite cyclic group not obtained')
  return {'group':'Z','H1':{'free_rank':1,'torsion':[]}}
 if role=='branch_downstairs':
  if len(gs)!=2:raise ValueError('downstairs generator count')
  # Original a=x*y^-2, b=y; inverse x=a*b^2, y=b.
  images={1:[1,-2,-2],2:[2]};inverse={1:[1,2,2],2:[2]}
  target=canonical_relator([1,1,-2,-2,-2])
  if {canonical_relator(substitute(w,images)) for w in rel}!={target}:raise ValueError('not certified trefoil relator')
  name='trefoil_group_<x,y|x^2=y^3>';H={'free_rank':1,'torsion':[]}
 elif role=='ramification_upstairs':
  if len(gs)!=2:raise ValueError('upstairs generator count')
  # Original a=v*u^-1, b=u; inverse u=b, v=a*b.
  images={1:[2,-1],2:[1]};inverse={1:[2],2:[1,2]}
  target=canonical_relator([1,1,1,-2,-2,-2])
  if {canonical_relator(substitute(w,images)) for w in rel}!={target}:raise ValueError('not certified index-two relator')
  name='index2_trefoil_<u,v|u^3=v^3>';H={'free_rank':1,'torsion':[3]}
 else:
  if len(gs)!=2:raise ValueError('twin generator count')
  canonical={canonical_relator(w) for w in rel}
  required={canonical_relator([2,2,2]),canonical_relator([1,2,-1,-2])}
  if not required<=canonical:raise ValueError('missing direct product relations')
  for w in rel:
   if sum((1 if a>0 else -1) for a in w if abs(a)==1)!=0 or sum((1 if a>0 else -1) for a in w if abs(a)==2)%3:raise ValueError('extra relation in direct product')
  return {'group':'Z x Z/3','H1':{'free_rank':1,'torsion':[3]},'proof':'contains commutator and cubic relation; all other relations vanish in this direct product'}
 for g in (1,2):
  if substitute(substitute([g],images),inverse)!=(g,) or substitute(substitute([g],inverse),images)!=(g,):raise ValueError('basis map not invertible')
 ans={'group':name,'H1':H,'old_basis_in_standard_basis':images,'standard_basis_in_old_basis':inverse}
 if 'meridian_reduced_word' in case:
  w=substitute(case['meridian_reduced_word'],ren);ans['meridian_in_standard_basis']=list(substitute(w,images))
 return ans

def augment(r):
 for role,c in r['cases'].items():c['identification']=identify(c,r['fixture'],role)
 if r['fixture']=='spun_trefoil':
  dn=r['cases']['branch_downstairs'];up=r['cases']['ramification_upstairs']
  pu,paths,wu=presentation(up['skeleton']);pd,pdpaths,wd=presentation(dn['skeleton'])
  emap=r['covering_meridian']['upstairs_to_downstairs_oriented_edge_map']
  # Images of surviving upstairs generators in the reduced downstairs presentation.
  if pu['root']!=pd['root']:raise ValueError('basepoints require an explicit whisker')
  image={}
  dg=dn['reduction']['generators'];rename={g:[i+1] for i,g in enumerate(dg)}
  std={int(k):v for k,v in dn['identification']['old_basis_in_standard_basis'].items()}
  for g in up['reduction']['generators']:
   e=pu['generator_edge_ids'][g-1];a,b=up['skeleton']['edge_endpoints'][e]
   loop=paths[a]+(e+1,)+inv(paths[b]);projected=[emap[abs(z)-1]*(1 if z>0 else -1) for z in loop]
   raw=wd(projected);v=rewrite(raw,dn['reduction']);image[g]=substitute(substitute(v,rename),std)
  for w in up['reduction']['relators']:
   if amalgam_nf(substitute(w,image),2,3)!=(0,()):raise ValueError('cover inclusion violates relation')
  mu=dn['identification']['meridian_in_standard_basis'];tilde=substitute(up['meridian_reduced_word'],image)
  if amalgam_nf(tilde+inv(tuple(mu)*2),2,3)!=(0,()):raise ValueError('based meridian square identity')
  r['covering_meridian']['remaining_upstairs_generators_in_trefoil_xy']={str(k):list(v) for k,v in image.items()}
  r['covering_meridian']['meridian_square_checked_in_exact_trefoil_normal_form']=True
 return r

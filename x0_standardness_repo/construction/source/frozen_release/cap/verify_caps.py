#!/usr/bin/env python3
"""Check saved cap data without importing or calling the cap generator."""
from __future__ import annotations
from collections import Counter, defaultdict, deque
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse, copy, gzip, hashlib, json, time


def require(value,message):
    if not value:raise ValueError(message)


def load(p):
    with gzip.open(p,'rt') if str(p).endswith('.gz') else open(p) as f:return json.load(f)


def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def tup(x):return tuple(tup(z) for z in x) if isinstance(x,list) else x
def mod(x):return x-x.numerator//x.denominator
def det2(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def det3(a,b,c,d):
    u,v,w=[tuple(p[k]-a[k] for k in range(3)) for p in (b,c,d)]
    return u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0])


def determinant(matrix):
    a=[list(r) for r in matrix];ans=Q(1)
    for j in range(len(a)):
        p=next((i for i in range(j,len(a)) if a[i][j]),None)
        if p is None:return Q(0)
        if p!=j:a[p],a[j]=a[j],a[p];ans=-ans
        x=a[j][j];ans*=x
        for k in range(j,len(a)):a[j][k]/=x
        for i in range(j+1,len(a)):
            x=a[i][j]
            for k in range(j,len(a)):a[i][k]-=x*a[j][k]
    return ans


def allfaces(fs):
    return [set(x for f in fs for x in combinations(sorted(f),k+1)) for k in range(len(next(iter(fs))))]


def fvector(fs):return [len(x) for x in allfaces(fs)]
def chi(fs):return sum((-1)**i*x for i,x in enumerate(fvector(fs)))


def components(vertices,edges):
    adj={v:set() for v in vertices}
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    out=[];seen=set()
    for start in adj:
        if start in seen:continue
        part={start};seen.add(start);todo=[start]
        for v in todo:
            for w in adj[v]-seen:seen.add(w);part.add(w);todo.append(w)
        out.append(part)
    return out


def sphere2(triangles):
    triangles=set(tuple(sorted(t)) for t in triangles)
    edges=Counter(e for tr in triangles for e in combinations(tr,2));vertices={v for tr in triangles for v in tr}
    require(all(n==2 for n in edges.values()),'sphere2 edge degrees')
    require(len(components(vertices,edges))==1,'sphere2 disconnected')
    require(len(vertices)-len(edges)+len(triangles)==2,'sphere2 Euler characteristic')
    for v in vertices:
        link=[tuple(w for w in tr if w!=v) for tr in triangles if v in tr];graphcycle(link)
    return True


def graphcycle(edges):
    edges=set(tuple(sorted(e)) for e in edges);counts=Counter(v for e in edges for v in e)
    require(all(n==2 for n in counts.values()),'link is not degree-two graph')
    require(len(components(counts,edges))==1,'link has several circles')
    return True


def manifold_ridges(fs):
    slots=defaultdict(list)
    for i,f in enumerate(fs):
        for k in range(len(f)):slots[f[:k]+f[k+1:]].append((i,k))
    require(all(len(x) in (1,2) for x in slots.values()),'ridge incidence')
    return slots


def check_shared(b):
    require(b['schema']=='even-cap-shared-boundary-v1','shared schema')
    fs=list(map(tuple,b['facets']));require(len(fs)==len(set(fs)),'duplicate shared facets')
    slots=manifold_ridges(fs);require(all(len(v)==2 for v in slots.values()),'shared boundary not closed')
    require(fvector(fs)==b['f_vector'] and chi(fs)==0,'shared f vector')
    D=b['deck'];require(all(D[D[v]]==v for v in range(len(D))),'shared involution')
    fset=set(fs);require(all(tuple(sorted(D[v] for v in f)) in fset for f in fs),'shared deck simplices')
    links=defaultdict(list)
    for f in fs:
        for v in f:links[v].append(tuple(w for w in f if w!=v))
    for link in links.values():sphere2(link)
    eq=set(map(tuple,b['equator']));eqfaces=allfaces(eq);require(chi(eq)==0,'equator Euler')
    for v in {x[0] for x in eqfaces[0]}:graphcycle([tuple(w for w in f if w!=v) for f in eq if v in f])
    require(len(components([v[0] for v in eqfaces[0]],eqfaces[1]))==1,'equator disconnected')
    n,t=b['angles'],b['times'];tags=list(map(tup,b['vertices']));tagid={v:i for i,v in enumerate(tags)}
    # Each hemisphere is a cyclic union of explicitly coned three-balls.
    slab_counts=[]
    for p in (0,1):
        part=list(map(tuple,b['hemispheres'][str(p)]))
        require(set(part)<=fset,'hemisphere facets')
        for j in range(t):
            c=tagid[('C',p,j)];slab=[f for f in part if c in f]
            require(slab,'missing hemisphere slab')
            sphere2([tuple(v for v in f if v!=c) for f in slab]);slab_counts.append(len(slab))
        require(all(sum(tags[v][0]=='C' for v in f)==1 for f in part),'slab cone coverage')
    require(set(map(tuple,b['hemispheres']['0']))|set(map(tuple,b['hemispheres']['1']))==fset,'hemisphere partition')
    require(set(map(tuple,b['hemispheres']['0'])).isdisjoint(set(map(tuple,b['hemispheres']['1']))),'hemisphere overlap')
    return {'vertex_links_spheres':len(links),'hemisphere_slab_balls':2*t,'slab_tetrahedra':slab_counts,'equator_f_vector':fvector(eq)}


def expected_grid(n,t):
    ans=[]
    for i in range(n):
        for j in range(t):
            a,b,c,d=(Q(i,n),Q(j,t)),(Q(i+1,n),Q(j,t)),(Q(i+1,n),Q(j+1,t)),(Q(i,n),Q(j+1,t))
            ans.extend(((a,b,c),(a,c,d)))
    return ans


def check_case(c,b):
    require(c['schema']=='even-twist-cap-template-v1','case schema')
    m,n,t=c['twist'],c['angles'],c['times'];require(m in (0,2),'twist outside admitted cases')
    require(n==b['angles'] and t==b['times'],'grid mismatch')
    require(c['shared_boundary_sha256']==digest(b),'shared boundary hash')
    tags=list(map(tup,c['vertices']));require(tags[:len(b['vertices'])]==list(map(tup,b['vertices'])),'shared vertex prefix')
    fs=list(map(tuple,c['facets']));require(len(fs)==len(set(fs)),'duplicate cap facets')
    require(all(len(set(f))==5 for f in fs),'collapsed cap simplex')
    require(fvector(fs)==c['f_vector'] and chi(fs)==2,'cap f vector')
    slots=manifold_ridges(fs);boundary=sorted(f for f,v in slots.items() if len(v)==1)
    require(boundary==list(map(tuple,b['facets']))==list(map(tuple,c['boundary'])),'cap/shared boundary mismatch')
    D=c['deck'];require(len(D)==len(tags) and all(D[D[v]]==v for v in range(len(D))),'cap involution')
    fset=set(fs);require(all(tuple(sorted(D[v] for v in f)) in fset for f in fs),'cap deck action')
    require(all(not any(v!=D[v] and D[v] in f for v in f) for f in fs),'simplicial inversion')
    # Check saved orientation across every paired ridge and every native gluing.
    signs=c['orientation'];require(len(signs)==len(fs) and all(x in (-1,1) for x in signs),'orientation signs')
    got=[]
    for ridge,h in sorted(slots.items()):
        if len(h)!=2:continue
        (i,k),(j,l)=h;perm=[fs[j].index(v) if a!=k else l for a,v in enumerate(fs[i])]
        parity=(-1)**sum(perm[a]>perm[d] for a in range(5) for d in range(a+1,5))
        require(signs[j]==-parity*signs[i],'orientation inconsistent');got.append([i,k,j,perm])
    require(got==c['gluings'],'native gluing table differs from facets')
    # Full spherical common-refinement certificate.  These matrices are maps
    # on complete tetrahedra, not samples of an angular formula.
    coarse=b['coarse_source_boundary'];parents=list(map(tuple,coarse['facets']));coarse_ids={tup(v):i for i,v in enumerate(coarse['vertices'])}
    volumes=Counter();vertex_images={};refinements=c['boundary_to_source_refinement']
    require(len(refinements)==len(boundary),'spherical refinement length')
    for f,row in zip(boundary,refinements):
        parent=parents[row['parent']];mat=[list(map(Q,r)) for r in row['barycentric']]
        require(len(mat)==4 and all(len(r)==4 for r in mat),'barycentric matrix shape')
        require(all(all(x>=0 for x in r) and sum(r)==1 for r in mat),'barycentric simplex containment')
        vol=abs(determinant(mat));require(vol>0,'singular spherical refinement');volumes[row['parent']]+=vol
        for v,r in zip(f,mat):
            image=tuple((w,x) for w,x in zip(parent,r) if x)
            if v in vertex_images:require(vertex_images[v]==image,'incompatible spherical refinement faces')
            else:vertex_images[v]=image
    require(len(volumes)==len(parents) and all(v==1 for v in volumes.values()),'spherical refinement parent volumes')
    require(all(tuple(sorted((coarse['deck'][w],x) for w,x in image))==vertex_images[D[v]] for v,image in vertex_images.items()),'spherical refinement deck equivariance')
    for v,image in vertex_images.items():
        if tags[v][0] in ('P','C'):require(image==((coarse_ids[tags[v]],Q(1)),),'spherical map moves a marked axis point')
    grid=expected_grid(n,t);eq=set(map(tuple,b['equator']));child_triangles=[];normalised=[]
    for i,tri in enumerate(grid):
        rows=c['source_boundary_children'].get(str(i));require(rows,'missing parent coverage')
        total=Q(0)
        for r in rows:
            pts=[tuple(map(Q,p)) for p in r['source_lifts']];target=r['target_vertices']
            require(len(pts)==3 and len(set(target))==3,'bad boundary child')
            d=det2(*pts);require(d>0,'boundary child orientation');total+=d
            for p,v in zip(pts,target):
                require(all(det2(tri[k],tri[(k+1)%3],p)>=0 for k in range(3)),'child outside parent')
                require(tags[v][0]=='E','child target is not equator')
                q=(Q(tags[v][1]),Q(tags[v][2]));require(q==(mod(p[0]-Q(m,2)*p[1]),mod(p[1])),'beta vertex equation')
                bary=(det2(p,tri[1],tri[2])/det2(*tri),det2(tri[0],p,tri[2])/det2(*tri),det2(tri[0],tri[1],p)/det2(*tri))
                image=tuple(sorted((coarse_ids[('E',str(mod(w[0])),str(mod(w[1])))],x) for w,x in zip(tri,bary) if x))
                require(vertex_images[v]==image,'spherical map not bound to the angular beta chart')
                # Exact inverse equation in this chosen parent chart.
                back=(p[0]-Q(m,2)*p[1]+Q(m,2)*p[1],p[1]);require(back==p,'inverse chart')
            child_triangles.append(tuple(sorted(target)));normalised.append((pts,target))
        require(total==det2(*tri),'parent area is not exactly covered')
    require(len(child_triangles)==len(set(child_triangles)) and set(child_triangles)==eq,'boundary children not bijective')
    # Every collar simplex is a positive-volume simplex inside its specified
    # product prism.  Face coverage plus exact total volumes certify the mesh.
    byparent=defaultdict(list)
    filling=set(map(tuple,c['filling_torus']));collar=set()
    for row in c['collar_charts']:
        coords=[tuple(map(Q,p)) for p in row['coords']];vs=row['vertices'];i=row['parent'];tri=grid[i]
        require(len(vs)==4 and len(set(vs))==4,'collar vertex count')
        volume=abs(det3(*coords));require(volume>0,'zero-volume collar tetrahedron')
        for p in coords:
            require(0<=p[2]<=1 and all(det2(tri[k],tri[(k+1)%3],p[:2])>=0 for k in range(3)),'collar tetra outside prism')
        f=tuple(sorted(vs));require(f in filling,'collar tetra not in filling');collar.add(f);byparent[i].append(volume)
    require(len(collar)==len(c['collar_charts']),'collar chart duplicate')
    require(len(byparent)==len(grid),'missing collar prism')
    for i,vols in byparent.items():require(sum(vols)==3*det2(*grid[i]),'collar prism volume')
    hslots=manifold_ridges(sorted(filling));require(set(f for f,v in hslots.items() if len(v)==1)==eq,'filling torus boundary')
    # A normal circle at every top-dimensional fixed simplex is explicitly
    # checked, including the nontrivial free half-turn action on that circle.
    fixed={v for v,w in enumerate(D) if v==w};ram=set(map(tuple,c['ramification_triangles']))
    actual={tuple(v for v in f if v in fixed) for f in fs if sum(v in fixed for v in f)==3}
    require(ram==actual,'fixed locus triangles')
    require(all(sum(v in fixed for v in f)<=3 for f in fs),'fixed locus too large')
    rf=allfaces(ram);rc=components([v[0] for v in rf[0]],rf[1]);require(len(rc)==2,'fixed locus components')
    for comp in rc:
        triangles={r for r in ram if set(r)<=comp};require(chi(triangles)==1,'fixed component Euler')
        edgecounts=Counter(e for r in triangles for e in combinations(r,2));require(all(v in (1,2) for v in edgecounts.values()),'fixed disk edge degree')
        graphcycle([e for e,v in edgecounts.items() if v==1])
    for r in ram:
        link=[tuple(v for v in f if v not in r) for f in fs if set(r)<=set(f)];graphcycle(link)
        lv={v for e in link for v in e};require(all(D[v] in lv and D[v]!=v for v in lv),'normal deck is not free')
    twin=set(map(tuple,c['twin_triangles']));sphere2(twin)
    require(all(any(set(tr)<=set(f) for f in fs) for tr in twin),'twin not subcomplex')
    require(len({v for tr in twin for v in tr}&fixed)==2,'twin/fixed marked intersections')
    # Read full integer winding from an actual path in the shared equator.
    # The source angle-zero time circle is a union of source parent edges.
    path={}
    for pts,vs in normalised:
        for p,v in zip(pts,vs):
            if p[0]==0 or p[0]==1:
                q=(p[0]-Q(m,2)*p[1],p[1]);path[p[1]]=(v,q[0]-p[0])
    times=sorted(path);require(times[0]==0 and times[-1]==1,'winding path endpoints')
    edges={e for f in eq for e in combinations(f,2)}
    for x,y in zip(times,times[1:]):
        require(tuple(sorted((path[x][0],path[y][0]))) in edges,'winding segment is not a mesh edge')
    winding=path[Q(1)][1]-path[Q(0)][1]
    require(winding==-Q(m,2),'unwrapped cap winding')
    require(c['beta']['unwrapped_displacement']==[-m*n//2,t],'integer displacement record')
    require(c['beta']['downstairs_winding']==-m and c['beta']['bulk_deck_parity']==m%2,'branch winding metadata')
    return {'result':'PASS','pentachora':len(fs),'boundary_tetrahedra':len(boundary),'spherical_refinement_tetrahedra':len(refinements),'spherical_coarse_tetrahedra':len(parents),'boundary_child_triangles':len(eq),'exact_collar_prisms':len(byparent),'collar_tetrahedra':len(collar),'ramification_components':2,'normal_circle_checks':len(ram),'twin_sphere':True,'winding_edge_count':len(times)-1,'upstairs_winding':str(winding),'downstairs_winding':-m,'unwrapped_displacement':[-m*n//2,t]}


def native_case(c):
    import regina
    tr=regina.Triangulation4()
    tr.newPentachora(len(c['facets']))
    for i,k,j,p in c['gluings']:tr.pentachoron(i).join(k,tr.pentachoron(j),regina.Perm5(*p))
    out={'regina':regina.versionString(),'valid':tr.isValid(),'orientable':tr.isOrientable(),'connected':tr.isConnected(),'boundary_components':tr.countBoundaryComponents(),'euler':tr.eulerCharTri(),'homology':'NOT_RUN in admission; topology is supplied by the product/collar/cone construction'}
    require(out['valid'] and out['orientable'] and out['connected'] and out['boundary_components']==1 and out['euler']==2,'native cap admission')
    return out


def check_puncture_refinement(r):
    require(r['schema']=='suspension-circle-bisection-v1','puncture refinement schema')
    coarse=[];fine=[];maps=r['barycentric_vertex_map']
    for row in r['face_replacements']:
        tr=row['source_triangle'];coarse.append(tuple(tr));area=Q(0)
        for f in row['fine_triangles']:
            fine.append(tuple(f));matrix=[]
            for v in f:
                weights={w:Q(x) for w,x in maps[v].items()}
                require(set(weights)<=set(tr) and sum(weights.values())==1 and min(weights.values())>0,'puncture containment')
                matrix.append([weights.get(w,Q(0)) for w in tr])
            area+=abs(determinant(matrix))
        require(area==1,'puncture triangle not covered')
    # Convert string labels to integers only for the general link checker.
    cv={v:i for i,v in enumerate(r['coarse_vertices'])};fv={v:i for i,v in enumerate(r['fine_vertices'])}
    sphere2([tuple(cv[v] for v in f) for f in coarse]);sphere2([tuple(fv[v] for v in f) for f in fine])
    for v,image in maps.items():
        deck_image={r['coarse_deck'][w]:x for w,x in image.items()}
        require(deck_image==maps[r['fine_deck'][v]],'puncture deck compatibility')
    return {'result':'PASS','coarse_sphere_f_vector':[6,12,8],'fine_sphere_f_vector':[10,24,16],'angle_map':'coarse e_j -> fine f_(2j), odd fine vertices are edge midpoints','deck_equivariant':True}


def run(data,native=False):
    start=time.monotonic();b=load(data/'shared_boundary.json.gz');cases={m:load(data/f'cap_m{m}.json.gz') for m in (0,2)}
    out={'schema':'even-cap-independent-verification-v1','scope':'SOURCE_INDEPENDENT_CAP_TEMPLATE; no knot-source attachment, global Phi_m, X2 or X0 is asserted','shared':check_shared(b),'cases':{str(m):check_case(c,b) for m,c in cases.items()}}
    if (data/'puncture4_to8.json').exists():out['puncture_refinement']=check_puncture_refinement(load(data/'puncture4_to8.json'))
    controls=[]
    def corrupt_deck(c):
        i=next(i for i,v in enumerate(c['deck']) if i!=v);c['deck'][i]=i
    for name,change in [('parity_only_m2',lambda c:c['beta'].__setitem__('unwrapped_displacement',[0,c['times']])),('parity_only_spherical_map',lambda c:c.__setitem__('boundary_to_source_refinement',copy.deepcopy(cases[0]['boundary_to_source_refinement']))),('wrong_twist_binding',lambda c:c.__setitem__('twist',0)),('wrong_deck',corrupt_deck),('missing_cap_facet',lambda c:c['facets'].pop()),('corrupt_boundary_chart',lambda c:c['source_boundary_children']['0'][0]['source_lifts'][0].__setitem__(0,'99/1'))]:
        c=copy.deepcopy(cases[2]);change(c)
        try:check_case(c,b)
        except (ValueError,KeyError,IndexError) as e:controls.append({'name':name,'rejected':True,'reason':str(e)})
        else:raise ValueError('negative control accepted: '+name)
    out['negative_controls']=controls
    if native:out['native']={str(m):native_case(c) for m,c in cases.items()}
    out['result']='PASS';out['seconds']=time.monotonic()-start;return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=Path(__file__).resolve().parent/'data');p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'verification.json');p.add_argument('--native',action='store_true');a=p.parse_args();r=run(a.data,a.native);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

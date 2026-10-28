#!/usr/bin/env python3
"""Independent saved-data check for the source-to-cap boundary collar."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse,copy,hashlib,json,time
from verify_caps import load,digest,tup,require,mod,determinant,fvector,chi,manifold_ridges,sphere2,graphcycle,components


def circle(theta):
    z=mod(theta)*4;i=int(z);s=z-i;corners=((1,0),(0,1),(-1,0),(0,-1))
    return tuple((1-s)*corners[i][k]+s*corners[(i+1)%4][k] for k in (0,1))


def expected_source(interface):
    T=interface['time_intervals'];vertices=sorted((v,j) for v in interface['poles']+interface['equator_cycle_diagram_positive'] for j in range(T));vi={v:i for i,v in enumerate(vertices)};fs=[]
    for a,b,c in interface['ordered_triangles']:
        for j in range(T):
            k=(j+1)%T
            for tet in (((a,j),(b,j),(c,j),(c,k)),((a,j),(b,j),(b,k),(c,k)),((a,j),(a,k),(b,k),(c,k))):fs.append(tuple(sorted(vi[v] for v in tet)))
    return vertices,sorted(fs)


def geometry(coarse,interface,cap=False):
    T=interface['time_intervals'];out={}
    if not cap:
        eq=interface['equator_cycle_diagram_positive'];poles=interface['poles'];xy={v:circle(Q(i,4)) for i,v in enumerate(eq)};xy.update({p:(Q(0),Q(0)) for p in poles})
        tags=list(map(tuple,coarse['vertices']))
        for i,f in enumerate(coarse['facets']):
            p=next(v for v,j in (tags[k] for k in f) if v in poles);hem=poles.index(p);js={tags[k][1] for k in f};j=next(j for j in js if (j+1)%T in js)
            points=tuple((*xy[tags[k][0]],Q(tags[k][1] if not (j==T-1 and tags[k][1]==0) else T,T)) for k in f);out[i]=(hem,j,points)
    else:
        tags=list(map(tup,coarse['vertices']))
        for i,f in enumerate(coarse['facets']):
            _,hem,j=next(tags[v] for v in f if tags[v][0]=='C');points=[]
            for v in f:
                x=tags[v]
                if x[0]=='E':
                    t=Q(x[2]);t=t+1 if j==T-1 and t==0 else t;points.append((*circle(Q(x[1])),t))
                elif x[0]=='P':
                    t=Q(x[2],T);t=t+1 if j==T-1 and t==0 else t;points.append((Q(0),Q(0),t))
                else:points.append((Q(0),Q(0),Q(2*j+1,2*T)))
            out[i]=(hem,j,tuple(points))
    return out


def check(b,cap,interface,cut):
    require(b['schema']=='source-cap-boundary-bridge-v1','bridge schema')
    require(b['source_interface']==interface,'source interface differs')
    require(b['cap_shared_canonical_sha256']==digest(cap),'cap boundary digest')
    require(b['source_cut_binding'] is not None and b['source_cut_binding']['puncture_canonical_sha256']==digest(cut['puncture']),'puncture source binding')
    poles=interface['poles'];eq=interface['equator_cycle_diagram_positive'];T=interface['time_intervals']
    require(poles==cut['puncture']['fixed_edge_cut_vertices'],'source pole binding')
    require(eq==[cut['puncture']['equator_cycle'][i] for i in (0,3,2,1)],'diagram-positive angle reversal binding')
    require({tuple(sorted(t)) for t in interface['ordered_triangles']}==set(map(tuple,cut['puncture']['boundary_triangles'])),'source boundary triangles')
    require(all(cut['deck'][v]==interface['deck'][str(v)] for v in poles+eq),'source deck binding')
    require(all(t==sorted(t,key=lambda v:min(v,cut['deck'][v])) for t in interface['ordered_triangles']),'source branching order')
    src=b['source_boundary'];vs,sourcefs=expected_source(interface)
    require(vs==list(map(tuple,src['vertices'])) and sourcefs==list(map(tuple,src['facets'])),'source staircase product')
    fine=b['common_boundary'];ffs=list(map(tuple,fine['facets']));require(fvector(ffs)==fine['f_vector'],'common f-vector')
    require(len(ffs)==len(set(ffs)),'common facet duplicate')
    ridges=manifold_ridges(ffs);require(all(len(v)==2 for v in ridges.values()),'common boundary not closed')
    links=defaultdict(list)
    for f in ffs:
        for v in f:links[v].append(tuple(w for w in f if w!=v))
    for l in links.values():sphere2(l)
    require(chi(ffs)==0,'common boundary Euler')
    fd=fine['deck'];require(all(fd[fd[v]]==v for v in range(len(fd))),'common involution')
    fset=set(ffs);require(all(tuple(sorted(fd[v] for v in f)) in fset for f in ffs),'common deck facets')
    coarse={'source':src,'cap':cap};geos={s:geometry(c,interface,s=='cap') for s,c in coarse.items()};images={}
    for side,co in coarse.items():
        total=Counter();img={}
        for f,row in zip(ffs,fine['refinements']):
            require(tuple(row['vertices'])==f,'common refinement facet index')
            parent=row[side+'_parent'];mat=[list(map(Q,r)) for r in row[side+'_barycentric']];pf=co['facets'][parent]
            require(len(mat)==4 and all(len(r)==4 and sum(r)==1 and min(r)>=0 for r in mat),'common barycentric containment')
            d=abs(determinant(mat));require(d>0,'common refinement singular');total[parent]+=d
            hem,slab,coords=geos[side][parent]
            require(geos['source'][row['source_parent']][:2]==geos['cap'][row['cap_parent']][:2],'common hemisphere/time mismatch')
            for v,w in zip(f,mat):
                im=tuple((q,x) for q,x in zip(pf,w) if x)
                if v in img:require(img[v]==im,'common face-map compatibility')
                img[v]=im
                p=tuple(sum(w[j]*coords[j][k] for j in range(4)) for k in range(3));tag=fine['vertices'][v]
                require((Q(tag[1]),Q(tag[2]),Q(tag[3]))==(p[0],p[1],mod(p[2])),'common geometric chart binding')
                require(tag[0]==(-1 if abs(p[0])+abs(p[1])==1 else hem),'common hemisphere vertex binding')
        require(len(total)==len(co['facets']) and all(v==1 for v in total.values()),'common parent volume coverage')
        require(all(tuple(sorted((co['deck'][v],x) for v,x in im))==img[fd[k]] for k,im in img.items()),'common map deck equivariance')
        images[side]=img
    tags=list(map(tup,b['vertices']));fs=list(map(tuple,b['facets']));require(len(fs)==len(set(fs)),'bridge pentachoron duplicate')
    require(fvector(fs)==b['f_vector'] and chi(fs)==0,'bridge f-vector/Euler')
    slots=manifold_ridges(fs);bd=sorted(f for f,h in slots.items() if len(h)==1)
    D=b['deck'];require(all(D[D[v]]==v for v in range(len(D))),'bridge involution');full=set(fs)
    require(all(tuple(sorted(D[v] for v in f)) in full for f in fs),'bridge deck simplices')
    require(not any(any(D[v]!=v and D[v] in f for v in f) for f in fs),'bridge simplex inversion')
    charts_all=set();portfacets=[];statistics={}
    for side,co in coarse.items():
        port=b['ports'][side];pvs=port['vertices'];expected=sorted(tuple(sorted(pvs[v] for v in f)) for f in co['facets'])
        require(expected==list(map(tuple,port['facets'])),'boundary port facets')
        require(all(tags[v]==('B',side,i) for i,v in enumerate(pvs)),'boundary port vertex labels');portfacets.extend(expected)
        count=Counter();got=set()
        for row in b['collar_charts'][side]:
            parent=co['facets'][row['parent']];vs=row['vertices'];pts=[list(map(Q,p)) for p in row['barycentric_height']]
            require(len(vs)==5 and len(set(vs))==5 and all(len(p)==5 and sum(p[:4])==1 and min(p[:4])>=0 and 0<=p[4]<=1 for p in pts),'collar containment')
            for v,p in zip(vs,pts):
                tag=tags[v]
                if tag[0]=='B':require(tag[1]==side and p==[Q(w==tag[2]) for w in parent]+[Q(0)],'bottom collar chart binding')
                elif tag[0]=='R':require(p==[dict(images[side][tag[1]]).get(w,Q(0)) for w in parent]+[Q(1)],'top collar chart binding')
                else:raise ValueError('unexpected vertex in pulling collar')
            matrix=[[pts[j][k]-pts[0][k] for k in (1,2,3,4)] for j in range(1,5)];v=abs(determinant(matrix));require(v>0,'zero-volume collar simplex');count[row['parent']]+=v
            f=tuple(sorted(vs));require(f in full,'collar simplex missing');require(f not in got,'collar simplex repeated');got.add(f)
        require(len(count)==len(co['facets']) and all(v==4 for v in count.values()),'collar prism volume coverage')
        require(got==set(map(tuple,b['parts'][side])),'collar part binding');require(not(got&charts_all),'collar interiors overlap combinatorially');charts_all|=got
        statistics[side]={'coarse_tetrahedra':len(co['facets']),'pentachora':len(got),'prism_volumes_verified':len(count)}
    require(charts_all==full,'unaccounted bridge simplex');require(sorted(portfacets)==bd,'bridge boundary mismatch')
    # Both ports are individually connected and closed, hence these are the
    # two complete boundary components.
    for side in ('source','cap'):
        pf=list(map(tuple,b['ports'][side]['facets']));pr=manifold_ridges(pf);require(all(len(x)==2 for x in pr.values()),'port not closed')
        vertices={v for f in pf for v in f};edges={e for f in pf for e in combinations(f,2)};require(len(components(vertices,edges))==1,'port disconnected')
    fixed={v for v in range(len(D)) if v==D[v]};ram=set(map(tuple,b['ramification_triangles']));actual={tuple(v for v in f if v in fixed) for f in fs if sum(v in fixed for v in f)==3}
    require(ram==actual and all(sum(v in fixed for v in f)<=3 for f in fs),'bridge fixed triangles')
    rf=[{x for r in ram for x in combinations(r,k+1)} for k in range(3)];comps=components([x[0] for x in rf[0]],rf[1]);require(len(comps)==2,'bridge fixed components')
    for comp in comps:
        tr={r for r in ram if set(r)<=comp};require(chi(tr)==0,'fixed annulus Euler');ec=Counter(e for r in tr for e in combinations(r,2));require(all(n in (1,2) for n in ec.values()),'fixed annulus edge degrees');be={e for e,n in ec.items() if n==1};bc=components({v for e in be for v in e},be);require(len(bc)==2,'fixed annulus boundary circles')
        for q in bc:graphcycle([e for e in be if set(e)<=q])
    for tr in ram:
        lk=[tuple(v for v in f if v not in tr) for f in fs if set(tr)<=set(f)];graphcycle(lk);vv={v for e in lk for v in e};require(all(D[v] in vv and D[v]!=v for v in vv),'fixed normal half-turn')
    return {'result':'PASS','pentachora':len(fs),'f_vector':b['f_vector'],'common_vertex_links_spheres':len(links),'common_tetrahedra':len(ffs),'source_angle_cycle':eq,'source_poles':poles,'source_order':'min(v,deck[v]); poles last','ports':{s:len(b['ports'][s]['facets']) for s in ('source','cap')},'collars':statistics,'fixed_annuli':2,'normal_circle_checks':len(ram)}


def native(b):
    import regina
    tr=regina.Triangulation4();tr.newPentachora(len(b['facets']))
    # Reconstruct from actual facets, then compare with saved native gluings.
    fs=list(map(tuple,b['facets']));slots=manifold_ridges(fs);glu=[]
    for ridge,h in sorted(slots.items()):
        if len(h)!=2:continue
        (i,k),(j,l)=h;p=[fs[j].index(v) if q!=k else l for q,v in enumerate(fs[i])];glu.append([i,k,j,p]);tr.pentachoron(i).join(k,tr.pentachoron(j),regina.Perm5(*p))
    require(glu==b['gluings'],'bridge native gluing table')
    r={'regina':regina.versionString(),'valid':tr.isValid(),'orientable':tr.isOrientable(),'connected':tr.isConnected(),'boundary_components':tr.countBoundaryComponents(),'euler':tr.eulerCharTri()};require(r['valid'] and r['orientable'] and r['connected'] and r['boundary_components']==2 and r['euler']==0,'native bridge admission');return r


if __name__=='__main__':
    root=Path(__file__).resolve().parent;p=argparse.ArgumentParser();p.add_argument('--bridge',type=Path,default=root/'data/source_cap_bridge.json.gz');p.add_argument('--cap',type=Path,default=root/'data/shared_boundary.json.gz');p.add_argument('--interface',type=Path,default=root.parent/'peripheral/puncture_product_interface.json');p.add_argument('--source-cut',type=Path,default=root.parent/'peripheral/cut_source.json.gz');p.add_argument('--output',type=Path,default=root/'bridge_verification.json');p.add_argument('--native',action='store_true');a=p.parse_args();start=time.monotonic();b=load(a.bridge);c=load(a.cap);i=load(a.interface);cut=load(a.source_cut)
    require(b['source_interface_file_sha256']==hashlib.sha256(a.interface.read_bytes()).hexdigest(),'interface byte hash');require(b['source_cut_binding']['file_sha256']==hashlib.sha256(a.source_cut.read_bytes()).hexdigest(),'source cut byte hash')
    result=check(b,c,i,cut);result['scope']='source-bound S2 x S1 refinement bridge; no closed X2 or X0';result['negative_controls']=[]
    for name,mutate in [('reverse_meridian_again',lambda x:x['source_interface']['equator_cycle_diagram_positive'].reverse()),('wrong_common_map',lambda x:x['common_boundary']['refinements'][0]['source_barycentric'][0].__setitem__(0,'7')),('missing_bridge_pentachoron',lambda x:x['facets'].pop())]:
        x=copy.deepcopy(b);mutate(x)
        try:check(x,c,i,cut)
        except (ValueError,IndexError,KeyError) as e:result['negative_controls'].append({'name':name,'rejected':True,'reason':str(e)})
        else:raise ValueError('accepted corrupt bridge: '+name)
    if a.native:result['native']=native(b)
    result['seconds']=time.monotonic()-start;a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

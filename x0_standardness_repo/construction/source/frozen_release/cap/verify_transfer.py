#!/usr/bin/env python3
"""Independent affine-cell verification of saved F2 and Phi collar maps."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict,Counter
from itertools import combinations
import argparse,copy,hashlib,json,time
from verify_caps import load,require,determinant,sphere2
from verify_bridge import geometry


def mat(x):return [list(map(Q,r)) for r in x]
def vec(x):return tuple(map(Q,x))
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def mul(a,b):return [[sum(x*y for x,y in zip(r,col)) for col in zip(*b)] for r in a]
def image(p,a):return tuple(sum(x*y for x,y in zip(p,col)) for col in zip(*a))


def matrix_rank(rows):
    rows=[list(r) for r in rows];rank=0
    for k in range(len(rows[0]) if rows else 0):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][k]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank];q=rows[rank][k];rows[rank]=[x/q for x in rows[rank]]
        for i in range(rank+1,len(rows)):
            q=rows[i][k];rows[i]=[x-q*y for x,y in zip(rows[i],rows[rank])]
        rank+=1
    return rank


def arank(points):return matrix_rank([[x-y for x,y in zip(p,points[0])] for p in points[1:]]) if points else -1


def check(data,shared,cap,bridge):
    require(data['schema']=='finite-F2-Phi-polyhedral-transfer-v1','transfer schema')
    require(data['lifts']['constant_lattice_adjustment']==1 and data['lifts']['bulk_D2_unwrapped_shift']==1,'full integer lift binding')
    require(data['bindings']['source_interface']==bridge['source_interface'],'transfer source interface')
    interface=bridge['source_interface'];g=data['gamma'];p=data['Phi'];gs=list(map(vec,g['source_points']));gt=list(map(vec,g['target_points']));gkeys=list(map(tuple,g['vertex_keys']));gi={k:i for i,k in enumerate(gkeys)};gd=g['deck']
    require(len(g['cells'])==len(shared['facets']),'gamma cell coverage')
    coarse=shared['coarse_source_boundary'];cocharts=geometry(coarse,interface,True);tcharts=geometry(shared,interface,True)
    gammas=[];radius_checks=0
    for i,row in enumerate(g['cells']):
        require(row['source_boundary_tetrahedron']==i,'gamma source simplex index')
        f=shared['facets'][i];cr=cap['boundary_to_source_refinement'][i];co=cocharts[cr['parent']][2];w=mat(cr['barycentric']);src=[tuple(sum(r[j]*co[j][k] for j in range(4)) for k in range(3)) for r in w];dst=tcharts[i][2]
        expected=[gi[(v,int(q[2]==1))] for v,q in zip(f,src)];require(row['vertices']==expected,'gamma vertex binding')
        require(all(gs[v]==x and gt[v]==y for v,x,y in zip(expected,src,dst)),'gamma coordinates differ from cap map')
        A,B=mat(row['affine']),mat(row['inverse']);require(mul(A,B)==eye(4)==mul(B,A),'gamma inverse');require(determinant(A)>0,'gamma orientation')
        require(all(image((*gs[v],Q(1)),A)==(*gt[v],Q(1)) for v in expected),'gamma affine vertex equation')
        require([r[2] for r in A]==[0,0,1,0] and [r[3] for r in A]==[0,0,0,1],'gamma time/affine coordinates')
        # The radial function is linear in the single quadrant occupied by
        # each source and target cell.  Check its full coefficient identity.
        signs=[]
        for points in (src,dst):
            sg=[]
            for k in (0,1):
                nonzero={1 if q[k]>0 else -1 for q in points if q[k]};require(len(nonzero)==1,'gamma crosses radial quadrant');sg.append(next(iter(nonzero)))
            signs.append(sg)
        target_radius=[[Q(signs[1][0])],[Q(signs[1][1])],[Q(0)],[Q(0)]]
        require(mul(A,target_radius)==[[Q(signs[0][0])],[Q(signs[0][1])],[Q(0)],[Q(0)]],'gamma radial preservation');radius_checks+=1
        F,Fi=mat(row['F2_affine']),mat(row['F2_inverse']);require(mul(F,Fi)==eye(4)==mul(Fi,F),'F2 inverse')
        for v in expected:
            x,y,s=gs[v];u,vv,ss=gt[v];require(image((x,y,1-s,Q(1)),F)==(u,vv,1-s,Q(1)),'F2 pullback equation')
        require(determinant(F)>0,'F2 orientation');gammas.append((A,B,F,Fi))
    require(all(gd[gd[v]]==v and gkeys[gd[v]]==(shared['deck'][gkeys[v][0]],gkeys[v][1]) for v in range(len(gd))),'gamma deck binding')
    require(all(gs[gd[v]]==(-gs[v][0],-gs[v][1],gs[v][2]) and gt[gd[v]]==(-gt[v][0],-gt[v][1],gt[v][2]) for v in range(len(gd))),'gamma deck coordinates')
    require(all(gs[v]==gt[v] for v in range(len(gs)) if gs[v][2] in (0,1)),'R0 or R1 is not identity')
    require(all(gs[v]==gt[v] for v in range(len(gs)) if gs[v][0]==gs[v][1]==0),'gamma moves a fixed pole')
    gfacets=Counter(tuple(sorted(f)) for row in g['cells'] for f in combinations(row['vertices'],3));require(all(n in (1,2) for n in gfacets.values()),'gamma facet incidence')
    for s in (Q(0),Q(1)):
        end=[f for f,n in gfacets.items() if n==1 and all(gs[v][2]==s for v in f)];sphere2(end)
    require(all(n==2 or any(all(gs[v][2]==s for v in f) for s in (0,1)) for f,n in gfacets.items()),'unexpected gamma boundary')
    gsimp={tuple(sorted(r['vertices'])):i for i,r in enumerate(g['cells'])};gpartner=[]
    DD=eye(4);DD[0][0]=DD[1][1]=Q(-1)
    for i,row in enumerate(g['cells']):
        j=gsimp[tuple(sorted(gd[v] for v in row['vertices']))];gpartner.append(j);require(mul(DD,gammas[j][0])==mul(gammas[i][0],DD),'gamma deck affine equation')
    ps=list(map(vec,p['source_points']));pt=list(map(vec,p['target_points']));keys=[(int(r[0]),Q(r[1]),Q(r[2])) for r in p['vertex_keys']];pi={k:i for i,k in enumerate(keys)};D=p['deck'];require(len(p['cells'])==2*len(g['cells']),'Phi branch coverage')
    for i,(v,q,t) in enumerate(keys):
        require(0<=q<=1 and 0<=t<=1 and gs[v][2]==min(1-q,t),'Phi vertex profile')
        require(ps[i]==(*gs[v][:2],q,t) and pt[i]==(*gt[v][:2],q,t),'Phi vertex coordinate binding')
        require(D[i]==pi[(gd[v],q,t)],'Phi vertex deck binding')
    facets=defaultdict(list);cellmap={};seams=0;determinants=[]
    for index,row in enumerate(p['cells']):
        i,branch=row['gamma_cell'],row['branch'];require(branch in ('A','B') and (i,branch) not in cellmap,'duplicate Phi branch');cellmap[(i,branch)]=index
        exp=set()
        for v in g['cells'][i]['vertices']:
            s=gs[v][2];qt=((Q(0),s),(1-s,s)) if branch=='A' else ((1-s,s),(1-s,Q(1)))
            exp.update(pi[(v,q,t)] for q,t in qt)
        require(sorted(exp)==row['vertices'],'Phi is not the complete affine pullback polytope')
        require(arank([ps[v] for v in row['vertices']])==4,'Phi cell dimension')
        A,B=mat(row['affine']),mat(row['inverse']);require(mul(A,B)==eye(5)==mul(B,A),'Phi affine inverse');determinants.append(determinant(A));require(determinants[-1]>0,'Phi orientation')
        require(all(image((*ps[v],Q(1)),A)==(*pt[v],Q(1)) for v in row['vertices']),'Phi affine map on polytope vertices')
        require([r[2] for r in A]==[0,0,1,0,0] and [r[3] for r in A]==[0,0,0,1,0],'Phi moves q or t')
        sourcefacets=[];targetfacets=[]
        for hs,points,out in ((mat(row['source_halfspaces']),ps,sourcefacets),(mat(row['target_halfspaces']),pt,targetfacets)):
            require(len(hs)==6,'Phi halfspace count')
            for h in hs:
                values={v:sum(x*y for x,y in zip((*points[v],Q(1)),h)) for v in row['vertices']};require(min(values.values())>=0,'Phi polytope violates halfspace')
                f=tuple(v for v in row['vertices'] if values[v]==0)
                if len(f)>=4 and arank([points[v] for v in f])==3:out.append(f)
        require(sorted(set(sourcefacets))==sorted(set(targetfacets))==list(map(tuple,row['facets'])),'Phi face-poset map')
        # Independently recover every facet from the affine extrusion of
        # the four Gamma faces and its two ends, without trusting halfspaces.
        expected_faces=set()
        for gf in combinations(g['cells'][i]['vertices'],3):
            face=tuple(v for v in row['vertices'] if keys[v][0] in gf)
            if len(face)>=4 and arank([ps[v] for v in face])==3:expected_faces.add(face)
        ends=((lambda q,t:q==0),(lambda q,t:q+t==1)) if branch=='A' else ((lambda q,t:q+t==1),(lambda q,t:t==1))
        for condition in ends:
            face=tuple(v for v in row['vertices'] if condition(ps[v][2],ps[v][3]))
            if len(face)>=4 and arank([ps[v] for v in face])==3:expected_faces.add(face)
        require(expected_faces==set(sourcefacets),'Phi extrusion facet completeness')
        for f in set(sourcefacets):facets[f].append(index)
        # Bind the affine map to R_s on the full cell.  Agreement on its
        # affinely spanning vertices implies equality everywhere.
        G=gammas[i][0]
        for v in row['vertices']:
            x,y,q,t=ps[v];s=t if branch=='A' else 1-q;z=image((x,y,s,Q(1)),G);require(pt[v]==(z[0],z[1],q,t),'Phi actual Gamma pullback')
    require(len(cellmap)==2*len(g['cells']),'Phi branch count')
    DD5=eye(5);DD5[0][0]=DD5[1][1]=Q(-1)
    for i,row in enumerate(g['cells']):
        aa=p['cells'][cellmap[(i,'A')]];bb=p['cells'][cellmap[(i,'B')]];A,B=mat(aa['affine']),mat(bb['affine']);gamma=gammas[i][0]
        # A-B is exactly (q+t-1) times the s-column of R_s.
        expected=[[Q(0) for _ in range(5)] for _ in range(5)]
        for j in (0,1):expected[2][j]=expected[3][j]=gamma[2][j];expected[4][j]=-gamma[2][j]
        require([[x-y for x,y in zip(r,s)] for r,s in zip(A,B)]==expected,'Phi branch seam equation');seams+=1
        # At t=1 the B matrix is exactly the F2 matrix, with q retained.
        F=gammas[i][2]
        restricted=[[B[k][j] if k!=3 else B[3][j]+B[4][j] for j in (0,1,2,4)] for k in (0,1,2,3)]
        require(restricted==F,'Phi final-time seam to F2')
        for branch in ('A','B'):
            C=mat(p['cells'][cellmap[(i,branch)]]['affine']);other=mat(p['cells'][cellmap[(gpartner[i],branch)]]['affine']);require(mul(DD5,other)==mul(C,DD5),'Phi deck affine equation')
    require([{'vertices':list(f),'incidences':v} for f,v in sorted(facets.items())]==p['facets'],'Phi global facet table')
    require(all(len(v) in (1,2) for v in facets.values()),'Phi facet multiplicity')
    boundary=[]
    for f,inc in facets.items():
        if len(inc)==1:
            sides=[label for coord,value,label in ((2,0,'q0'),(2,1,'q1'),(3,0,'t0'),(3,1,'t1')) if all(ps[v][coord]==value for v in f)];require(len(sides)==1,'unexpected Phi boundary')
            boundary.append({'vertices':list(f),'cell':inc[0],'side':sides[0]})
    require(sorted(boundary,key=lambda r:r['vertices'])==p['boundary'],'Phi boundary table')
    # Endpoint checks on full boundary polyhedra follow from their vertices.
    require(all(ps[v]==pt[v] for v in range(len(ps)) if ps[v][2]==1 or ps[v][3]==0),'Phi bulk/initial boundary not identity')
    # At q=0 the construction uses exactly the actual R_t cap map.
    for v,(gv,q,t) in enumerate(keys):
        if q==0:require(gs[gv][2]==t and ps[v][:2]==gs[gv][:2] and pt[v][:2]==gt[gv][:2],'Phi puncture restriction')
    lifts=data['lifts'];require(lifts['constant_lattice_adjustment']==1 and lifts['bulk_D2_unwrapped_shift']==1 and lifts['cap_displacement']==[-8,4] and lifts['cap_q0_upstairs_winding']==-1 and lifts['cap_q0_downstairs_winding']==-2,'integer lift/winding data')
    # Coefficients in (theta,q,1): D2 + Phi_(t=1) = F2 exactly.
    require([Q(1),Q(1),Q(-1)+lifts['bulk_D2_unwrapped_shift']]==[Q(1),Q(1),Q(0)],'unwrapped seam equation')
    return {'result':'PASS','scope':data['scope'],'source_global_identification':data['source_global_identification'],'gamma_tetrahedra':len(g['cells']),'gamma_radius_identities':radius_checks,'R0_R1_identity':True,'F2_affine_maps_and_inverses':len(g['cells']),'Phi_affine_4_cells':len(p['cells']),'Phi_affine_inverses':len(p['cells']),'Phi_facet_incidence_checks':len(facets),'Phi_branch_seam_identities':seams,'Phi_final_time_F2_identities':seams,'Phi_boundary_facets_by_side':dict(Counter(r['side'] for r in boundary)),'deck_equivariant':True,'initial_time_and_bulk_identity':True,'puncture_restriction_is_actual_beta2':True,'unwrapped_D2_Phi1_equals_F2':True,'source_angle_cycle':interface['equator_cycle_diagram_positive'],'cap_displacement':[-8,4],'simplicial_expansion':'Not materialized; canonical barycentric face chains give an explicit equivariant simplicial realization'}


if __name__=='__main__':
    root=Path(__file__).resolve().parent;p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=root);p.add_argument('--output',type=Path,default=root/'transfer_verification.json');a=p.parse_args();root=a.root;start=time.monotonic();data=load(root/'data/finite_F2_Phi_transfer.json.gz');shared=load(root/'data/shared_boundary.json.gz');cap=load(root/'data/cap_m2.json.gz');bridge=load(root/'data/source_cap_bridge.json.gz')
    for key,fn in (('shared_boundary_file_sha256','shared_boundary.json.gz'),('cap_m2_file_sha256','cap_m2.json.gz'),('source_bridge_file_sha256','source_cap_bridge.json.gz')):require(data['bindings'][key]==hashlib.sha256((root/'data'/fn).read_bytes()).hexdigest(),'transfer dependency bytes: '+fn)
    result=check(data,shared,cap,bridge);result['negative_controls']=[]
    for name,mutate in [('drop_bulk_full_turn',lambda x:x['lifts'].__setitem__('bulk_D2_unwrapped_shift',0)),('corrupt_Phi_affine_cell',lambda x:x['Phi']['cells'][0]['affine'][0].__setitem__(0,'99')),('omit_Phi_cell',lambda x:x['Phi']['cells'].pop())]:
        x=copy.deepcopy(data);mutate(x)
        try:check(x,shared,cap,bridge)
        except (ValueError,KeyError,IndexError) as e:result['negative_controls'].append({'name':name,'rejected':True,'reason':str(e)})
        else:raise ValueError('accepted corrupt transfer: '+name)
    result['seconds']=time.monotonic()-start;a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

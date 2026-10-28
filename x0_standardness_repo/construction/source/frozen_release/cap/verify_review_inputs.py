#!/usr/bin/env python3
"""Independent saved-data check of the repaired localization review inputs.

Reads the mixed-profile matrices as stored; does not regenerate or import the
mixed-profile producer/checker. No native or large homology calls are used.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from itertools import combinations
import copy,gzip,hashlib,json,time

ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent

def require(ok,message):
    if not ok:raise ValueError(message)
def load(path):
    b=path.read_bytes()
    return json.loads(gzip.decompress(b) if path.suffix=='.gz' else b)
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def tup(x):return tuple(map(tup,x)) if isinstance(x,(tuple,list)) else x
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def det(a):
    a=[list(r) for r in a];out=Q(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:return Q(0)
        if i!=j:a[i],a[j]=a[j],a[i];out=-out
        pivot=a[i][i];out*=pivot
        for j in range(i+1,len(a)):
            factor=a[j][i]/pivot
            for k in range(i,len(a)):a[j][k]-=factor*a[i][k]
    return out

def check_profile(rec,central,model,cut):
    require(rec['central_certificate_sha256']==digest(central),'mixed central hash')
    require(rec['adaptor_model_sha256']==digest(model),'mixed adaptor hash')
    require(rec['source_cut_sha256']==digest(cut),'mixed source hash')
    require((rec['roll'],rec['twist'])==(1,2),'mixed coefficients')
    require(rec['source_winding_upstairs']==[1,1] and rec['source_winding_downstairs']==[1,2],'mixed integer winding')
    require(len(rec['cells'])==6,'mixed profile cell count')
    edges=defaultdict(list);area=Q(0)
    for i,row in enumerate(rec['cells']):
        k=i//2
        want=((Q(k),Q(0)),(Q(k+1),Q(0)),(Q(k+1),Q(1))) if i%2==0 else ((Q(k),Q(0)),(Q(k+1),Q(1)),(Q(k),Q(1)))
        points=tuple(tuple(map(Q,p)) for p in row['vertices']);require(points==want,'mixed profile partition')
        a,b,c=map(Q,row['J_affine_coefficients'])
        for r,t in points:
            expected=Q(0) if t==0 else r/3-min(Q(1),max(Q(0),r-1))
            require(a*r+b*t+c==expected,'mixed profile coefficient')
        F=[list(map(Q,r)) for r in row['forward_homogeneous_matrix']]
        G=[list(map(Q,r)) for r in row['inverse_homogeneous_matrix']]
        expected=eye(5)
        for j in (0,1):expected[j][2:5]=[a,b,c]
        require(F==expected,'mixed saved map does not shift both u and v')
        require(mm(F,G)==mm(G,F)==eye(5),'mixed saved affine inverse')
        require(det(F)==1,'mixed affine orientation')
        for e in combinations(points,2):edges[tuple(sorted(e))].append((a,b,c))
        p,q,r=points;area+=((q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0]))/2
    require(area==3 and all(len(v) in (1,2) for v in edges.values()),'mixed profile exact coverage')
    boundary_count=0
    for edge,coefficients in edges.items():
        if len(coefficients)==2:
            require(all(len({a*r+b*t+c for a,b,c in coefficients})==1 for r,t in edge),'mixed shared edge')
        else:
            boundary_count+=1;a,b,c=coefficients[0]
            require(all(r==0 for r,t in edge) or all(r==3 for r,t in edge) or all(t==0 for r,t in edge) or all(t==1 for r,t in edge),'unexpected profile boundary')
            for r,t in edge:
                if r in (0,3) or t==0:require(a*r+b*t+c==0,'mixed profile moves puncture boundary')
                if t==1:require(a*r+b*t+c+min(Q(1),max(Q(0),r-1))==r/3,'mixed full seam identity')
    return {'stored_affine_maps':6,'stored_affine_inverses':6,'exact_domain_area':'3','boundary_segments':boundary_count,'both_angular_components_checked':True,'puncture_fixed_and_integer_seam_checked':True}

def check_enlargement(model,cut,chart):
    equator=[56,150,58,148];poles=[204,228]
    require(cut['puncture']['fixed_edge_cut_vertices']==poles,'puncture poles')
    raw=cut['puncture']['equator_cycle'];require([raw[0],*reversed(raw[1:])]==equator,'puncture positive cycle')
    fixed=set(cut['fixed_source_vertices']);levels=[]
    for v in equator:
        p=cut['points'][v];require(p[0]=='m','equator cut midpoint')
        levels.append(Q(sum(x in fixed for x in p[1:]),2))
    require(levels==[Q(1,2)]*4,'original puncture torus contact')
    uv={v:(Q(u),Q(w)) for v,u,w in chart['upstairs_vertex_coordinates']}
    require(len({uv[v][0] for v in equator})==2,'original equator nonconstant longitude')
    full=list(map(tup,model['stages']['N']))+[tuple(('N','o',v) for v in t) for t in cut['puncture']['removed_N_tetrahedra']]
    require(len({frozenset(t) for t in full})==len(full),'filled enlarged tube duplicate')
    faces=Counter(frozenset(f) for t in full for f in combinations(t,3))
    require(all(n in (1,2) for n in faces.values()),'filled enlarged tube facet incidence')
    boundary={f for f,n in faces.items() if n==1}
    require(len(boundary)==24 and all(v[0:2]==('N','g') for f in boundary for v in f),'enlarged tube outer grid')
    sphere={frozenset(('N','o',v) for v in f) for f in cut['puncture']['boundary_triangles']}
    require(len(sphere)==8 and all(faces[f]==2 for f in sphere),'puncture not interior after filling')
    require(not {v for f in sphere for v in f}&{v for f in boundary for v in f},'puncture meets enlarged tube boundary')
    return {'filled_N_plus_tetrahedra':len(full),'outer_grid_triangles':len(boundary),'interior_puncture_triangles':len(sphere),'puncture_disjoint_from_outer_grid':True,'original_equator_levels':['1/2']*4,'original_equator_u_values':sorted({str(uv[v][0]) for v in equator}),'positive_equator':equator,'poles':poles}

if __name__=='__main__':
    start=time.monotonic()
    names={'mixed_profile':'central_roll/mixed_profile_transfer.json','central':'central_roll/certificate.json.gz','adaptor':'adaptor/model.json.gz','cut':'peripheral/cut_source.json.gz','chart':'peripheral/diagram_chart.json.gz','proof':'central_roll/LOCALIZATION_LEMMA.txt'}
    data={k:load(PARENT/v) for k,v in names.items() if k!='proof'}
    profile=check_profile(data['mixed_profile'],data['central'],data['adaptor'],data['cut'])
    enlarged=check_enlargement(data['adaptor'],data['cut'],data['chart'])
    wrong=copy.deepcopy(data['mixed_profile']);wrong['cells'][0]['forward_homogeneous_matrix'][1][3]='0'
    try:check_profile(wrong,data['central'],data['adaptor'],data['cut'])
    except ValueError as e:negative={'name':'drop_v_profile_correction','rejected':True,'reason':str(e)}
    else:raise ValueError('mixed-profile corruption escaped')
    out={'schema':'independent-localization-branch-review-v1','finite_inputs_result':'PASS','bindings':{k:{'path':v,'file_sha256':hashlib.sha256((PARENT/v).read_bytes()).hexdigest()} for k,v in names.items()},'profile':profile,'enlarged_neighborhood':enlarged,'negative_control':negative,'geometric_review':'Revised localization argument accepted with the stated relative PL normal-neighborhood and marked boundary-extension theorem dependency; no global finite normal atlas is claimed.','branch_model_review':'Piecewise equivariant product/collar/cone geometry plus relative deck-paired ball replacements proves all local branch models, conditional on parent actual assembly port/incidence/deck checks.','not_proved':'PL standardness of closed X2; entirely finite global pointwise conjugacy.','seconds':time.monotonic()-start}
    (ROOT/'geometric_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

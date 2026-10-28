#!/usr/bin/env python3
"""Finite rational polyhedral F2/Phi transfer from the actual spherical loop.

The output keeps the local collar as affine polytopes rather than expanding a
large barycentric 4D triangulation.  Every cell has its map, inverse, vertices,
halfspaces and facets.  Canonical barycentric subdivision makes it simplicial.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict,Counter
import argparse,hashlib,json,time
from build_caps import save,digest
from build_bridge import load,cap_charts,detmatrix


def inverse(a):
    n=len(a);a=[list(r)+[Q(i==j) for j in range(n)] for i,r in enumerate(a)]
    for j in range(n):
        p=next((i for i in range(j,n) if a[i][j]),None)
        if p is None:raise ValueError('singular affine matrix')
        a[j],a[p]=a[p],a[j];s=a[j][j];a[j]=[x/s for x in a[j]]
        for i in range(n):
            if i==j:continue
            s=a[i][j];a[i]=[x-s*y for x,y in zip(a[i],a[j])]
    return [r[n:] for r in a]


def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def apply(v,a):return [sum(x*y for x,y in zip(v,col)) for col in zip(*a)]


def rank(vs):
    if not vs:return -1
    a=[[x-y for x,y in zip(v,vs[0])] for v in vs[1:]];r=0
    for j in range(len(vs[0])):
        p=next((i for i in range(r,len(a)) if a[i][j]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];s=a[r][j];a[r]=[x/s for x in a[r]]
        for i in range(r+1,len(a)):
            s=a[i][j];a[i]=[x-s*y for x,y in zip(a[i],a[r])]
        r+=1
    return r


def phi_matrix(a,branch):
    out=[[Q(0) for _ in range(5)] for _ in range(5)]
    for j in (0,1):
        out[0][j]=a[0][j];out[1][j]=a[1][j]
        if branch=='A':out[3][j]=a[2][j];out[4][j]=a[3][j]
        else:out[2][j]=-a[2][j];out[4][j]=a[3][j]+a[2][j]
    for j in (2,3,4):out[j][j]=Q(1)
    return out


def f_matrix(a):
    out=[[Q(0) for _ in range(4)] for _ in range(4)]
    for j in (0,1):out[0][j]=a[0][j];out[1][j]=a[1][j];out[2][j]=-a[2][j];out[3][j]=a[2][j]+a[3][j]
    out[2][2]=out[3][3]=Q(1);return out


def pull_inequalities(tet,branch):
    k=inverse([list(v)+[Q(1)] for v in tet]);out=[]
    for i in range(4):
        if branch=='A':out.append([k[0][i],k[1][i],Q(0),k[2][i],k[3][i]])
        else:out.append([k[0][i],k[1][i],-k[2][i],Q(0),k[2][i]+k[3][i]])
    if branch=='A':out.extend([[Q(0),Q(0),Q(1),Q(0),Q(0)],[Q(0),Q(0),Q(-1),Q(-1),Q(1)]])
    else:out.extend([[Q(0),Q(0),Q(1),Q(1),Q(-1)],[Q(0),Q(0),Q(0),Q(-1),Q(1)]])
    return out


def build(root):
    root=Path(root);shared_path=root/'data/shared_boundary.json.gz';cap_path=root/'data/cap_m2.json.gz';bridge_path=root/'data/source_cap_bridge.json.gz'
    b=load(shared_path);c=load(cap_path);bridge=load(bridge_path);coarse=dict(b['coarse_source_boundary']);coarse['times']=b['times'];cocharts=cap_charts(coarse);target=cap_charts(b)
    gpoints={};rows=[]
    for i,(f,row,tchart) in enumerate(zip(b['facets'],c['boundary_to_source_refinement'],target)):
        co=cocharts[row['parent']]['points'];mat=[list(map(Q,r)) for r in row['barycentric']]
        src=[tuple(sum(r[j]*co[j][k] for j in range(4)) for k in range(3)) for r in mat];dst=tchart['points'];keys=[]
        for v,p,z in zip(f,src,dst):
            if p[2]!=z[2]:raise ValueError('gamma does not preserve time')
            key=(v,int(p[2]==1));keys.append(key)
            if key in gpoints and gpoints[key]!=(p,z):raise ValueError('gamma point mismatch')
            gpoints[key]=(p,z)
        hs=[list(p)+[Q(1)] for p in src];ht=[list(p)+[Q(1)] for p in dst];a=mm(inverse(hs),ht);ai=mm(inverse(ht),hs)
        if detmatrix(a)<=0:raise ValueError('gamma orientation')
        rows.append({'source_boundary_tetrahedron':i,'vertices':keys,'affine':a,'inverse':ai,'F2_affine':f_matrix(a),'F2_inverse':f_matrix(ai)})
    gkeys=sorted(gpoints);gi={v:i for i,v in enumerate(gkeys)};gs=[gpoints[k][0] for k in gkeys];gt=[gpoints[k][1] for k in gkeys];gd=[gi[(b['deck'][v],lift)] for v,lift in gkeys]
    for row in rows:row['vertices']=[gi[v] for v in row['vertices']]
    vertices={};cells=[]
    for i,row in enumerate(rows):
        for branch in ('A','B'):
            keys=[]
            for v in row['vertices']:
                s=gs[v][2]
                qt=((Q(0),s),(1-s,s)) if branch=='A' else ((1-s,s),(1-s,Q(1)))
                for q,t in qt:
                    key=(v,q,t);vertices[key]=((*gs[v][:2],q,t),(*gt[v][:2],q,t));keys.append(key)
            keys=sorted(set(keys));cells.append({'gamma_cell':i,'branch':branch,'vertices':keys,'affine':phi_matrix(row['affine'],branch),'inverse':phi_matrix(row['inverse'],branch),'source_halfspaces':pull_inequalities([gs[v] for v in row['vertices']],branch),'target_halfspaces':pull_inequalities([gt[v] for v in row['vertices']],branch)})
    keys=sorted(vertices);vi={v:i for i,v in enumerate(keys)};ps=[vertices[k][0] for k in keys];pt=[vertices[k][1] for k in keys];D=[vi[(gd[v],q,t)] for v,q,t in keys];facets=defaultdict(list)
    for i,row in enumerate(cells):
        row['vertices']=[vi[v] for v in row['vertices']];ff=set()
        for h in row['source_halfspaces']:
            f=tuple(v for v in row['vertices'] if sum(x*y for x,y in zip((*ps[v],Q(1)),h))==0)
            if len(f)>=4 and rank([ps[v] for v in f])==3:ff.add(f)
        row['facets']=sorted(ff)
        for f in ff:facets[f].append(i)
    if any(len(x) not in (1,2) for x in facets.values()):raise ValueError('pullback facet incidence')
    boundary=[]
    for f,inc in facets.items():
        if len(inc)==1:
            sides=[label for coord,val,label in ((2,Q(0),'q0'),(2,Q(1),'q1'),(3,Q(0),'t0'),(3,Q(1),'t1')) if all(ps[v][coord]==val for v in f)]
            if len(sides)!=1:raise ValueError(('unexpected transfer boundary',f,sides))
            boundary.append({'vertices':f,'cell':inc[0],'side':sides[0]})
    result={'schema':'finite-F2-Phi-polyhedral-transfer-v1','scope':'Source-bound local collar PL maps; canonical barycentric subdivision is specified but not expanded into 4-simplices','source_global_identification':'The actual new puncture collar and its maps are meshed. Identifying this relative F2 with the original torus-supported twist uses the oriented normal-disk-bundle localization theorem; no global fiber-action atlas for N is exported here.','bindings':{'shared_boundary_file_sha256':hashlib.sha256(shared_path.read_bytes()).hexdigest(),'cap_m2_file_sha256':hashlib.sha256(cap_path.read_bytes()).hexdigest(),'source_bridge_file_sha256':hashlib.sha256(bridge_path.read_bytes()).hexdigest(),'source_interface':bridge['source_interface']},'coordinates':['x','y','q','t'],'homogeneous_matrix_convention':'row vector [coordinates,1] times matrix','gamma':{'definition':'R_s is the actual beta2 time-family from the spherical cap refinement','vertex_keys':gkeys,'source_points':gs,'target_points':gt,'deck':gd,'cells':rows},'Phi':{'definition':'Phi(x,q,t)=(R_min(1-q,t)(x),q,t)','vertex_keys':keys,'source_points':ps,'target_points':pt,'deck':D,'cells':cells,'facets':[{'vertices':f,'incidences':inc} for f,inc in sorted(facets.items())],'boundary':sorted(boundary,key=lambda r:r['vertices'])},'F2':{'definition':'F2(x,q)=(R_(1-q)(x),q)','cells':'gamma.cells supply F2_affine and F2_inverse; use coordinates (x,y,q=1-s)','source_collaring':'Adjoin the explicit product collar of the recorded source puncture sphere. Extend F2 by identity outside; this is relative to both collar ends.'},'lifts':{'accepted_F2_equatorial_lift':'theta + q','gamma_expression_equatorial_lift':'theta + q - 1','constant_lattice_adjustment':1,'Phi_A_equatorial_lift':'theta - t','Phi_B_equatorial_lift':'theta + q - 1','bulk_D2_unwrapped_shift':1,'cap_q0_upstairs_winding':-1,'cap_q0_downstairs_winding':-2,'cap_displacement':[-8,4]},'PL_realization':'Each source/target cell is a rational convex polytope in its hemisphere chart. Affine bijections preserve the face poset. Barycentres of all faces and chains give an explicit canonical equivariant simplicial refinement, without materializing the larger unrelated bulk.'}
    save(root/'data/finite_F2_Phi_transfer.json.gz',result)
    summary={'schema':'finite-transfer-generation-v1','gamma_vertices':len(gkeys),'gamma_tetrahedra':len(rows),'Phi_vertices':len(keys),'Phi_4_cells':len(cells),'Phi_3_facets':len(facets),'boundary_facets_by_side':dict(Counter(r['side'] for r in boundary)),'source_angle_cycle':bridge['source_interface']['equator_cycle_diagram_positive'],'canonical_sha256':digest(result)};save(root/'data/transfer_generation.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args();build(a.root)

#!/usr/bin/env python3
"""Finite longitudinal/time collar for the marked even, positive-roll bulk.

This builds only a three-dimensional T^2 x I.  The fourth, meridional,
coordinate is passive and is a four-edge circle in the bulk assembly.
All rational charts and the source-relative PL profile transfer are saved.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import combinations
import argparse,gzip,hashlib,json,math,sys,time

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from exact_flat_mesh import grid_faces,overlay,cylinder_half

def load(p):
    p=Path(p);b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)

def tup(x):return tuple(tup(z) for z in x) if isinstance(x,(list,tuple)) else x

def conv(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,set):return sorted((conv(v) for v in x),key=repr)
    if isinstance(x,(list,tuple)):return [conv(v) for v in x]
    if isinstance(x,dict):return {str(k):conv(v) for k,v in x.items() if k!='key'}
    return x

def digest(x):return hashlib.sha256(json.dumps(conv(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    b=(json.dumps(conv(x),sort_keys=True,separators=(',',':'))+'\n').encode()
    p.write_bytes(gzip.compress(b,mtime=0) if p.suffix=='.gz' else b)

def inverse(A):
    n=len(A);a=[[Q(v) for v in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        p=next(i for i in range(k,n) if a[i][k]);a[k],a[p]=a[p],a[k]
        d=a[k][k];a[k]=[v/d for v in a[k]]
        for i in range(n):
            if i!=k:
                d=a[i][k];a[i]=[x-d*y for x,y in zip(a[i],a[k])]
    return [r[n:] for r in a]

def determinant3(co):
    a,b,c=[tuple(p[k]-co[0][k] for k in range(3)) for p in co[1:]]
    return a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+c[0]*(a[1]*b[2]-a[2]*b[1])

def bary(P,co):
    M=[[co[j][i] for j in range(4)] for i in range(3)]+[[Q(1)]*4]
    I=inverse(M);w=[sum(a*b for a,b in zip(row,(*P,Q(1)))) for row in I]
    return w if min(w)>=0 else None

def check_branching(tets):
    ed={};seen=set()
    for t in tets:
        if len(set(t))!=4 or frozenset(t) in seen:raise ValueError('folded/duplicate central tetrahedron')
        seen.add(frozenset(t))
        for a,b in combinations(t,2):
            e=frozenset((a,b))
            if e in ed and ed[e]!=(a,b):raise ValueError('central branching conflict')
            ed[e]=(a,b)

def centtag(v,j):
    if v[0]=='a':return ('B',('N','g',v[1],j),v[2])
    if v[0]=='b':return ('B',('E','g',v[1],j),v[2])
    return ('C',v,j)

def profile_transfer():
    cells=[]
    for i in range(3):
        a,b=Q(i),Q(i+1)
        for P in (((a,Q(0)),(b,Q(0)),(b,Q(1))),((a,Q(0)),(b,Q(1)),(a,Q(1)))):
            values=[Q(0) if t==0 else r/3-min(Q(1),max(Q(0),r-1)) for r,t in P]
            M=[[r,t,Q(1)] for r,t in P];I=inverse(M)
            abc=[sum(x*y for x,y in zip(row,values)) for row in I]
            cells.append({'vertices':P,'values':values,'affine_coefficients':abc})
    return {'schema':'relative-PL-source-roll-profile-transfer-v1','roll':1,
            'source_radial_coordinate':'r=3*s_input','source_profile':'H0(r)=r/3',
            'assembled_profile':'H1(r)=clamp(r-1,0,1)',
            'map':'(u,v,r,t)->(u+J(r,t),v,r,t)',
            'inverse':'(u,v,r,t)->(u-J(r,t),v,r,t)',
            'J_triangles':cells,'boundary_identity':['r=0','r=3','t=0'],
            'seam_equation':'J(r,1)+H1(r)=H0(r)',
            'cap_and_puncture_boundary_pointwise_fixed':True,
            'PL_basis':'affine interpolation on six declared triangles, never bilinear t*(H0-H1)',
            'extension':'identity outside the attached source torus collar'}

def straightening():
    return {'schema':'two-piece-PL-roll-straightening-v1',
            'coordinate_order':['u','v','h','t'],
            'map':'(u,v,h,t)->(u+min(h,t),v,h,t)',
            'inverse':'(U,v,h,t)->(U-min(h,t),v,h,t)',
            'cells':[
                {'h_t_vertices':[['0','0'],['1','0'],['1','1']], 'condition':'t<=h',
                 'forward_matrix':[[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
                 'inverse_matrix':[[1,0,0,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]]},
                {'h_t_vertices':[['0','0'],['1','1'],['0','1']], 'condition':'h<=t',
                 'forward_matrix':[[1,0,1,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
                 'inverse_matrix':[[1,0,-1,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]}],
            'inner_attachment':'u_N=U at h=0',
            'outer_attachment':'u_E=U-t at h=1',
            'time_seam':'(u,v,h,1)~(u+h,v,h,0)',
            'assembled_straightened_time_seam':'(U,v,h,1)~(U,v,h,0)',
            'meridian_deck':'(U,v,h,t)->(U,v+1/2,h,t)',
            'full_roll_winding':1}

def locate_central(c,U,t,h,v=Q(1,7)):
    for cell in c['universal_charts']:
        co=[tuple(map(Q,p)) for p in cell['coords']]
        if not min(p[1] for p in co)<=t<=max(p[1] for p in co):continue
        if not min(p[2] for p in co)<=h<=max(p[2] for p in co):continue
        for sh in range(math.ceil(min(p[0] for p in co)-U),math.floor(max(p[0] for p in co)-U)+1):
            w=bary((U+sh,t,h),co)
            if w is None:continue
            base=tup(conv(cell['vertices']));j=int(4*v)%4;eta=4*v-j
            for k in range(4):
                top=eta-sum(w[k+1:])
                if 0<=top<=w[k]:
                    path=[*((base[i],j) for i in range(k+1)),*((base[i],(j+1)%4) for i in range(k,4))]
                    ww=w[:k]+[w[k]-top,top]+w[k+1:]
                    return {'part':'central','base_tetrahedron':base,'base_barycentric_weights':w,
                            'coordinates':[U,t,h],'meridian_coordinate':v,'longitude_lattice_translation':sh,
                            'expected_product_carrier':[centtag(a,b) for a,b in path],
                            'product_barycentric_weights':ww}
    raise ValueError('central source witness has no exact carrier')

def source_witnesses(c,model,coarsening):
    endpoint={s:{frozenset(tup(t)):tup(t) for t in ts} for s,ts in model['stages'].items()}
    protected=coarsening['protected_source_witnesses'];out=[];u,v=Q(1,17),Q(1,7)
    for s in (Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)):
        r=3*s;H=min(Q(1),max(Q(0),r-1));entry={'input_s':s,'input_u':u,'input_v':v,
            'source_coordinate_seam_displacement':[s,Q(0)],'full_unwrapped_source_winding':s,'carriers':[]}
        for tm,up in ((1,u+s-H),(0,u+s)):
            if 1<r<2:
                cc=locate_central(c,up+min(H,Q(tm)),Q(tm),H,v)
                cc.update({'source_u_after_profile_transfer':up,'time':tm})
            else:
                cc0=next(w for w in protected if w['roll']==1 and Q(w['input_s'])==s and w['time']==tm)
                side=cc0['side'];base=endpoint[side][frozenset(tup(cc0['vertices']))]
                old=dict(zip(map(tup,cc0['vertices']),map(Q,cc0['barycentric_weights'])));w=[old[a] for a in base]
                if tm==0:pv=[('B',a,0) for a in base]+[('B',base[-1],1)];ww=w+[Q(0)]
                else:pv=[('B',base[0],3)]+[('B',a,0) for a in base];ww=[Q(0)]+w
                cc={'part':side,'base_tetrahedron':base,'base_barycentric_weights':w,
                    'source_u_after_profile_transfer':up,'time':tm,'source_adaptor_witness':cc0,
                    'expected_product_carrier':pv,'product_barycentric_weights':ww}
            entry['carriers'].append(cc)
        entry['seam_lattice_residue']=Q(-1) if r>=2 else Q(0)
        out.append(entry)
    return out

def build(output=ROOT/'central_roll'):
    start=time.monotonic();out=Path(output)
    if (out/'model.json.gz').exists():raise FileExistsError('use a fresh central output directory')
    cut=load(ROOT/'peripheral/cut_source.json.gz');chart=load(ROOT/'peripheral/diagram_chart.json.gz')
    adaptor=load(ROOT/'adaptor/model.json.gz');coarse=load(ROOT/'adaptor/coarsening_certificate.json.gz')
    if adaptor['diagram_chart_sha256']!=digest(chart) or adaptor['source_cut_sha256']!=digest(cut):raise ValueError('source binding')
    A=grid_faces(3,4);B=grid_faces(3,4,1);O=overlay(A,B)
    aa,ca=cylinder_half(A,O,'a',lambda z:('a',*z),lambda p:('r',*p),lambda i:('ia',i))
    bb,cb=cylinder_half(B,O,'b',lambda z:('b',*z),lambda p:('r',*p),lambda i:('ib',i))
    tets=aa+bb;check_branching(tets);uc=[]
    for half,cc in (('a',ca),('b',cb)):
        for cell in cc:
            co=[(p[0],p[1],p[2]/2 if half=='a' else 1-p[2]/2) for p in cell['coords']]
            M=[[co[j][i] for j in range(4)] for i in range(3)]+[[Q(1)]*4]
            uc.append({'half':half,'parent':cell['parent'],'vertices':cell['vertices'],'coords':co,
                'barycentric_inverse':inverse(M),'determinant':determinant3(co)})
    C={'schema':'source-bound-even-positive-roll-central-collar-v1','L':3,'T':4,'V':4,'roll':1,
       'A':A,'B':B,'overlay':O,'charts_a':ca,'charts_b':cb,'tetrahedra':tets,'universal_charts':uc,
       'physical_radial_coordinate':{'a':'h=local_h/2','b':'h=1-local_h/2','source':'r=1+h'},
       'coordinate_order':['U','time','h'],'meridional_coordinate':'v is a separate four-edge circle',
       'product_tags':{'a':"('B',('N','g',i,jv),time)",'b':"('B',('E','g',i,jv),time)",'interior':"('C',central_vertex,jv)"},
       'deck_rule':'central_vertex unchanged; jv -> (jv+2) mod 4',
       'time_seam':'identity on N and E products and on straightened central coordinates; no meridian shift',
       'source_bindings':{'source_cut_sha256':digest(cut),'diagram_chart_sha256':digest(chart),'adaptor_model_sha256':digest(adaptor),
                          'adaptor_coarsening_certificate_sha256':digest(coarse)},
       'coordinate_meaning':chart['coordinate_meaning'],'diagram_sign_transport':chart['diagram_sign_transport'],
       'straightening':straightening(),'source_profile_transfer':profile_transfer(),
       'scope':'ordered 3D central collar, exact PL charts and expected product carriers; actual 4D assembly membership is checked separately'}
    C['source_seam_witnesses']=source_witnesses(C,adaptor,coarse)
    incidence=Counter(frozenset(f) for t in tets for f in combinations(t,3))
    want={frozenset(('a',*v) for v in f['ids']) for f in A}|{frozenset(('b',*v) for v in f['ids']) for f in B}
    if {f for f,m in incidence.items() if m==1}!=want or any(m not in (1,2) for m in incidence.values()):raise ValueError('central collar boundary')
    model={'schema':'ordered-central-roll-collar-model-v1','L':3,'T':4,'V':4,'roll':1,'tetrahedra':tets,
           'central_certificate_sha256':digest(C),'source_bindings':C['source_bindings'],
           'boundary_a':[tuple(('a',*v) for v in f['order']) for f in A],
           'boundary_b':[tuple(('b',*v) for v in f['order']) for f in B],
           'deck_on_meridian_indices':[2,3,0,1]}
    save(out/'certificate.json.gz',C);save(out/'model.json.gz',model)
    summary={'schema':'central-roll-build-summary-v1','roll':1,'tetrahedra':len(tets),
             'overlay_polygons':len(O['pieces']),'overlay_triangles':len(O['triangles']),
             'boundary_triangles':len(want),'meridian_product_pentachora':16*len(tets),
             'model_sha256':digest(model),'certificate_sha256':digest(C),'seconds':time.monotonic()-start}
    save(out/'summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    return C,model

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'central_roll')
    args=p.parse_args();build(args.output)

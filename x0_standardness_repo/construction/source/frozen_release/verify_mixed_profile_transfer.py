#!/usr/bin/env python3
"""Exact simultaneous u/v profile transfer, bound to the central certificate.

This data-only checker does not import either geometric producer.  The
original six-triangle J is read literally and its affine equations checked
independently.  No new subdivision or bulk construction is performed.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import combinations
import gzip,hashlib,json,time

ROOT=Path(__file__).resolve().parent
def req(ok,msg):
    if not ok:raise ValueError(msg)
def load(p):
    p=ROOT/p;b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def dig(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def mul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]

def verify():
    start=time.monotonic();C=load('central_roll/certificate.json.gz');M=load('adaptor/model.json.gz');A=load('adaptor/adaptor.json.gz');cut=load('peripheral/cut_source.json.gz')
    req(C['source_bindings']['adaptor_model_sha256']==dig(M) and M['adaptor_sha256']==dig(A),'mixed profile source binding')
    req(C['source_bindings']['source_cut_sha256']==dig(cut),'mixed profile actual puncture source')
    cells=C['source_profile_transfer']['J_triangles'];req(len(cells)==6,'mixed profile cells')
    records=[]
    for i,cell in enumerate(cells):
        k=i//2;a,b=Q(k),Q(k+1)
        want=((a,0),(b,0),(b,1)) if i%2==0 else ((a,0),(b,1),(a,1))
        ps=tuple(tuple(map(Q,p)) for p in cell['vertices']);req(ps==want,'mixed profile rectangle subdivision')
        aa,bb,cc=map(Q,cell['affine_coefficients'])
        for (r,t),val in zip(ps,map(Q,cell['values'])):
            H1=min(Q(1),max(Q(0),r-1));expected=Q(0) if t==0 else r/3-H1
            req(aa*r+bb*t+cc==val==expected,'mixed profile affine equations')
            if r in (0,3):req(val==0,'mixed profile moves original collar boundary')
        # Coordinates are homogeneous (u,v,r,t,1).  Both angular rows get
        # the same affine scalar, and r,t are unchanged.
        F=eye(5);G=eye(5)
        for row in (0,1):F[row][2:5]=[aa,bb,cc];G[row][2:5]=[-aa,-bb,-cc]
        req(mul(F,G)==mul(G,F)==eye(5),'mixed profile inverse')
        req(all(F[j][k]==Q(j==k) for j in (2,3,4) for k in range(5)),'mixed profile changes radial/time coordinate')
        records.append({'vertices':cell['vertices'],'J_affine_coefficients':cell['affine_coefficients'],
                        'forward_homogeneous_matrix':[[str(z) for z in row] for row in F],
                        'inverse_homogeneous_matrix':[[str(z) for z in row] for row in G],
                        'jacobian_determinant':1})
    # Pin the collar enlargement to actual cells.  Filling the puncture back
    # into the compact saved N stage leaves precisely the outer grid torus.
    tup=lambda x:tuple(tup(z) for z in x) if isinstance(x,(tuple,list)) else x
    tets=list(map(tup,M['stages']['N']));removed=[tuple(('N','o',v) for v in t) for t in cut['puncture']['removed_N_tetrahedra']]
    keys={frozenset(t) for t in tets};req(not any(frozenset(t) in keys for t in removed),'puncture was not removed')
    full=tets+removed;inc=Counter(frozenset(f) for t in full for f in combinations(t,3));req(all(n in (1,2) for n in inc.values()),'filled enlarged tube face incidence')
    bd={f for f,n in inc.items() if n==1};grid={frozenset(('N','g',*v) for v in f['ids']) for f in A['G']}
    req(bd==grid,'filled enlarged tube boundary is not the saved outer grid')
    sphere={frozenset(('N','o',v) for v in f) for f in cut['puncture']['boundary_triangles']}
    req(all(inc[f]==2 for f in sphere),'puncture sphere is not internal after collar enlargement')
    sphere_vertices={v for f in sphere for v in f};outer_vertices={v for f in bd for v in f}
    req(not sphere_vertices&outer_vertices,'puncture meets outer grid boundary')
    original_torus={v for f in cut['T'] for v in f};equator=cut['puncture']['equator_cycle'];req(set(equator)<=original_torus,'original puncture/torus contact record')
    return {'schema':'source-bound-mixed-profile-transfer-v1','passed':True,'roll':1,'twist':2,
            'central_certificate_sha256':dig(C),'adaptor_model_sha256':dig(M),'source_cut_sha256':dig(cut),
            'coordinate_order':['u','v_upstairs','r','t','1'],'source_profile':'H0=r/3','central_profile':'H1=clamp(r-1,0,1)',
            'source_monodromy':'(u,v,r)->(u+H0(r),v+H0(r),r)',
            'central_monodromy':'(u,v,r)->(u+H1(r),v+H1(r),r)',
            'map':'Psi_t(u,v,r)=(u+J(r,t),v+J(r,t),r)',
            'seam':'F_central o Psi_1 = Psi_0 o F_source, since J(r,1)+H1(r)=H0(r)',
            'cells':records,'joint_PL':True,'deck_equivariant':True,'time_0_identity':True,
            'puncture_pointwise_fixed':'Psi fixes original N and all of original torus r=0; the actual puncture ball lies in their union',
            'source_winding_upstairs':[1,1],'source_winding_downstairs':[1,2],
            'enlarged_N':{'definition':'original N plus its saved N-adaptor, with the puncture temporarily filled back in; relative coarsening extends by identity over the puncture',
                'tetrahedra':len(full),'outer_boundary_triangles':len(bd),'puncture_boundary_triangles':len(sphere),
                'original_puncture_equator_on_original_torus':True,'puncture_disjoint_from_outer_grid':True,
                'all_puncture_boundary_faces_internal_after_filling':True},
            'composition':'first Psi to the central mixed profile; then K using the framed chart of enlarged N; no commutation of these two global transfers is asserted',
            'scope':'finite source-profile change and actual collar-enlargement incidences; existence of the framed global N-plus chart remains a PL theorem dependency',
            'seconds':time.monotonic()-start}

if __name__=='__main__':
    ans=verify();(ROOT/'central_roll/mixed_profile_transfer.json').write_text(json.dumps(ans,indent=2)+'\n')
    print(json.dumps({k:v for k,v in ans.items() if k!='cells'},indent=2),flush=True)

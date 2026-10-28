#!/usr/bin/env python3
"""Independent exact check of the central roll collar and source transfer.

Does not import the central or adaptor producers.  Polygon intersection,
affine prism and ball-replay checks are separate from the producing routines.
The actual assembled 4D facet membership remains an explicit later gate.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from itertools import combinations
import argparse,copy,gzip,hashlib,json,sys,time

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from verify_adaptor_mesh import overlay_check,half_check,same_coarse,grid,det3

def req(ok,msg):
    if not ok:raise ValueError(msg)
def tup(x):return tuple(tup(z) for z in x) if isinstance(x,(tuple,list)) else x
def load(p):
    p=Path(p);b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def dig(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def mul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def affine3(P,co):
    p=co[0];cols=[tuple(q[k]-p[k] for k in range(3)) for q in co[1:]];D=det3(cols)
    req(D!=0,'singular carrier');rhs=tuple(P[k]-p[k] for k in range(3));a=[]
    for k in range(3):z=cols.copy();z[k]=rhs;a.append(det3(z)/D)
    return [1-sum(a),*a]
def path_product(base,j,k):
    return [tuple([*((base[i],j) for i in range(r+1)),*((base[i],k) for i in range(r,len(base)))]) for r in range(len(base))]
def centtag(v,j):
    if v[0]=='a':return ('B',('N','g',v[1],j),v[2])
    if v[0]=='b':return ('B',('E','g',v[1],j),v[2])
    return ('C',v,j)

def source_profile_check(rec):
    req(rec['roll']==1 and rec['source_radial_coordinate']=='r=3*s_input','source profile parameter')
    req(rec['map']=='(u,v,r,t)->(u+J(r,t),v,r,t)' and rec['inverse']=='(u,v,r,t)->(u-J(r,t),v,r,t)','source transfer formulas')
    expected=[]
    for i in range(3):
        expected.extend([((Q(i),Q(0)),(Q(i+1),Q(0)),(Q(i+1),Q(1))),
                         ((Q(i),Q(0)),(Q(i+1),Q(1)),(Q(i),Q(1)))])
    req(len(rec['J_triangles'])==6,'source transfer triangle count')
    edges=defaultdict(list);area=Q(0)
    for cell,ps in zip(rec['J_triangles'],expected):
        got=tuple(tuple(map(Q,p)) for p in cell['vertices']);req(got==ps,'source transfer rectangle subdivision')
        a,b,c=map(Q,cell['affine_coefficients']);values=list(map(Q,cell['values']))
        req(len(values)==3,'source transfer values count')
        for (r,t),val in zip(ps,values):
            want=Q(0) if t==0 else r/3-min(Q(1),max(Q(0),r-1))
            req(val==want==a*r+b*t+c,'wrong affine profile transfer at vertex')
        for p,q in combinations(ps,2):edges[tuple(sorted((p,q)))].append((a,b,c))
        (r0,t0),(r1,t1),(r2,t2)=ps;area+=((r1-r0)*(t2-t0)-(r2-r0)*(t1-t0))/2
        F=eye(5);G=eye(5);F[0][2:5]=[a,b,c];G[0][2:5]=[-a,-b,-c]
        req(mul(F,G)==mul(G,F)==eye(5),'source transfer inverse')
    req(area==3 and all(len(v) in (1,2) for v in edges.values()),'source transfer domain coverage')
    for edge,aa in edges.items():
        if len(aa)==2:
            for r,t in edge:req(len({a*r+b*t+c for a,b,c in aa})==1,'source transfer discontinuity')
        else:
            a,b,c=aa[0]
            for r,t in edge:
                if r in (0,3) or t==0:req(a*r+b*t+c==0,'source transfer boundary not fixed')
                if t==1:req(a*r+b*t+c+min(Q(1),max(Q(0),r-1))==r/3,'source transfer seam equation')
    req(rec['cap_and_puncture_boundary_pointwise_fixed'] is True,'source transfer support')
    return {'triangles':6,'domain_exact_area':'3','affine_inverses':6,'all_shared_edges_compatible':True,
            'boundary_and_seam_equations':'identities on affine boundary segments, not samples'}

def straightening_check(rec):
    req(rec['coordinate_order']==['u','v','h','t'],'straightening coordinate order')
    req(rec['map']=='(u,v,h,t)->(u+min(h,t),v,h,t)' and rec['inverse']=='(U,v,h,t)->(U-min(h,t),v,h,t)','straightening formula')
    P=[[[0,0],[1,0],[1,1]],[[0,0],[1,1],[0,1]]]
    req(len(rec['cells'])==2,'straightening chart count')
    for k,c in enumerate(rec['cells']):
        req([[Q(z) for z in p] for p in c['h_t_vertices']]==P[k],'straightening triangles')
        A=[[Q(v) for v in row] for row in c['forward_matrix']];B=[[Q(v) for v in row] for row in c['inverse_matrix']]
        want=eye(4);want[0][3 if k==0 else 2]=1
        req(A==want and mul(A,B)==mul(B,A)==eye(4),'straightening affine matrix or inverse')
    req(rec['inner_attachment']=='u_N=U at h=0' and rec['outer_attachment']=='u_E=U-t at h=1','straightening boundary charts')
    req(rec['time_seam']=='(u,v,h,1)~(u+h,v,h,0)' and rec['full_roll_winding']==1,'integer roll seam')
    # The two affine corrections t and h agree on h=t; they vanish at h=0,
    # equal t at h=1, and their top-minus-bottom value is identically h.
    return {'affine_pieces':2,'jacobian_determinants':[1,1],'inverse_compositions':'identity',
            'boundary_and_seam_equations':'exact affine identities','full_integer_winding':1}

def branching(tets):
    ed={};seen=set()
    for t in tets:
        req(len(set(t))==4 and frozenset(t) not in seen,'folded or repeated central tetrahedron');seen.add(frozenset(t))
        for a,b in combinations(t,2):
            e=frozenset((a,b));req(e not in ed or ed[e]==(a,b),'inconsistent central branching');ed[e]=(a,b)

def interface_check(C,adaptor):
    out={}
    for side,letter in [('N','a'),('E','b')]:
        # Direct product of the actual saved grid boundary with the time circle.
        counts=Counter(frozenset(f) for t in map(tup,adaptor['stages'][side]) for f in combinations(t,3))
        boundary={f for f,n in counts.items() if n==1 and all(v[1]=='g' for v in f)}
        G=grid(3,4);req(boundary=={frozenset((side,'g',*v) for v in f['ids']) for f in G},'actual grid endpoint')
        left=set()
        for f in G:
            base=tuple((side,'g',*v) for v in f['order'])
            for j in range(4):
                left.update(tuple(('B',v,k) for v,k in p) for p in path_product(base,j,(j+1)%4))
        right=set()
        for f in C['A' if letter=='a' else 'B']:
            base=tuple((letter,*v) for v in f['order'])
            for j in range(4):right.update(tuple(centtag(v,k) for v,k in p) for p in path_product(base,j,(j+1)%4))
        req(left==right,'central/adaptor product triangulations do not match as ordered tetrahedra')
        out[side]=len(left)
    return out

def witness_check(C,adaptor,coarse):
    endpoint={s:{frozenset(tup(t)):tup(t) for t in ts} for s,ts in adaptor['stages'].items()}
    chartbytet={tup(z['vertices']):z for z in C['universal_charts']};nchecked=0
    wits=C['source_seam_witnesses'];req(len(wits)==5,'source winding witness count')
    for entry,s in zip(wits,(Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1))):
        r=3*s;H=min(Q(1),max(Q(0),r-1));u,v=Q(1,17),Q(1,7)
        req(Q(entry['input_s'])==s and Q(entry['input_u'])==u and Q(entry['input_v'])==v,'source witness coordinates')
        req(list(map(Q,entry['source_coordinate_seam_displacement']))==[s,0] and Q(entry['full_unwrapped_source_winding'])==s,'integer source winding witness')
        req(Q(entry['seam_lattice_residue'])==(-1 if r>=2 else 0),'source seam lattice residue')
        req(len(entry['carriers'])==2,'source carrier pair')
        for cc,tm in zip(entry['carriers'],(1,0)):
            up=u+s-(H if tm else 0);req(cc['time']==tm and Q(cc['source_u_after_profile_transfer'])==up,'profile-transferred source point')
            base=tup(cc['base_tetrahedron']);lam=list(map(Q,cc['base_barycentric_weights']));ww=list(map(Q,cc['product_barycentric_weights']))
            req(len(lam)==4 and min(lam)>=0 and sum(lam)==1 and len(ww)==5 and min(ww)>=0 and sum(ww)==1,'carrier barycentric simplex')
            pv=tup(cc['expected_product_carrier'])
            if 1<r<2:
                req(cc['part']=='central' and base in chartbytet,'central witness parent')
                co=[tuple(map(Q,p)) for p in chartbytet[base]['coords']];U=up+min(H,Q(tm));sh=cc['longitude_lattice_translation']
                req(all(sum(lam[i]*co[i][k] for i in range(4))==x for k,x in enumerate((U+sh,Q(tm),H))),'central witness affine carrier')
                req(list(map(Q,cc['coordinates']))==[U,tm,H] and Q(cc['meridian_coordinate'])==v,'central witness coordinates')
                j=int(4*v);eta=4*v-j;possible=[]
                for k,p in enumerate(path_product(base,j,(j+1)%4)):
                    top=eta-sum(lam[k+1:])
                    if 0<=top<=lam[k]:possible.append((tuple(centtag(a,b) for a,b in p),lam[:k]+[lam[k]-top,top]+lam[k+1:]))
                req((pv,ww) in possible,'central witness product carrier')
            else:
                side='N' if r<=1 else 'E';req(cc['part']==side and endpoint[side].get(frozenset(base))==base,'adapter witness endpoint order')
                old=cc['source_adaptor_witness'];req(old in coarse['protected_source_witnesses'],'unbound adapter witness')
                req(old['roll']==1 and Q(old['input_s'])==s and old['time']==tm,'wrong protected adapter witness')
                weights=dict(zip(map(tup,old['vertices']),map(Q,old['barycentric_weights'])))
                req(lam==[weights[z] for z in base],'reordered adapter witness weights')
                co=[tuple(map(Q,p)) for p in old['coordinates']];oldw=list(map(Q,old['barycentric_weights']));sh=old['lattice_translation']
                req(all(sum(oldw[i]*co[i][k] for i in range(4))==x for k,x in enumerate((up+sh[0],v+sh[1],r))),'adapter witness point')
                want=tuple([('B',base[0],3)]+[('B',a,0) for a in base]) if tm else tuple([('B',a,0) for a in base]+[('B',base[-1],1)])
                req(pv==want and ww==([Q(0)]+lam if tm else lam+[Q(0)]),'adapter witness time product')
            nchecked+=1
    return {'source_levels':['0','1/4','1/2','3/4','1'],'exact_carriers':nchecked,
            'quarter_displacement':['1/4','0'],'full_unwrapped_displacement':['1','0'],
            'actual_4D_facet_membership':'pending bulk assembler replay; expected carriers verified against the product recipe'}

def verify(C,model,adaptor,coarse,cut,chart):
    req(C['schema']=='source-bound-even-positive-roll-central-collar-v1','central schema')
    req((C['L'],C['V'],C['T'],C['roll'])==(3,4,4,1),'central grid/roll scope')
    bindings={'source_cut_sha256':dig(cut),'diagram_chart_sha256':dig(chart),
              'adaptor_model_sha256':dig(adaptor),'adaptor_coarsening_certificate_sha256':dig(coarse)}
    req(C['source_bindings']==model['source_bindings']==bindings,'central source bindings')
    req(model['central_certificate_sha256']==dig(C),'central certificate hash')
    req(C['diagram_sign_transport']==chart['diagram_sign_transport'] and C['coordinate_meaning']==chart['coordinate_meaning'],'corrected source chart')
    same_coarse(C['A'],grid(3,4));same_coarse(C['B'],grid(3,4,1));oc=overlay_check(C['A'],C['B'],C['overlay'])
    aa=half_check(C['A'],C['overlay'],oc,'a',C['charts_a'],lambda z:('a',*z),lambda p:('r',str(p[0]),str(p[1])),lambda i:('ia',i))
    bb=half_check(C['B'],C['overlay'],oc,'b',C['charts_b'],lambda z:('b',*z),lambda p:('r',str(p[0]),str(p[1])),lambda i:('ib',i))
    ts=aa+bb;req(tup(C['tetrahedra'])==tup(model['tetrahedra'])==tuple(ts),'central ordered tetrahedra');branching(ts)
    faces=Counter(frozenset(f) for t in ts for f in combinations(t,3))
    want={frozenset(('a',*v) for v in f['ids']) for f in C['A']}|{frozenset(('b',*v) for v in f['ids']) for f in C['B']}
    req(all(n in (1,2) for n in faces.values()) and {f for f,n in faces.items() if n==1}==want,'central full face incidence/boundary')
    req(len(C['universal_charts'])==len(ts),'physical chart count');volume=Q(0)
    for out,(half,cell) in zip(C['universal_charts'],[('a',c) for c in C['charts_a']]+[('b',c) for c in C['charts_b']]):
        co=[(Q(p[0]),Q(p[1]),Q(p[2])/2 if half=='a' else 1-Q(p[2])/2) for p in cell['coords']]
        req(out['half']==half and out['parent']==cell['parent'] and out['vertices']==cell['vertices'],'physical chart carrier')
        req([tuple(map(Q,p)) for p in out['coords']]==co,'physical radial rescale')
        M=[[co[j][i] for j in range(4)] for i in range(3)]+[[Q(1)]*4];I=[[Q(v) for v in row] for row in out['barycentric_inverse']]
        req(mul(I,M)==mul(M,I)==eye(4),'physical chart barycentric inverse')
        p=co[0];d=det3([tuple(q[k]-p[k] for k in range(3)) for q in co[1:]])
        req(d!=0 and Q(out['determinant'])==d,'physical chart determinant');volume+=abs(d)/6
    req(volume==1,'physical central torus collar volume')
    req(model['deck_on_meridian_indices']==[2,3,0,1] and C['deck_rule']=='central_vertex unchanged; jv -> (jv+2) mod 4','central deck action')
    req(C['time_seam']=='identity on N and E products and on straightened central coordinates; no meridian shift','even time seam')
    interfaces=interface_check(C,adaptor);transfer=source_profile_check(C['source_profile_transfer']);straight=straightening_check(C['straightening'])
    witnesses=witness_check(C,adaptor,coarse)
    return {'schema':'independent-central-roll-geometry-check-v1','passed':True,'roll':1,'tetrahedra':len(ts),
            'overlay_polygons':oc['polygons'],'overlay_triangles':oc['refined_triangles'],
            'physical_exact_volume':str(volume),'boundary_triangles':len(want),
            'ordered_adaptor_product_interface_tetrahedra':interfaces,
            'source_profile_transfer':transfer,'straightening':straight,'source_seam_witnesses':witnesses,
            'deck':'passive four-edge meridian factor, shift 2; all cyclic edge products preserved',
            'model_sha256':dig(model),'certificate_sha256':dig(C),
            'scope':'central collar and source-relative roll geometry only; no proof here of global twist localization, X2 standardness or X0 construction'}

def mutation_checks(C):
    cases=[]
    bad=copy.deepcopy(C['source_profile_transfer']);bad['J_triangles'][2]['affine_coefficients'][1]='1/3'
    try:source_profile_check(bad)
    except ValueError:cases.append('wrong source-profile affine transfer rejected')
    else:raise ValueError('profile mutation escaped')
    bad=copy.deepcopy(C['straightening']);bad['cells'][0]['inverse_matrix'][0][3]=1
    try:straightening_check(bad)
    except ValueError:cases.append('wrong straightening inverse rejected')
    else:raise ValueError('inverse mutation escaped')
    bad=copy.deepcopy(C['B']);bad[0]['points'][2][0]='1/3'
    try:same_coarse(bad,grid(3,4,1))
    except ValueError:cases.append('outer shear reduced to zero rejected')
    else:raise ValueError('winding mutation escaped')
    return cases

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mutations',action='store_true');args=p.parse_args();start=time.monotonic()
    C=load(ROOT/'central_roll/certificate.json.gz');model=load(ROOT/'central_roll/model.json.gz')
    ans=verify(C,model,load(ROOT/'adaptor/model.json.gz'),load(ROOT/'adaptor/coarsening_certificate.json.gz'),load(ROOT/'peripheral/cut_source.json.gz'),load(ROOT/'peripheral/diagram_chart.json.gz'))
    if args.mutations:ans['mutation_checks']=mutation_checks(C)
    ans['seconds']=time.monotonic()-start;(ROOT/'central_roll/independent_verification.json').write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2),flush=True)

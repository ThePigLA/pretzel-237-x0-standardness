#!/usr/bin/env python3
"""Exact marked adapters from the actual source torus to a 3x4 flat grid.

Builds only three-dimensional collars and optional certified relative
coarsenings. The original cut and raw/corrected peripheral charts are not
modified. The two pieces have identical collar atlases but separate tags.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from itertools import combinations
import argparse,gzip,hashlib,json,math,sys,time

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from exact_flat_mesh import grid_faces,overlay,cylinder_half,mod
from verify_adaptor_mesh import overlay_check,half_check,same_coarse,grid
from compact3 import Complex,branching,bary,tup
from verify_compact3 import replay

def load(p):
    p=Path(p);b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)

def convert(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,set):return sorted((convert(v) for v in x),key=repr)
    if isinstance(x,(tuple,list)):return [convert(v) for v in x]
    if isinstance(x,dict):return {str(k):convert(v) for k,v in x.items() if k!='key'}
    return x

def digest(x):return hashlib.sha256(json.dumps(convert(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    b=(json.dumps(convert(x),sort_keys=True,separators=(',',':'))+'\n').encode()
    p.write_bytes(gzip.compress(b,mtime=0) if p.suffix=='.gz' else b)

def check_branching(tets):
    edges={};seen=set()
    for t in tets:
        if len(set(t))!=4 or frozenset(t) in seen:raise ValueError('duplicate or folded tetrahedron')
        seen.add(frozenset(t))
        for a,b in combinations(t,2):
            e=frozenset((a,b))
            if e in edges and edges[e]!=(a,b):raise ValueError('incompatible branching')
            edges[e]=(a,b)

def make_adaptor(L=3,V=4):
    if (L,V)!=(3,4):raise ValueError('only declared 3x4 grid admitted')
    cut=load(ROOT/'peripheral/cut_source.json.gz');chart=load(ROOT/'peripheral/diagram_chart.json.gz')
    if chart['cut_sha256']!=digest(cut):raise ValueError('corrected chart cut source hash')
    xy={v:(Q(x),Q(y)) for v,x,y in chart['upstairs_vertex_coordinates']};D=cut['deck']
    per={tuple(e):tuple(p) for e,p in zip(chart['torus_upstairs']['edge_representatives'],chart['upstairs_periods'])}
    rank=lambda v:min(v,D[v])
    A=[]
    for t in cut['T']:
        if len({rank(v) for v in t})!=3:raise ValueError('source torus triangle repeats deck orbit')
        a=t[0];coords=[]
        for v in t:
            p=(0,0) if v==a else per[a,v]
            coords.append(tuple(xy[v][k]+p[k] for k in (0,1)))
        A.append({'ids':tuple(t),'points':tuple(coords),'order':tuple(sorted(t,key=rank))})
    G=grid_faces(L,V);O=overlay(A,G,Q(1,2))
    ai={frozenset(f['ids']):i for i,f in enumerate(A)};gi={frozenset(f['ids']):i for i,f in enumerate(G)}
    da=[ai[frozenset(D[v] for v in f['ids'])] for f in A]
    dg=[gi[frozenset((i,(j+V//2)%V) for i,j in f['ids'])] for f in G]
    e=tuple(cut['puncture']['fixed_edge_cut_vertices']);removed=list(map(tuple,cut['puncture']['removed_N_tetrahedra']))
    if {t for t in map(tuple,cut['N']) if set(e)<=set(t)}!=set(removed):raise ValueError('puncture source star mismatch')
    bd=list(map(tuple,cut['puncture']['boundary_triangles']));stages={};deck={};charts={}
    for side in ('N','E'):
        cv=lambda v:(side,'o',v);rv=lambda p:(side,'r',*p);gv=lambda q:(side,'g',*q)
        ca=lambda i:(side,'a',i);cb=lambda i:(side,'b',i)
        orig=[tuple(cv(v) for v in sorted(t,key=rank)) for t in cut[side] if side!='N' or tuple(t) not in set(removed)]
        ta,cha=cylinder_half(A,O,'a',cv,rv,ca);tb,chb=cylinder_half(G,O,'b',gv,rv,cb)
        stages[side]=orig+ta+tb;check_branching(stages[side]);charts[side]={'original':orig,'a':cha,'b':chb}
        for v in {v for t in stages[side] for v in t}:
            if v[1]=='o':w=cv(D[v[2]])
            elif v[1]=='r':w=rv(mod((v[2],v[3]+Q(1,2))))
            elif v[1]=='g':w=gv((v[2],(v[3]+V//2)%V))
            elif v[1]=='a':w=ca(da[v[2]])
            else:w=cb(dg[v[2]])
            deck[v]=w
        ts={frozenset(t):t for t in stages[side]}
        for t in stages[side]:
            im=tuple(deck[v] for v in t)
            if im!=ts.get(frozenset(im)):raise ValueError('deck does not preserve simplex order')
        counts=Counter(frozenset(f) for t in stages[side] for f in combinations(t,3))
        boundary={f for f,c in counts.items() if c==1}
        want={frozenset(gv(v) for v in f['ids']) for f in G}
        if side=='N':want|={frozenset(cv(v) for v in f) for f in bd}
        if boundary!=want or any(c not in (1,2) for c in counts.values()):raise ValueError('adapter boundary mismatch')
    ac={'schema':'marked-source-rectangular-torus-adaptor-v1','A':A,'G':G,'overlay':O,'charts':charts,
        'cut_sha256':digest(cut),'chart_sha256':digest(chart),'puncture_edge':e,'removed_tets':removed,'boundary':bd,
        'L':L,'V':V,'source_order':'increasing min(cut_vertex,deck(cut_vertex))',
        'chart_meaning':chart['coordinate_meaning'],'diagram_sign_transport':chart['diagram_sign_transport'],
        'radial_coordinates':{'N_a':'r=h/2','N_b':'r=1-h/2','E_a':'r=3-h/2','E_b':'r=2+h/2'},
        'inverse':'on each recorded affine tetrahedron, exact barycentric coordinates recover the same (u,v,h); common-chart coordinate is unchanged modulo Z2',
        'scope':'3D source-bound adapters only; no mapping torus or four-dimensional cap'}
    return stages,deck,ac,cut

def verify_adaptor(stages,deck,ac,cut):
    # Serialize before checking: the independent checker reads the same exact
    # rational data as a future replay, not in-memory producer closures.
    record=convert(ac);chart=load(ROOT/'peripheral/diagram_chart.json.gz')
    if ac['cut_sha256']!=digest(cut) or ac['chart_sha256']!=digest(chart):raise ValueError('adapter source hashes')
    xy={v:(Q(x),Q(y)) for v,x,y in chart['upstairs_vertex_coordinates']}
    E={tuple(e):tuple(p) for e,p in zip(chart['torus_upstairs']['edge_representatives'],chart['upstairs_periods'])}
    actual=[]
    for tri in cut['T']:
        a=tri[0];co=[]
        for v in tri:
            period=(0,0) if v==a else E[a,v]
            co.append(tuple(xy[v][k]+period[k] for k in (0,1)))
        actual.append({'ids':tuple(tri),'points':tuple(co),'order':tuple(sorted(tri,key=lambda v:min(v,cut['deck'][v])))})
    same_coarse(record['A'],actual);same_coarse(record['G'],grid(3,4))
    oc=overlay_check(record['A'],record['G'],record['overlay'],Q(1,2))
    wanted={};D=cut['deck'];rank=lambda v:min(v,D[v]);removed=set(map(tuple,cut['puncture']['removed_N_tetrahedra']))
    for side in ('N','E'):
        cv=lambda v:(side,'o',v);rv=lambda p:(side,'r',str(p[0]),str(p[1]));gv=lambda q:(side,'g',*q)
        orig=[tuple(cv(v) for v in sorted(t,key=rank)) for t in cut[side] if side!='N' or tuple(t) not in removed]
        if tup(record['charts'][side]['original'])!=tuple(orig):raise ValueError('serialized source altered')
        aa=half_check(record['A'],record['overlay'],oc,'a',record['charts'][side]['a'],cv,rv,lambda i:(side,'a',i))
        bb=half_check(record['G'],record['overlay'],oc,'b',record['charts'][side]['b'],gv,rv,lambda i:(side,'b',i))
        wanted[side]=orig+aa+bb
        if wanted[side]!=list(map(tup,stages[side])):raise ValueError('adapter tetrahedra differ from exact checked geometry')
    return {'schema':'exact-source-adaptor-geometry-verification-v1','passed':True,
            'source_torus_triangles':len(ac['A']),'grid_triangles':len(ac['G']),
            'overlay_polygons':oc['polygons'],'overlay_triangles':oc['refined_triangles'],
            'tetrahedra_before_coarsening':{s:len(ts) for s,ts in stages.items()},
            'each_half_collar_exact_volume':'1','inverse_scope':'nonzero exact affine determinants on all child tetrahedra, exact face compatibility and nonoverlapping complete parent partitions',
            'source_sha256':ac['cut_sha256'],'chart_sha256':ac['chart_sha256']}

def support_tags(ac,chart,deck):
    witnesses=[];protect=set()
    for n in (0,1):
        for s in (Q(0),Q(1,4),Q(3,4),Q(1)):
            r=3*s;side='N' if r<=1 else 'E';h=min(Q(1),max(Q(0),r-1));u,v=Q(1,17),Q(1,7)
            for tm,x in ((1,u+n*(s-h)),(0,u+n*s)):
                found=None
                for half in ('a','b'):
                    for cell in ac['charts'][side][half]:
                        co=[tuple(map(Q,p)) for p in cell['coords']]
                        if side=='N':co=[(p[0],p[1],p[2]/2 if half=='a' else 1-p[2]/2) for p in co]
                        else:co=[(p[0],p[1],3-p[2]/2 if half=='a' else 2+p[2]/2) for p in co]
                        if not min(z[2] for z in co)<=r<=max(z[2] for z in co):continue
                        ir=range(math.ceil(min(z[0] for z in co)-x),math.floor(max(z[0] for z in co)-x)+1)
                        jr=range(math.ceil(min(z[1] for z in co)-v),math.floor(max(z[1] for z in co)-v)+1)
                        for i in ir:
                            for j in jr:
                                w=bary((x+i,v+j,r),co)
                                if w is not None:
                                    found={'roll':n,'input_s':str(s),'time':tm,'side':side,'half':half,
                                           'vertices':tup(cell['vertices']),'coordinates':[[str(z) for z in p] for p in co],
                                           'barycentric_weights':list(map(str,w)),'lattice_translation':[i,j]};break
                            if found:break
                        if found:break
                    if found:break
                if not found:raise ValueError('unlocated exact source-quarter witness')
                witnesses.append(found);protect.update(found['vertices'])
    per={('E','o',chart['torus_downstairs']['vertices'][0])}
    for loop in chart['peripheral_loops'].values():
        for e in loop['torus_edge_loop']:
            per.update(('E','o',v) for v in chart['torus_downstairs']['edge_representatives'][abs(e)-1])
    protect|=per;protect|={deck[v] for v in list(protect)}
    return protect,witnesses,sorted(per,key=repr)

def coarsen_adapter(stages,deck,ac,chart,seconds=120):
    stages={s:list(map(tup,ts)) for s,ts in stages.items()};deck={tup(a):tup(b) for a,b in deck.items()}
    protected,wits,per=support_tags(ac,chart,deck)
    cert={'schema':'equivariant-relative-3D-coarsening-v1','parts':{},'protected_source_witnesses':wits,
          'protected_source_peripheral_vertices':per,
          'topology_rule':'interior-star ball replacements with identical boundary; fixed arc, puncture sphere, grid and selected source paths/witnesses fixed'}
    for side,ts in stages.items():
        tags=sorted({v for t in ts for v in t},key=repr);ix={v:i for i,v in enumerate(tags)};D=[ix[deck[v]] for v in tags]
        initial=[tuple(ix[v] for v in t) for t in ts];K={tuple(sorted(t)) for t in initial}
        faces=Counter(f for t in K for f in combinations(t,3));bd={f for f,c in faces.items() if c==1}
        frozen={v for f in bd for v in f}|{i for i,d in enumerate(D) if i==d}|{ix[v] for v in protected if v in ix};frozen|={D[v] for v in list(frozen)}
        worker=Complex(K,D,frozen);events=worker.run(seconds=seconds);orders,stats=branching(worker.K,initial,D,frozen)
        cert['parts'][side]={'tags':tags,'deck':D,'protected':sorted(frozen),'initial':sorted(K),'initial_ordered':initial,
                            'initial_sha256':digest(sorted(K)),'endpoint':sorted(worker.K),'endpoint_sha256':digest(sorted(worker.K)),
                            'ordered_endpoint':orders,'steps':events,'branching_search':stats}
        print('ADAPTER COARSEN',side,len(K),'->',len(worker.K),'pairs',len(events),flush=True)
    compact,deck,receipt=replay(convert(cert),stages,deck,convert(ac),chart)
    return compact,deck,cert,receipt

def build_bundle(output=ROOT/'adaptor',coarsen=True,seconds=120):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if (output/'model.json.gz').exists():raise FileExistsError('use a fresh output directory')
    start=time.monotonic();stages,deck,ac,cut=make_adaptor();chart=load(ROOT/'peripheral/diagram_chart.json.gz')
    print('ADAPTER BUILT', {s:len(t) for s,t in stages.items()},flush=True)
    geometry=verify_adaptor(stages,deck,ac,cut);save(output/'geometry_verified.json',geometry)
    save(output/'adaptor.json.gz',ac);save(output/'initial_stages.json.gz',{'stages':stages,'deck':list(deck.items())})
    if coarsen:
        stages,deck,cert,receipt=coarsen_adapter(stages,deck,ac,chart,seconds)
        save(output/'coarsening_certificate.json.gz',cert);save(output/'coarsening_verified.json',receipt)
    model={'schema':'marked-pretzel-rectangular-adaptor-model-v1','L':3,'V':4,'stages':stages,'deck':list(deck.items()),
           'puncture_boundary':[[('N','o',v) for v in sorted(f,key=lambda v:min(v,cut['deck'][v]))] for f in ac['boundary']],
           'source_cut_sha256':digest(cut),'diagram_chart_sha256':digest(chart),'adaptor_sha256':digest(ac),
           'coarsened':coarsen,'coarsening_certificate_sha256':digest(cert) if coarsen else None,
           'source_orientation':'J-oriented source; corrected chart is (u_old,-v_old); full source identities retained'}
    save(output/'model.json.gz',model)
    result={'passed':True,'tetrahedra':{s:len(ts) for s,ts in stages.items()},'coarsened':coarsen,
            'naive_four_time_product_pentachora':16*sum(map(len,stages.values())),
            'model_sha256':digest(model),'seconds':time.monotonic()-start}
    save(output/'summary.json',result);print(json.dumps(result,indent=2),flush=True)
    return stages,deck,ac,cut

def read_bundle(output=ROOT/'adaptor'):
    output=Path(output);x=load(output/'model.json.gz');ac=load(output/'adaptor.json.gz');cut=load(ROOT/'peripheral/cut_source.json.gz')
    return {s:list(map(tup,ts)) for s,ts in x['stages'].items()},{tup(a):tup(b) for a,b in x['deck']},ac,cut

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'adaptor');p.add_argument('--no-coarsen',action='store_true');p.add_argument('--seconds',type=int,default=120)
    args=p.parse_args();build_bundle(args.output,not args.no_coarsen,args.seconds)

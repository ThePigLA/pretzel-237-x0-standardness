#!/usr/bin/env python3
"""Assemble marked X0 from the identical audited X2 bulk and the audited m0 cap.

This module never calls or weakens the old odd-twist builder. It requires the
actual pretzel adaptor, a longitudinal/time shear collar, the source-bound
spherical bridge, and the zero-winding m=0 cap. A source-bound X2 acceptance receipt is required.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter, defaultdict
import argparse, gzip, hashlib, json, sys, time

HERE=Path(__file__).resolve().parents[1]
ROOT=HERE/'frozen_release'
sys.path.insert(0,str(ROOT))
from torus_adaptor import read_bundle, load, save, digest, convert

def tup(x):
    return tuple(map(tup,x)) if isinstance(x,list) else x

def stair(a,b):
    for seq in combinations(range(len(a)+len(b)-2),len(a)-1):
        seq=set(seq);i=j=0;p=[(a[0],b[0])]
        for k in range(len(a)+len(b)-2):
            if k in seq:i+=1
            else:j+=1
            p.append((a[i],b[j]))
        yield tuple(p)

def parity(p):
    return (-1)**sum(p[a]>p[b] for a in range(len(p)) for b in range(a+1,len(p)))

def boundary(fs):
    hits=Counter(f for t in fs for f in combinations(sorted(t),len(t)-1))
    if any(c not in (1,2) for c in hits.values()):raise ValueError('nonmanifold facet incidence')
    return sorted(f for f,c in hits.items() if c==1)

def face_data(facets,closed=True):
    pending={};seen=set();gluings=[];adj=[[] for _ in facets]
    for i,f in enumerate(facets):
        for k in range(5):
            face=f[:k]+f[k+1:]
            if face in seen:raise ValueError('face appears more than twice')
            if face not in pending:pending[face]=(i,k)
            else:
                j,l=pending.pop(face);seen.add(face)
                p=[facets[j].index(v) if a!=k else l for a,v in enumerate(f)]
                gluings.append((i,k,j,p));s=-parity(p)
                adj[i].append((j,s));adj[j].append((i,s))
    if closed and pending:raise ValueError(('unpaired facets',len(pending)))
    signs={0:1};todo=[0]
    for i in todo:
        for j,s in adj[i]:
            q=signs[i]*s
            if j in signs:
                if signs[j]!=q:raise ValueError('nonorientable')
            else:signs[j]=q;todo.append(j)
    if len(signs)!=len(facets):raise ValueError('disconnected')
    return gluings,[signs[i] for i in range(len(facets))],sorted(pending)

def quotient_data(facets,deck,gluings):
    lookup={f:i for i,f in enumerate(facets)}
    cell_deck=[lookup[tuple(sorted(deck[v] for v in f))] for f in facets]
    if any(i==j for i,j in enumerate(cell_deck)):raise ValueError('invariant top simplex')
    reps=[i for i,j in enumerate(cell_deck) if i<j];qi={i:k for k,i in enumerate(reps)}
    face_adj={}
    for i,k,j,p in gluings:
        face_adj[i,k]=(j,p)
        inv=[p.index(a) for a in range(5)];face_adj[j,p[k]]=(i,inv)
    qgl=[]
    for i in reps:
        for k in range(5):
            j,p=face_adj[i,k];r=min(j,cell_deck[j]);q=qi[r]
            pp=[facets[r].index(deck[facets[j][p[a]]] if j!=r else facets[j][p[a]]) for a in range(5)]
            if (qi[i],k)<(q,pp[k]):qgl.append((qi[i],k,q,pp))
    if len(qgl)*2!=5*len(reps):raise ValueError('quotient facet count')
    covering=[]
    for i,f in enumerate(facets):
        r=min(i,cell_deck[i])
        p=[facets[r].index(deck[v] if i!=r else v) for v in f]
        covering.append((qi[r],p))
    return {'schema':'orbit-cell-quotient-v1','simplex_count':len(reps),'gluings':qgl,
            'cover_cell_representatives':reps,'covering_simplex_maps':covering,
            'cell_deck':cell_deck,'vertex_orbit_labels':[min(v,deck[v]) for v in range(len(deck))]}

def surface_check(triangles,closed=True):
    triangles=list(map(tuple,triangles));vs={v for t in triangles for v in t}
    ec=Counter(e for t in triangles for e in combinations(sorted(t),2))
    if any(c not in ((2,) if closed else (1,2)) for c in ec.values()):raise ValueError('surface edge links')
    for v in vs:
        es=[tuple(x for x in t if x!=v) for t in triangles if v in t];g=defaultdict(set)
        for a,b in es:g[a].add(b);g[b].add(a)
        seen={next(iter(g))};todo=list(seen)
        for a in todo:
            for b in g[a]:
                if b not in seen:seen.add(b);todo.append(b)
        if seen!=set(g) or (closed and any(len(n)!=2 for n in g.values())):raise ValueError('surface vertex link')
    g=defaultdict(set)
    for a,b in ec:g[a].add(b);g[b].add(a)
    seen={next(iter(vs))};todo=list(seen)
    for a in todo:
        for b in g[a]:
            if b not in seen:seen.add(b);todo.append(b)
    if seen!=vs:raise ValueError('surface disconnected')
    return [len(vs),len(ec),len(triangles)]

def anchor_orientation(facets,tags,signs,stages,cut):
    """Anchor Y x positively oriented time to an unchanged source tetrahedron."""
    co=load(ROOT/'peripheral/cut_orientation.json');idx={tup(v):i for i,v in enumerate(tags)}
    fmap={tuple(f):i for i,f in enumerate(facets)}
    ci={tuple(t):i for i,t in enumerate(cut['N'])}
    te=next(t for t in stages['N'] if all(v[1]=='o' for v in t) and tuple(sorted(v[2] for v in t)) in ci)
    original=tuple(v[2] for v in te);so=tuple(sorted(original));k=ci[so]
    source_sign=co['N_tetrahedron_orientations'][k]
    base_sign=source_sign*parity([so.index(v) for v in original])
    pp=[('B',v,0) for v in te]+[('B',te[-1],1)]
    ordered=tuple(idx[v] for v in pp);f=tuple(sorted(ordered));i=fmap[f]
    target=base_sign*parity([f.index(v) for v in ordered]);mult=target*signs[i]
    result=[mult*s for s in signs]
    anchor={'schema':'source-oriented-product-anchor-v1','source_ambient':'positive signed-diagram ambient = canonical Seifert master',
            'four_orientation':'source three-orientation followed by positive time',
            'source_N_tetrahedron_index':k,'source_cut_vertices_in_product_order':original,'source_sorted_chain_sign':source_sign,
            'branch_order_sign':base_sign,'time_edge':[0,1],'ordered_product_vertices':ordered,
            'pentachoron_index':i,'canonical_pentachoron':f,'canonical_orientation':target,
            'cut_orientation_sha256':digest(co),'argument':'The last-vertical-edge staircase simplex has determinant +1 in source-then-time product coordinates; propagate across every glued facet.'}
    return result,anchor

def build(max_pents=250000):
    gate_path=HERE/'x2_acceptance.json'
    gate=load(gate_path)
    required='af21a7dd107b765d5fef0a257dcb4df6d6df8689c594f584e4318d15c49aa709'
    if gate.get('gate_open_for_X0') is not True or gate.get('original_master_file_sha256')!=required:
        raise ValueError('Independent source-bound X2 standardization gate is not open')
    if hashlib.sha256((ROOT/'control/marked_control.json.gz').read_bytes()).hexdigest()!=required:
        raise ValueError('Frozen X2 master is not the accepted control')
    start=time.monotonic();out=HERE/'candidate';out.mkdir(exist_ok=False)
    stages,deck3,ac,cut=read_bundle()
    central=load(ROOT/'central_roll/model.json.gz')
    ct=list(map(tup,central['tetrahedra']))
    bridge=load(ROOT/'cap/data/source_cap_bridge.json.gz');cap=load(ROOT/'cap/data/cap_m0.json.gz')
    if cap['twist']!=0 or cap['beta']['cap_full_turns']!=0:raise ValueError('m0 identity cap required')
    L,V,T=3,4,4
    predicted=4*T*sum(map(len,stages.values()))+4*V*len(ct)+len(bridge['facets'])+len(cap['facets'])
    print('PREDICTED',predicted,flush=True)
    if predicted>max_pents:raise RuntimeError('declared construction cap exceeded before allocation')
    parts={k:[] for k in ('N','E','central','bridge','cap')};D={}
    def prod(v,j):return ('B',v,j%T)
    for side in ('N','E'):
        for te in stages[side]:
            for j in range(T):parts[side].extend(tuple(prod(v,h) for v,h in p) for p in stair(te,(j,(j+1)%T)))
        for v in {v for t in stages[side] for v in t}:
            for j in range(T):D[prod(v,j)]=prod(deck3[v],j)
    def cent(v,j):
        j%=V
        if v[0]=='a':return prod(('N','g',v[1],j),v[2])
        if v[0]=='b':return prod(('E','g',v[1],j),v[2])
        return ('C',v,j%V)
    for te in ct:
        for j in range(V):parts['central'].extend(tuple(cent(v,h) for v,h in p) for p in stair(te,(j,(j+1)%V)))
    for v in {v for t in ct for v in t}:
        for j in range(V):D[cent(v,j)]=cent(v,j+V//2)
    bulkbd=boundary(parts['N']+parts['E']+parts['central'])
    sv=bridge['source_boundary']['vertices'];smap={b:prod(('N','o',v),j) for b,(v,j) in zip(bridge['ports']['source']['vertices'],sv)}
    sourceport=sorted(tuple(sorted(smap[v] for v in f)) for f in bridge['ports']['source']['facets'])
    if bulkbd!=sourceport:raise ValueError(('bulk/source port mismatch',len(bulkbd),len(sourceport)))
    cmap={b:('K',i) for i,b in enumerate(bridge['ports']['cap']['vertices'])}
    bm={v:smap[v] if v in smap else cmap[v] if v in cmap else ('H',v) for v in range(len(bridge['vertices']))}
    for v,w in enumerate(bridge['deck']):
        if bm[v] in D and D[bm[v]]!=bm[w]:raise ValueError('bridge deck port mismatch')
        D[bm[v]]=bm[w]
    parts['bridge']=[tuple(bm[v] for v in f) for f in bridge['facets']]
    parts['cap']=[tuple(('K',v) for v in f) for f in cap['facets']]
    if sorted(tuple(sorted(cmap[v] for v in f)) for f in bridge['ports']['cap']['facets'])!=sorted(tuple(('K',v) for v in f) for f in cap['boundary']):raise ValueError('cap port mismatch')
    for v,w in enumerate(cap['deck']):
        if ('K',v) in D and D['K',v]!=('K',w):raise ValueError('cap deck port mismatch')
        D['K',v]=('K',w)
    allord=[f for pp in parts.values() for f in pp];tags=sorted({v for f in allord for v in f},key=repr);idx={v:i for i,v in enumerate(tags)}
    ids=lambda f:tuple(sorted(idx[v] for v in f))
    facets=sorted(map(ids,allord));deck=[idx[D[v]] for v in tags];fset=set(facets)
    if len(facets)!=predicted or len(fset)!=predicted:raise ValueError('folded/duplicate pentachora')
    if any(deck[deck[v]]!=v for v in range(len(deck))):raise ValueError('deck square')
    for f in facets:
        if len(set(f))!=5 or tuple(sorted(deck[v] for v in f)) not in fset:raise ValueError('deck facet mismatch')
        if len({min(v,deck[v]) for v in f})!=5:raise ValueError('simplex repeats a deck orbit')
    gl,sg,bd=face_data(facets)
    sg,orientation_anchor=anchor_orientation(facets,tags,sg,stages,cut)
    f_index={f:i for i,f in enumerate(facets)}
    for i,f in enumerate(facets):
        image=tuple(deck[v] for v in f);j=f_index[tuple(sorted(image))]
        if sg[j]!=sg[i]*parity([facets[j].index(v) for v in image]):raise ValueError('deck reverses orientation')
    fixed={v for v,w in enumerate(deck) if v==w}
    if max(sum(v in fixed for v in f) for f in facets)!=3:raise ValueError('fixed dimension')
    ram=sorted({tuple(v for v in f if v in fixed) for f in facets if sum(v in fixed for v in f)==3})
    twin=sorted(ids([('K',v) for v in f]) for f in cap['twin_triangles'])
    rf=surface_check(ram);tf=surface_check(twin)
    if sum((-1)**i*n for i,n in enumerate(rf))!=2 or sum((-1)**i*n for i,n in enumerate(tf))!=2:raise ValueError('marked sphere Euler characteristic')
    ramset=set(ram);links=defaultdict(list)
    for f in facets:
        r=tuple(v for v in f if v in fixed)
        if len(r)==3:links[r].append(tuple(v for v in f if v not in fixed))
    for r,es in links.items():
        g=defaultdict(set)
        for a,b in es:g[a].add(b);g[b].add(a)
        seen={next(iter(g))};todo=list(seen)
        for a in todo:
            for b in g[a]:
                if b not in seen:seen.add(b);todo.append(b)
        if len(set(es))!=len(es) or seen!=set(g) or any(len(n)!=2 for n in g.values()):raise ValueError('ramification normal link not circle')
        if any(deck[v]==v for v in g) or any(tuple(sorted((deck[a],deck[b])))==tuple(sorted((a,b))) for a,b in es):raise ValueError('normal deck not free')
    q=quotient_data(facets,deck,gl);sha=digest(facets)
    model={'schema':'marked-pretzel-X0-candidate-v1','twist':0,'roll':1,'grid':[L,V,T],
           'status':'ASSEMBLED_PENDING_X0_ADMISSION_AND_INVARIANTS','facets':facets,'facets_sha256':sha,'vertices':tags,'deck':deck,'orientation':sg,
           'ramification_triangles':ram,'twin_triangles':twin,'orientation_anchor':orientation_anchor,'parts':{s:sorted(map(ids,fs)) for s,fs in parts.items()},
           'input_bindings':{'adaptor':digest(load(ROOT/'adaptor/model.json.gz')),'central':digest(central),'bridge':digest(bridge),'cap':digest(cap)},
           'construction_limits':{'max_pentachora':max_pents},'bulk_boundary':list(map(ids,bulkbd)),
           'source_port_map':[[v,idx[bm[v]]] for v in bridge['ports']['source']['vertices']],
           'cap_port_map':[[i,idx[('K',i)]] for i in range(len(bridge['ports']['cap']['vertices']))]}
    rg={'simplex_count':len(facets),'source_facets_sha256':sha,'gluings':gl}
    q['source_facets_sha256']=sha
    save(out/'marked_candidate.json.gz',model);save(out/'cover_regina.json.gz',rg);save(out/'quotient_regina.json.gz',q)
    receipt={'schema':'marked-X0-assembly-receipt-v1','passed':True,'pentachora':len(facets),'vertices':len(tags),
             'parts':{s:len(fs) for s,fs in parts.items()},'quotient_pentachora':q['simplex_count'],
             'ramification_f_vector':rf,'twin_f_vector':tf,'normal_branch_circles_checked':len(links),
             'closed':True,'connected':True,'orientable':True,'deck_orientation_preserving':True,
             'bulk_boundary_tetrahedra':len(bulkbd),'facet_sha256':sha,'seconds':time.monotonic()-start,
             'scope':'Full finite assembly and orbit quotient; native link/invariant admission, global relative-transfer audit and standardization are separate gates.'}
    save(out/'assembly_receipt.json',receipt);print(json.dumps(receipt),flush=True)
    return model,rg,q

if __name__=='__main__':build()

#!/usr/bin/env python3
"""Finite local hypotheses for the source-linked twist-localization lemma.

Checks the actual fixed subcomplex, normal links, branched vertex links,
epsilon-neighborhood cut and framed puncture.  The regular-neighborhood
and relative PL extension theorems are stated explicitly as trust boundaries.
This is not a producer or a global normal-disk-bundle trivialization.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from itertools import combinations
import gzip,hashlib,json,time

ROOT=Path(__file__).resolve().parent
def req(ok,msg):
    if not ok:raise ValueError(msg)
def load(p):
    p=ROOT/p;b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def dig(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def orbit(face,D):
    f=tuple(sorted(face));return min(f,tuple(sorted(D[v] for v in f)))
def cycle(edges):
    edges={tuple(sorted(e)) for e in edges};adj=defaultdict(set)
    for a,b in edges:req(a!=b,'loop in normal link');adj[a].add(b);adj[b].add(a)
    req(adj and all(len(v)==2 for v in adj.values()),'link graph is not 2-regular')
    start=min(adj);out=[start];prev=None;cur=start
    while True:
        nxt=min(adj[cur]-({prev} if prev is not None else set()))
        if nxt==start:break
        req(nxt not in out,'link graph has multiple cycles');out.append(nxt);prev,cur=cur,nxt
    req(len(out)==len(adj)==len(edges),'link graph disconnected');return out
def sphere(tris):
    tris={tuple(sorted(t)) for t in tris};edges=Counter(e for t in tris for e in combinations(t,2));V={v for t in tris for v in t}
    req(all(n==2 for n in edges.values()),'vertex link has boundary or singular edges')
    for v in V:cycle(tuple(x for x in t if x!=v) for t in tris if v in t)
    adj=defaultdict(set)
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    seen={min(V)};todo=list(seen)
    for a in todo:
        for b in adj[a]:
            if b not in seen:seen.add(b);todo.append(b)
    req(seen==V and len(V)-len(edges)+len(tris)==2,'vertex link not a 2-sphere')
    return V,edges,tris
def parity(order):return (-1)**sum(order[i]>order[j] for i in range(len(order)) for j in range(i+1,len(order)))

def verify():
    start=time.monotonic();S=load('coarsened/source.json.gz');C=load('peripheral/cut_source.json.gz');H=load('peripheral/diagram_chart.json.gz');ORI=load('peripheral/cut_orientation.json')
    MIX=load('central_roll/mixed_profile_transfer.json');req(MIX['passed'] and MIX['source_cut_sha256']==dig(C),'collar enlargement/profile receipt')
    K=[tuple(t) for t in S['closed_cover_facets']];D=S['deck'];F={v for v,d in enumerate(D) if v==d}
    req(all(D[D[v]]==v for v in range(len(D))),'source deck is not involution')
    req(all(len({min(v,D[v]) for v in t})==4 for t in K),'source simplex folds a deck orbit')
    KS=set(K);req(all(tuple(sorted(D[v] for v in t)) in KS for t in K),'source deck fails tetrahedra')
    req(all(len(set(t)&F)<=2 for t in K),'fixed locus not one-dimensional')
    fixed_edges={tuple(v for v in t if v in F) for t in K if len(set(t)&F)==2};core=cycle(fixed_edges)
    req(set(core)==F,'fixed subcomplex has isolated vertices');edge_rows=[];vertex_rows=[]
    for e in sorted(fixed_edges):
        le={tuple(v for v in t if v not in e) for t in K if set(e)<=set(t)};cyc=cycle(le);n=len(cyc)
        req(n%2==0 and all(D[v]==cyc[(i+n//2)%n] for i,v in enumerate(cyc)),'normal deck is not a free half-turn')
        req(all(D[v]!=v for v in cyc),'fixed normal direction')
        edge_rows.append({'edge':e,'normal_link_cycle':cyc,'deck_cyclic_shift':n//2,'quotient_circle_cells':[n//2,n//2]})
    for v in sorted(F):
        L=[tuple(x for x in t if x!=v) for t in K if v in t];vv,ee,tt=sphere(L)
        fixed_dirs=vv&F;req(len(fixed_dirs)==2 and fixed_dirs=={x for e in fixed_edges if v in e for x in e if x!=v},'branch vertex directions')
        qv={orbit((x,),D) for x in vv};qe={orbit(e,D) for e in ee};qt={orbit(t,D) for t in tt}
        req(len(qv)-len(qe)+len(qt)==2,'branched quotient link Euler characteristic')
        # Away from the two fixed directions the action is free.  At either
        # fixed direction its link is the already checked free normal circle.
        # Thus the quotient is a connected closed PL surface with chi=2.
        req(len(qt)*2==len(tt) and len(qe)*2==len(ee),'unexpected fixed link cell')
        vertex_rows.append({'vertex':v,'link_f':[len(vv),len(ee),len(tt)],'fixed_directions':sorted(fixed_dirs),
                            'quotient_link_f':[len(qv),len(qe),len(qt)],'quotient_link':'connected closed PL surface, Euler characteristic 2'})
    req(C['source_sha256']==dig(S) and H['cut_sha256']==dig(C),'cut/chart source hash')
    req(C['level']=='sum(fixed-vertex barycentric coordinates)=1/2','regular-neighborhood level')
    req(set(C['fixed_source_vertices'])==F,'cut fixed vertices')
    req(ORI['source_sha256']==dig(S) and ORI['cut_sha256']==dig(C) and ORI['passed'],'cut exact-partition receipt')
    P=C['points'];DP=C['deck'];levels=[]
    for p in P:
        if p[0]=='v':level=Q(p[1] in F)
        else:
            req(p[0]=='m' and len(p)==3,'unknown cut point');level=Q((p[1] in F)+(p[2] in F),2)
        levels.append(level)
    for side in ('N','E'):
        for t in C[side]:req(all(levels[v]>=Q(1,2) for v in t) if side=='N' else all(levels[v]<=Q(1,2) for v in t),'wrong epsilon-neighborhood side')
    punct=C['puncture'];e=punct['fixed_edge_cut_vertices'];raw=punct['equator_cycle'];positive=[raw[0],*reversed(raw[1:])]
    req(all(DP[v]==v for v in e) and len(raw)==4,'puncture axis/circle')
    joins={tuple(sorted((*e,a,b))) for a,b in zip(raw,raw[1:]+raw[:1])}
    req(joins==set(map(tuple,punct['removed_N_tetrahedra']))=={tuple(t) for t in C['N'] if set(e)<=set(t)},'actual puncture is not complete edge join')
    req(all(DP[v]==raw[(i+2)%4] for i,v in enumerate(raw)),'puncture deck half-turn')
    req(set(map(tuple,punct['boundary_triangles']))=={tuple(sorted((p,a,b))) for p in e for a,b in zip(raw,raw[1:]+raw[:1])},'puncture sphere suspension')
    signs=dict(zip(map(tuple,C['N']),ORI['N_tetrahedron_orientations']));normal_signs=[]
    for a,b in zip(positive,positive[1:]+positive[:1]):
        join=(*e,a,b);normal_signs.append(signs[tuple(sorted(join))]*parity(join))
    req(normal_signs==[1]*4,'positive puncture normal orientation')
    core_points=H['oriented_core_point_cycle'];nc=dict(zip(map(tuple,H['N_skeleton']['edge_representatives']),H['N_cocycle']))
    degree=sum(nc[tuple(sorted((a,b)))]*(1 if a<b else -1) for a,b in zip(core_points,core_points[1:]))
    req(degree==1 and core_points[:2]==e,'puncture core does not follow positive longitude')
    req(H['diagram_sign_transport']['ambient_orientation_multiplier_vs_Seifert_master']==1 and H['diagram_sign_transport']['boundary_N_vs_du_dv']==-1,'ambient/normal sign calibration')
    return {'schema':'source-twist-localization-inputs-v1','passed':True,
            'source_sha256':dig(S),'cut_sha256':dig(C),'diagram_chart_sha256':dig(H),
            'mixed_profile_transfer_sha256':dig(MIX),'enlarged_N':MIX['enlarged_N'],
            'fixed_circle_vertices':len(F),'fixed_circle_edges':len(fixed_edges),'fixed_core_cycle':core,
            'fixed_locus_is_full_subcomplex':True,'no_simplex_contains_two_vertices_in_one_deck_orbit':True,
            'fixed_edge_normal_links':edge_rows,'fixed_vertex_links':vertex_rows,
            'regular_neighborhood':'literal epsilon-neighborhood f>=1/2 of the full fixed-circle subcomplex; exact parent partitions bound by cut_orientation receipt',
            'puncture':{'axis':e,'equator_positive':positive,'join_tetrahedron_signs':normal_signs,'core_longitude_degree':degree,
                        'ball':'actual fixed-edge join its four-edge normal circle','boundary':'actual suspension of the normal circle'},
            'theorem_trust_boundary':['PL regular-neighborhood and relative uniqueness theorem',
                'first apply the separately checked simultaneous u/v profile transfer and enlarge N through its saved adaptor',
                'oriented normal D2 bundle over a circle is trivial; prescribe the now-interior local puncture chart',
                'meridian-preserving OUTER GRID pair chart extends, with boundary isotopy supported away from puncture; original cut torus is not simultaneously straightened',
                'lift the framed downstairs pair homeomorphism through the marked double branched cover'],
            'not_claimed':'no global finite D2-bundle atlas, extension cell list or conjugating map is exported by this checker',
            'seconds':time.monotonic()-start}

if __name__=='__main__':
    ans=verify();(ROOT/'central_roll/localization_inputs.json').write_text(json.dumps(ans,indent=2)+'\n')
    print(json.dumps({k:v for k,v in ans.items() if k not in ('fixed_edge_normal_links','fixed_vertex_links')},indent=2),flush=True)

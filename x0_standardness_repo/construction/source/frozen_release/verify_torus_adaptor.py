#!/usr/bin/env python3
"""Replay saved torus-adaptor geometry and coarsening without the producer."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import combinations
import gzip,json,hashlib,sys,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from verify_adaptor_mesh import overlay_check,half_check,same_coarse,grid
from verify_compact3 import replay

def load(p):return json.loads(gzip.decompress(p.read_bytes()))
def tup(x):return tuple(map(tup,x)) if isinstance(x,list) else x
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def req(x,m):
    if not x:raise ValueError(m)

def main():
    start=time.monotonic();cut=load(ROOT/'peripheral/cut_source.json.gz');chart=load(ROOT/'peripheral/diagram_chart.json.gz')
    ac=load(ROOT/'adaptor/adaptor.json.gz');initial=load(ROOT/'adaptor/initial_stages.json.gz');model=load(ROOT/'adaptor/model.json.gz')
    D=cut['deck'];rank=lambda v:min(v,D[v]);xy={v:(Q(x),Q(y)) for v,x,y in chart['upstairs_vertex_coordinates']}
    E={tuple(e):tuple(p) for e,p in zip(chart['torus_upstairs']['edge_representatives'],chart['upstairs_periods'])}
    req(ac['cut_sha256']==digest(cut) and ac['chart_sha256']==digest(chart),'source hash')
    A=[]
    for t in cut['T']:
        a=t[0];co=[]
        for v in t:
            p=(0,0) if v==a else E[a,v]
            co.append(tuple(xy[v][k]+p[k] for k in (0,1)))
        A.append({'ids':tuple(t),'order':tuple(sorted(t,key=rank)),'points':tuple(co)})
    same_coarse(ac['A'],A);G=grid(3,4);same_coarse(ac['G'],G);oc=overlay_check(ac['A'],ac['G'],ac['overlay'],Q(1,2))
    removed=set(map(tuple,cut['puncture']['removed_N_tetrahedra']));bd=list(map(tuple,cut['puncture']['boundary_triangles']))
    req(set(map(tuple,ac['removed_tets']))==removed and list(map(tuple,ac['boundary']))==bd,'source puncture')
    stages={};deck={tup(a):tup(b) for a,b in initial['deck']}
    ai={frozenset(f['ids']):i for i,f in enumerate(A)};gi={frozenset(f['ids']):i for i,f in enumerate(G)}
    da=[ai[frozenset(D[v] for v in f['ids'])] for f in A]
    dg=[gi[frozenset((i,(j+2)%4) for i,j in f['ids'])] for f in G]
    for side in ('N','E'):
        cv=lambda v:(side,'o',v);rv=lambda p:(side,'r',str(p[0]),str(p[1]));gv=lambda q:(side,'g',*q)
        orig=[tuple(cv(v) for v in sorted(t,key=rank)) for t in cut[side] if side!='N' or tuple(t) not in removed]
        req(tup(ac['charts'][side]['original'])==tuple(orig),'source original piece')
        aa=half_check(ac['A'],ac['overlay'],oc,'a',ac['charts'][side]['a'],cv,rv,lambda i:(side,'a',i))
        bb=half_check(ac['G'],ac['overlay'],oc,'b',ac['charts'][side]['b'],gv,rv,lambda i:(side,'b',i))
        ts=orig+aa+bb;stages[side]=ts;req(ts==list(map(tup,initial['stages'][side])),'serialized initial stages')
        for v in {v for t in ts for v in t}:
            kind=v[1]
            if kind=='o':w=(side,'o',D[v[2]])
            elif kind=='r':
                z=Q(v[3])+Q(1,2);z-=z.numerator//z.denominator;w=(side,'r',v[2],str(z))
            elif kind=='g':w=(side,'g',v[2],(v[3]+2)%4)
            elif kind=='a':w=(side,'a',da[v[2]])
            elif kind=='b':w=(side,'b',dg[v[2]])
            else:raise ValueError('unknown source tag')
            req(deck[v]==w,'serialized deck differs from geometry')
        lookup={frozenset(t):t for t in ts};edges={}
        for t in ts:
            req(tuple(deck[v] for v in t)==lookup.get(frozenset(deck[v] for v in t)),'deck/branching failure')
            for a,b in combinations(t,2):
                e=frozenset((a,b));req(e not in edges or edges[e]==(a,b),'inconsistent simplex order');edges[e]=(a,b)
        fc=Counter(frozenset(f) for t in ts for f in combinations(t,3));actual={f for f,c in fc.items() if c==1}
        want={frozenset(gv(v) for v in f['ids']) for f in G}
        if side=='N':want|={frozenset(cv(v) for v in f) for f in bd}
        req(actual==want and all(c in (1,2) for c in fc.values()),'piece boundary/incidence mismatch')
    req(model['source_cut_sha256']==digest(cut) and model['diagram_chart_sha256']==digest(chart) and model['adaptor_sha256']==digest(ac),'model source provenance')
    req(model['L']==3 and model['V']==4 and model['coarsened'],'model grid and scope')
    cert=load(ROOT/'adaptor/coarsening_certificate.json.gz');req(model['coarsening_certificate_sha256']==digest(cert),'coarsening hash')
    compact,_,receipt=replay(cert,stages,deck,ac,chart)
    for side in ('N','E'):req(compact[side]==list(map(tup,model['stages'][side])),'model endpoint mismatch')
    req({tup(a):tup(b) for a,b in model['deck']}==deck,'model deck mismatch')
    out={'schema':'independent-torus-adaptor-replay-v1','passed':True,'source_torus_triangles':len(A),'grid_triangles':len(G),
         'overlay_polygons':oc['polygons'],'refined_triangles':oc['refined_triangles'],
         'half_collar_volume':'1 each, exact','tetrahedra':{s:len(ts) for s,ts in compact.items()},
         'coarsening':receipt,'model_sha256':digest(model),'seconds':time.monotonic()-start}
    (ROOT/'adaptor/independent_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()

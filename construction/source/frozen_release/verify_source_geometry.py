#!/usr/bin/env python3
"""Independent orientation, clipping, and exact torus-chart data checker.

Imports no construction, coarsening, or chart producer. Fractions and integer
chains are used throughout. Prior ball-replacement topology is independently
replayed elsewhere; here the actual replacement boundary chains transport the
geometrically pinned orientation.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict,Counter
from itertools import combinations
import argparse,copy,gzip,hashlib,json,time

ROOT=Path(__file__).resolve().parent

def load(path):
    b=path.read_bytes()
    return json.loads(gzip.decompress(b) if path.suffix=='.gz' else b)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def req(x,msg):
    if not x:raise ValueError(msg)
def sign(x):return 1 if x>0 else -1 if x<0 else 0
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def order_sign(p):return parity([sorted(p).index(v) for v in p])
def mod(x):return x-x.numerator//x.denominator

def determinant(rows):
    a=[list(map(Q,row)) for row in rows];n=len(a);v=Q(1)
    for k in range(n):
        p=next((i for i in range(k,n) if a[i][k]),None)
        if p is None:return Q(0)
        if p!=k:a[p],a[k]=a[k],a[p];v=-v
        pivot=a[k][k];v*=pivot
        for i in range(k+1,n):
            ratio=a[i][k]/pivot
            for j in range(k+1,n):a[i][j]-=ratio*a[k][j]
    return v

def boundary(orientations):
    out=Counter()
    for t,s in orientations.items():
        for k in range(len(t)):out[t[:k]+t[k+1:]]+=s*(-1)**k
    return {t:s for t,s in out.items() if s}

def check_closed_orientation(ori,deck=None):
    req(not boundary(ori),'closed orientation chain has nonzero boundary')
    if deck is not None:
        for t,s in ori.items():
            image=[deck[v] for v in t];u=tuple(sorted(image))
            req(u in ori and ori[u]==s*order_sign(image),'orientation/deck mismatch')

def transport(master,mo,cert,source):
    mf=list(map(tuple,master['closed_cover_facets']))
    req(digest(mf)==mo['facets_sha256']==master['facets_sha256'],'master orientation hash')
    req(len(mf)==len(mo['tetrahedron_orientations']),'master orientation length')
    ori=dict(zip(mf,mo['tetrahedron_orientations']))
    req(set(ori.values())<={-1,1},'master orientation signs')
    check_closed_orientation(ori,master['deck'])
    tags=cert['tags'];req(len(tags)==len(set(tags)),'repeated certificate tags')
    initial={tuple(sorted(tags[v] for v in t)) for t in cert['initial']}
    req(initial==set(ori),'certificate initial source mismatch')
    req(digest(cert['initial'])==cert['initial_sha256'],'certificate initial hash')
    individual=0;inserted=0
    def step(u,v):
        nonlocal individual,inserted
        old={t:s for t,s in ori.items() if u in t}
        req(old and any(v in t for t in old),'replacement edge missing')
        new={}
        for t,s in old.items():
            f=tuple(w for w in t if w!=u)
            if v in f:continue
            q=tuple(sorted((v,*f)))
            req(q not in new,'duplicate new tetrahedron')
            req(q not in ori or q in old,'new tetrahedron meets outside')
            new[q]=s*(-1)**(t.index(u)+q.index(v))
        req(new and boundary(old)==boundary(new),'oriented replacement boundary differs')
        for t in old:del ori[t]
        ori.update(new);individual+=1;inserted+=len(new)
        return len(old),len(new)
    for event in cert['steps']:
        req(len(ori)==event['before'],'recorded replacement source size')
        u,v,U,V=[tags[event[k]] for k in ('u','v','deck_u','deck_v')]
        req(master['deck'][u]==U and master['deck'][v]==V,'recorded replacement deck pair')
        a,b=step(u,v);c,d=step(U,V)
        req([a,c]==event['deleted_per_side'] and [b,d]==event['inserted_per_side'],'recorded replacement sizes')
        req(len(ori)==event['after'],'recorded replacement target size')
    final={tuple(sorted(tags[v] for v in t)) for t in cert['endpoint']}
    req(set(ori)==final and digest(cert['endpoint'])==cert['endpoint_sha256'],'certificate endpoint mismatch')
    parent=source['parent_vertices'];endpoint={}
    for t in map(tuple,source['closed_cover_facets']):
        lifted=[parent[v] for v in t];old=tuple(sorted(lifted))
        req(old in ori,'renumbered endpoint missing')
        endpoint[t]=ori[old]*order_sign(lifted)
    req(len(endpoint)==len(ori),'renumbered endpoint multiplicity')
    check_closed_orientation(endpoint,source['deck'])
    result={'schema':'source-orientation-transport-v1','passed':True,'master_facets_sha256':mo['facets_sha256'],
            'source_facets_sha256':digest(source['closed_cover_facets']),
            'tetrahedron_orientations':[endpoint[tuple(t)] for t in source['closed_cover_facets']],
            'coarsening_deck_pair_steps':len(cert['steps']),'individual_replacements':individual,
            'new_tetrahedron_orientation_assignments':inserted,
            'orientation_source':mo['orientation_source'],
            'method':'exact equality of oriented replacement boundary chains; retained orientation on unchanged tetrahedra'}
    return endpoint,result

def clip_orient(source,ori,cut):
    facets=list(map(tuple,source['closed_cover_facets']));pts=cut['points'];F={i for i,v in enumerate(source['deck']) if i==v}
    req(set(cut['fixed_source_vertices'])==F,'cut fixed set')
    req(cut['source_sha256']==digest(source),'cut source hash')
    req(len(pts)==len({tuple(p) for p in pts}),'repeated cut point')
    pi={tuple(p):i for i,p in enumerate(pts)}
    for i,p in enumerate(pts):
        if p[0]=='v':image=('v',source['deck'][p[1]])
        elif p[0]=='m':
            req(len(p)==3 and (p[1] in F)!=(p[2] in F),'nontransverse midpoint')
            image=('m',*sorted(source['deck'][v] for v in p[1:]))
        else:raise ValueError('unknown cut point')
        req(pi.get(image)==cut['deck'][i],'cut point/deck mismatch')
    parent_children=defaultdict(list);out={'N':{},'E':{}};volumes=defaultdict(Q)
    def bary(i,t):
        p=pts[i];row=[Q(0)]*4
        if p[0]=='v':
            req(p[1] in t,'cut vertex outside carrier');row[t.index(p[1])]=1
        else:
            req(set(p[1:])<=set(t),'cut midpoint outside carrier')
            for v in p[1:]:row[t.index(v)]=Q(1,2)
        return row
    carriers_seen=set()
    for carrier in cut['top_carriers']:
        side=carrier['side'];i=carrier['source_tetrahedron_index'];req(side in out and 0<=i<len(facets),'cut carrier index')
        req((side,i) not in carriers_seen,'duplicate top carrier');carriers_seen.add((side,i));t=facets[i]
        for child in carrier['children']:
            q=tuple(child);req(q==tuple(sorted(q)) and len(set(q))==4,'cut child ordering')
            req(q not in out[side],'duplicate cut child');B=[bary(v,t) for v in q];d=determinant(B)
            req(d!=0,'degenerate cut child')
            levels=[sum(r[k] for k,v in enumerate(t) if v in F) for r in B]
            req(all(x>=Q(1,2) for x in levels) if side=='N' else all(x<=Q(1,2) for x in levels),'cut child on wrong side')
            out[side][q]=ori[t]*sign(d);volumes[i]+=abs(d);parent_children[i].append((q,sign(d),B))
    req(set(parent_children)==set(range(len(facets))),'cut omitted parent tetrahedron')
    face_partitions=0
    for i,t in enumerate(facets):
        req(volumes[i]==1,'cut child volume sum differs from parent')
        local={q:s for q,s,B in parent_children[i]};bd=boundary(local);faceareas=defaultdict(Q)
        for f,s in bd.items():
            B=[bary(v,t) for v in f];zero=[k for k in range(4) if all(r[k]==0 for r in B)]
            req(len(zero)==1,'cut parent has spurious internal boundary')
            k=zero[0];d=determinant([[r[j] for j in range(4) if j!=k] for r in B])
            req(d and s*sign(d)==(-1)**k,'cut parent boundary orientation')
            faceareas[k]+=abs(d)
        req(all(faceareas[k]==1 for k in range(4)),'cut parent face area partition');face_partitions+=4
    for side in out:req(set(out[side])==set(map(tuple,cut[side])),'cut output differs from carried children')
    nb,eb=boundary(out['N']),boundary(out['E']);T=set(map(tuple,cut['T']))
    req(set(nb)==set(eb)==T,'tube/exterior boundary differs from torus')
    req(all(nb[f]==-eb[f] and abs(nb[f])==1 for f in T),'tube/exterior boundary orientation')
    check_closed_orientation({**out['N'],**out['E']},cut['deck'])
    rec={'schema':'exact-clipped-source-orientations-v1','passed':True,'source_sha256':digest(source),'cut_sha256':digest(cut),
         'N_tetrahedron_orientations':[out['N'][tuple(t)] for t in cut['N']],
         'E_tetrahedron_orientations':[out['E'][tuple(t)] for t in cut['E']],
         'T_boundary_N_orientations':[nb[tuple(t)] for t in cut['T']],
         'T_boundary_E_orientations':[eb[tuple(t)] for t in cut['T']],
         'parent_tetrahedra_checked':len(facets),'child_tetrahedra_checked':sum(map(len,out.values())),
         'parent_face_partitions_checked':face_partitions,
         'method':'exact child barycentric determinants, unit parent volumes, oriented parent-face unit-area partitions, and global boundary chains'}
    return nb,rec

def expected_skeleton(facets,D,quotient):
    rows=[]
    for k in (1,2,3):
        cells={tuple(sorted(f)) for t in facets for f in combinations(t,k)}
        if quotient:cells={min(f,tuple(sorted(D[v] for v in f))) for f in cells}
        rows.append(sorted(cells))
    ei={e:i+1 for i,e in enumerate(rows[1])}
    def edge(a,b):
        r=tuple(sorted((a,b)))
        if quotient:r=min(r,tuple(sorted((D[a],D[b]))))
        sg=1 if (a,b)==r or (quotient and (D[a],D[b])==r) else -1
        return sg*ei[r]
    return {'vertices':[x[0] for x in rows[0]],'edge_representatives':[list(x) for x in rows[1]],
            'edge_endpoints':[[min(a,D[a]),min(b,D[b])] if quotient else [a,b] for a,b in rows[1]],
            'triangle_representatives':[list(x) for x in rows[2]],
            'triangle_boundaries':[[edge(a,b),edge(b,c),edge(c,a)] for a,b,c in rows[2]]}

def chart_check(cut,chart,tori):
    D=cut['deck'];req(chart['cut_sha256']==digest(cut),'chart cut hash')
    for name,side,q in [('torus_downstairs','T',True),('torus_upstairs','T',False),('N_skeleton','N',True),('E_skeleton','E',True)]:
        req(chart[name]==expected_skeleton(cut[side],D,q),'chart skeleton differs from actual cut: '+name)
    for side in ('N','E'):
        sk=chart[side+'_skeleton'];co=chart[side+'_cocycle']
        req(len(co)==len(sk['edge_endpoints']) and all(isinstance(x,int) for x in co),'invalid '+side+' cocycle')
        for tr in sk['triangle_boundaries']:req(sum((1 if e>0 else -1)*co[abs(e)-1] for e in tr)==0,'nonclosed '+side+' cocycle')
    skd=chart['torus_downstairs'];sku=chart['torus_upstairs']
    coN=dict(zip(map(tuple,chart['N_skeleton']['edge_representatives']),chart['N_cocycle']))
    coE=dict(zip(map(tuple,chart['E_skeleton']['edge_representatives']),chart['E_cocycle']))
    req(chart['downstairs_periods']==[[coN[tuple(e)],coE[tuple(e)]] for e in skd['edge_representatives']],'periods not restricted tube/exterior cocycles')
    summary={};edge_displacements={};coords={}
    for name,sk,quotient in [('downstairs',skd,True),('upstairs',sku,False)]:
        periods=chart[name+'_periods'];req(len(periods)==len(sk['edge_endpoints']),'period count')
        req(all(len(p)==2 and all(isinstance(x,int) for x in p) for p in periods),'noninteger lattice period')
        xy={v:(Q(x),Q(y)) for v,x,y in chart[name+'_vertex_coordinates']}
        req(set(xy)==set(sk['vertices']),'vertex coordinates incomplete');coords[name]=xy
        disp=[]
        for (a,b),p in zip(sk['edge_endpoints'],periods):disp.append(tuple(xy[b][k]+p[k]-xy[a][k] for k in range(2)))
        edge_displacements[name]=disp;dets=[];signed=Q(0)
        for tr,bd in zip(sk['triangle_representatives'],sk['triangle_boundaries']):
            req(all(sum((1 if e>0 else -1)*periods[abs(e)-1][k] for e in bd)==0 for k in range(2)),'nonclosed torus periods')
            ab=tuple((1 if bd[0]>0 else -1)*z for z in disp[abs(bd[0])-1])
            ac=tuple(-(1 if bd[2]>0 else -1)*z for z in disp[abs(bd[2])-1])
            d=ab[0]*ac[1]-ab[1]*ac[0];req(d!=0,'degenerate '+name+' chart triangle')
            s=tori[tuple(sorted(tr))]*order_sign(tr);d*=s;dets.append(d);signed+=d/2
        req(len({sign(d) for d in dets})==1,'folded '+name+' chart')
        req(abs(signed)==1,'chart degree/area differs from one')
        summary[name]={'triangles':len(dets),'orientation_relative_to_boundary_N':sign(dets[0]),
            'signed_total_area':str(signed),'absolute_total_area':'1','all_determinants_nonzero':True,
            'minimum_absolute_double_area':str(min(map(abs,dets))),
            'maximum_absolute_double_area':str(max(map(abs,dets)))}
    dxy,uxy=coords['downstairs'],coords['upstairs']
    for v,xy in uxy.items():
        w=D[v];req(w in uxy,'upstairs deck vertex missing')
        req(tuple(mod(uxy[w][k]-xy[k]) for k in range(2))==(Q(0),Q(1,2)),'deck is not the stated half shift')
        r=min(v,D[v]);req(mod(dxy[r][0]-xy[0])==0 and mod(dxy[r][1]-2*xy[1])==0,'quotient coordinate factor two')
    ui={tuple(e):i for i,e in enumerate(sku['edge_representatives'])};di={tuple(e):i for i,e in enumerate(skd['edge_representatives'])}
    for j,(a,b) in enumerate(sku['edge_representatives']):
        pair=(D[a],D[b]);e=tuple(sorted(pair));sg=1 if pair==e else -1
        req(tuple(sg*z for z in edge_displacements['upstairs'][ui[e]])==edge_displacements['upstairs'][j],'deck edge derivative differs from identity')
        r=min((a,b),e);sg=1 if (a,b)==r or pair==r else -1;k=di[r]
        u=edge_displacements['upstairs'][j];d=tuple(sg*z for z in edge_displacements['downstairs'][k])
        req(d==(u[0],2*u[1]),'cover projection edge derivative')
    for name,want in [('longitude',[1,0]),('meridian',[0,1])]:
        rec=chart['peripheral_loops'][name];word=rec['torus_edge_loop'];req(word,'empty peripheral loop')
        vals=[0,0];first=None;end=None
        for e in word:
            a,b=skd['edge_endpoints'][abs(e)-1]
            if e<0:a,b=b,a
            if first is None:first=a
            if end is not None:req(end==a,'disconnected peripheral edge path')
            end=b
            for k in range(2):vals[k]+=(1 if e>0 else -1)*chart['downstairs_periods'][abs(e)-1][k]
        req(end==first,'open peripheral loop');req(vals==rec['periods']==want,'wrong peripheral loop periods')
    return {'schema':'independent-exact-torus-chart-check-v1','passed':True,'charts':summary,
            'deck_vertices_checked':len(uxy),'deck_edges_checked':len(sku['edge_representatives']),
            'deck_translation':['0','1/2'],'cover_projection':[[1,0],[0,2]],
            'source_tube_and_exterior_cocycles_closed':True,'peripheral_basis_periods_verified':True}

def mutations(master,mo,cert,source,cut,chart,ori,tori):
    done=[]
    def reject(label,fn):
        try:fn()
        except (ValueError,KeyError,IndexError):done.append(label)
        else:raise ValueError('mutation accepted: '+label)
    changed=dict(ori);k=next(iter(changed));changed[k]*=-1
    reject('one transported tetrahedron orientation flipped',lambda:check_closed_orientation(changed,source['deck']))
    bad=copy.deepcopy(cut);bad['top_carriers'][0]['source_tetrahedron_index']=(bad['top_carriers'][0]['source_tetrahedron_index']+1)%len(source['facets'])
    reject('cut child assigned wrong original tetrahedron',lambda:clip_orient(source,ori,bad))
    bad=copy.deepcopy(chart);bad['upstairs_periods'][0][0]+=1
    reject('one upstairs lattice period altered',lambda:chart_check(cut,bad,tori))
    bad=copy.deepcopy(chart);bad['upstairs_vertex_coordinates'][0][2]=str(Q(bad['upstairs_vertex_coordinates'][0][2])+Q(1,3))
    reject('one upstairs coordinate breaks covering/deck map',lambda:chart_check(cut,bad,tori))
    bad=copy.deepcopy(chart);bad['peripheral_loops']['longitude']['periods']=[0,1]
    reject('longitude and meridian periods confused',lambda:chart_check(cut,bad,tori))
    return done

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--mutations',action='store_true');args=ap.parse_args();r=args.root;t=time.monotonic()
    paths={'master':r/'coarsened/master_source.json.gz','mo':r/'seifert/data/master_orientation.json','cert':r/'coarsened/coarsening_certificate.json.gz',
           'source':r/'coarsened/source.json.gz','cut':r/'peripheral/cut_source.json.gz','chart':r/'peripheral/peripheral_chart.json.gz'}
    a={k:load(p) for k,p in paths.items()}
    ori,o=transport(a['master'],a['mo'],a['cert'],a['source']);o['input_file_sha256']={k:sha(paths[k]) for k in ('master','mo','cert','source')}
    (r/'coarsened/source_orientation.json').write_text(json.dumps(o,indent=2)+'\n')
    print('orientation transport PASS',len(ori),flush=True)
    tori,c=clip_orient(a['source'],ori,a['cut']);c['source_orientation_file_sha256']=sha(r/'coarsened/source_orientation.json')
    (r/'peripheral/cut_orientation.json').write_text(json.dumps(c,indent=2)+'\n');print('cut orientation PASS',c['child_tetrahedra_checked'],flush=True)
    ch=chart_check(a['cut'],a['chart'],tori);ch['input_file_sha256']={k:sha(paths[k]) for k in ('cut','chart')}
    if args.mutations:ch['mutation_checks_rejected']=mutations(a['master'],a['mo'],a['cert'],a['source'],a['cut'],a['chart'],ori,tori)
    ch['seconds']=time.monotonic()-t
    (r/'peripheral/chart_independent_verification.json').write_text(json.dumps(ch,indent=2)+'\n')
    print(json.dumps(ch,indent=2),flush=True)

if __name__=='__main__':main()

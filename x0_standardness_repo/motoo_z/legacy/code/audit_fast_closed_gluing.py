#!/usr/bin/env python3
"""Independent raw-facet audit of the saved Figure2 closed gluing.

No construction helper is imported. No4D skeleton, homology, or vertex
recognition is requested. This verifies the exact gluing combinatorics;
manifold admission and standardness are separate gates.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import regina

def require(ok,msg):
    if not ok:raise ValueError(msg)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def boundary_pairings(t):
    places=[];lookup={}
    for s in t.simplices():
        for f in range(5):
            if s.adjacentSimplex(f) is None:
                lookup[s.index(),f]=len(places)
                places.append((s.index(),f,tuple(v for v in range(5) if v!=f)))
    pairings={};maxwalk=0
    for i,(index,outer,vertices) in enumerate(places):
        for omitted in range(4):
            current=t.simplex(index);a=outer;b=vertices[omitted]
            carry={j:vertices[j] for j in range(4) if j!=omitted};seen=set()
            while current.adjacentSimplex(b) is not None:
                key=(current.index(),a,b)
                require(key not in seen,'codimension-two walk is cyclic without a boundary end');seen.add(key)
                p=current.adjacentGluing(b);other=current.adjacentSimplex(b)
                carry={j:p[v] for j,v in carry.items()};a,b=p[b],p[a];current=other
            maxwalk=max(maxwalk,len(seen))
            j=lookup[current.index(),b];target_vertices=places[j][2]
            target_omit=target_vertices.index(a)
            permutation=tuple(target_omit if k==omitted else target_vertices.index(carry[k]) for k in range(4))
            require(sorted(permutation)==list(range(4)),'bad induced boundary permutation')
            pairings[i,omitted]=(j,permutation)
    for (i,f),(j,p) in pairings.items():
        k,q=pairings[j,p[f]]
        require(k==i and all(q[p[v]]==v for v in range(4)),'boundary pairing is not reciprocal')
    return places,pairings,maxwalk

def audit(folder,output):
    folder,output=Path(folder),Path(output);now=time.monotonic()
    receipt=json.loads((folder/'RAW_CLOSED_GLUING.json').read_text())
    first_path=Path(receipt['first_source']);second_path=Path(receipt['second_source']);closed_path=folder/'closed_candidate.rga'
    require(sha(first_path)==receipt['first_sha256'],'first input SHA differs')
    require(sha(second_path)==receipt['second_sha256'],'second input SHA differs')
    require(sha(closed_path)==receipt['closed_rga_sha256'],'closed input SHA differs')
    first=regina.Triangulation4(regina.open(str(first_path)))
    second=regina.Triangulation4(regina.open(str(second_path)))
    closed=regina.Triangulation4(regina.open(str(closed_path)))
    offset=first.size();require(closed.size()==first.size()+second.size(),'incorrect total simplex count')
    maps=json.loads((folder/'boundary_gluing.json').read_text())
    by_source={(m['source_pentachoron'],m['source_facet']):m for m in maps}
    by_target={(m['target_pentachoron'],m['target_facet']):m for m in maps}
    require(len(by_source)==len(by_target)==len(maps),'repeated source or target attaching facet')
    counts={'preserved_internal_directed_facets':0,'checked_attaching_directed_facets':0}
    for piece,shift in [(first,0),(second,offset)]:
        for s in piece.simplices():
            cs=closed.simplex(shift+s.index())
            for f in range(5):
                old=s.adjacentSimplex(f);new=cs.adjacentSimplex(f)
                require(new is not None,'an unpaired facet remains in the closed file')
                actual=tuple(cs.adjacentGluing(f)[v] for v in range(5))
                if old is not None:
                    require(new.index()==shift+old.index(),'interior simplex adjacency changed')
                    require(actual==tuple(s.adjacentGluing(f)[v] for v in range(5)),'interior vertex permutation changed')
                    counts['preserved_internal_directed_facets']+=1
                elif shift==0:
                    m=by_source[s.index(),f]
                    require(new.index()==offset+m['target_pentachoron'],'wrong attaching target')
                    require(actual==tuple(m['permutation']),'wrong attaching permutation')
                    counts['checked_attaching_directed_facets']+=1
                else:
                    m=by_target[s.index(),f];p=m['permutation']
                    inv=tuple(p.index(v) for v in range(5))
                    require(new.index()==m['source_pentachoron'] and actual==inv,'wrong reciprocal attaching map')
                    counts['checked_attaching_directed_facets']+=1
    require(sum(counts.values())==5*closed.size(),'not all closed facets were checked')
    ap,aa,aw=boundary_pairings(first);bp,bb,bw=boundary_pairings(second)
    require(len(ap)==len(bp)==len(maps),'boundary tetrahedron count differs')
    simplex_map={};vertex_map={}
    for m in maps:
        i=m['source_boundary_tetrahedron'];j=m['target_boundary_tetrahedron']
        ai,af,av=ap[i];bi,bf,bv=bp[j]
        require((ai,af)==(m['source_pentachoron'],m['source_facet']),'source boundary indexing differs')
        require((bi,bf)==(m['target_pentachoron'],m['target_facet']),'target boundary indexing differs')
        p=m['permutation'];require(p[af]==bf,'attaching map misses target tetrahedron')
        simplex_map[i]=j;vertex_map[i]=tuple(bv.index(p[v]) for v in av)
    require(set(simplex_map)==set(range(len(ap))) and set(simplex_map.values())==set(range(len(bp))),
            'boundary simplex map is not a bijection')
    triangles=0
    for (i,f),(other,p) in aa.items():
        j=simplex_map[i];image=vertex_map[i];other_j,q=bb[j,image[f]]
        require(other_j==simplex_map[other],'boundary map does not preserve triangle adjacencies')
        other_image=vertex_map[other]
        require(all(other_image[p[v]]==q[image[v]] for v in range(4)),
                'boundary map does not commute with a triangle vertex gluing')
        triangles+=1
    report={'status':'PASS_INDEPENDENT_RAW_CLOSED_GLUING','first_sha256':sha(first_path),
            'second_sha256':sha(second_path),'closed_rga_sha256':sha(closed_path),
            'closed_tight_sha256':hashlib.sha256(closed.tightEncoding().encode()).hexdigest(),
            'first_pentachora':first.size(),'second_pentachora':second.size(),'closed_pentachora':closed.size(),
            'boundary_tetrahedra':len(ap),'boundary_directed_triangles_checked':triangles,
            'all_closed_facet_pairings_checked':True,'every_original_internal_pairing_preserved':True,
            'attaching_map_is_explicit_combinatorial_boundary_isomorphism':True,
            'maximum_boundary_triangle_walk':[aw,bw],**counts,
            'four_dimensional_skeleton_used':False,'manifold_recognition_performed':False,
            'standardness_asserted':False,'seconds':time.monotonic()-now,'checker_sha256':sha(__file__)}
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report),flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();audit(a.folder,a.out)

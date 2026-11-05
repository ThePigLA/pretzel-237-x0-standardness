#!/usr/bin/env python3
"""Independently reconstruct and replay the 24-to-9 standard upper cap.

The constructor, search driver and its move helper are not imported.
The source is a boundary sum of two explicit Delta3-times-circle products.
All fifteen moves in this certificate are guarded edge collapses.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import regina


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def tight(t):return hashlib.sha256(t.tightEncoding().encode()).hexdigest()


def construct_upper24(proof):
    # The three edges are oriented 0->1, 1->2, 2->0.  Over each edge,
    # use the four ordered staircase cells of Delta3 x interval.
    cells=[[(layer,j) for j in range(split+1)]+
           [((layer+1)%3,j) for j in range(split,4)]
           for layer in range(3) for split in range(4)]
    require([[list(v) for v in c] for c in cells]==proof['product_simplex_coordinates'],
            'recorded product coordinates differ from the explicit staircase')
    product=regina.Triangulation4();product.newSimplices(len(cells))
    facets=defaultdict(list)
    for i,cell in enumerate(cells):
        require(len(set(cell))==5,'degenerate product cell')
        for f in range(5):
            facets[frozenset(cell[:f]+cell[f+1:])].append((i,f))
    for key,occurrences in facets.items():
        require(len(occurrences) in (1,2),'product facet has invalid incidence')
        if len(occurrences)==1:continue
        (i,f),(j,g)=occurrences
        permutation=[g if v==f else cells[j].index(cells[i][v]) for v in range(5)]
        product.simplex(i).join(f,product.simplex(j),regina.Perm5(*permutation))
    require(product.size()==12 and product.isValid() and not product.isIdeal()
            and product.isConnected() and product.countBoundaryComponents()==1,
            'independent product construction is not a finite connected manifold')
    connecting=proof['embedded_connecting_tetrahedron']
    i,f=connecting['pentachoron'],connecting['facet']
    require(type(i) is int and 0<=i<12 and type(f) is int and 0<=f<5,
            'invalid connecting facet')
    require(product.simplex(i).adjacentSimplex(f) is None,'connecting facet is not on the boundary')
    require(len({product.simplex(i).vertex(v).index() for v in range(5) if v!=f})==4,
            'connecting boundary tetrahedron is not embedded')
    upper=regina.Triangulation4(product);upper.insertTriangulation(product)
    upper.simplex(i).join(f,upper.simplex(i+12),regina.Perm5())
    require(upper.size()==24 and upper.isValid() and not upper.isIdeal()
            and upper.isConnected() and upper.isOrientable(),
            'independent boundary sum is not a finite orientable manifold')
    return upper


def verify(directory,target_path,output):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina 7.4.1')
    started=time.monotonic();directory,target_path,output=map(Path,(directory,target_path,output))
    files={'construction':directory/'SOURCE_CONSTRUCTION.json',
           'result':directory/'RESULT.json','moves':directory/'moves.jsonl',
           'source':directory/'constructed_upper24.rga',
           'endpoint':directory/'guarded_upper9.rga','target':target_path}
    hashes={key:sha(path) for key,path in files.items()};checker_hash=sha(__file__)
    proof=json.loads(files['construction'].read_text());claim=json.loads(files['result'].read_text())
    require(hashes['source']==proof['source_sha256']==claim['source_sha256'],'source file binding differs')
    require(hashes['moves']==claim['moves_sha256'],'move log binding differs')
    require(hashes['target']==proof['target_sha256']==claim['target_file_sha256'],'saved target binding differs')
    require(hashes['endpoint']==claim['endpoint_sha256'],'endpoint file binding differs')
    t=construct_upper24(proof)
    source=regina.Triangulation4(regina.open(str(files['source'])))
    target=regina.Triangulation4(regina.open(str(files['target'])))
    endpoint=regina.Triangulation4(regina.open(str(files['endpoint'])))
    require(t==source and tight(t)==proof['source_tight_sha256']==claim['source_tight_sha256'],
            'independent geometric source differs from the recorded labelled source')
    require(target.size()==9 and target.isoSig()=='jafbGbafqajaadiLaNa',
            'unexpected saved upper target')
    require(all(sha(path)==hashes[key] for key,path in files.items()),'input changed during loading')
    steps=0
    for line in files['moves'].read_text().splitlines():
        row=json.loads(line)
        require(row['index']==steps and row['before_sha256']==tight(t),'discontinuous move chain')
        operation=row['operation']
        require(isinstance(operation,list) and len(operation)==3 and operation[:2]==['C',1],
                'this independent certificate accepts only the recorded edge collapses')
        index=operation[2]
        require(type(index) is int and 0<=index<t.countEdges(),'invalid collapse edge')
        old_size,old_vertices=t.size(),t.countVertices()
        edge=t.edge(index)
        require(t.hasCollapseEdge(edge),'edge-collapse guard failed')
        require(t.collapseEdge(edge),'guarded edge collapse execution failed')
        require(t.size()<old_size and t.countVertices()==old_vertices-1,'wrong collapse size change')
        require(t.size()==row['pentachora'] and tight(t)==row['after_sha256'],
                'replayed labelled state differs')
        steps+=1
    require(steps==15==claim['moves'],'wrong number of guarded collapses')
    require(t==endpoint and t==target and tight(t)==claim['endpoint_tight_sha256'],
            'replay does not reach the exact saved labelled target')
    iso=claim['isomorphism_to_saved_target']
    require(iso['simplices']==list(range(9)) and iso['permutations']==[list(range(5)) for _ in range(9)],
            'claimed endpoint identification is not the verified identity')
    require(all(sha(path)==hashes[key] for key,path in files.items()),'input changed during verification')
    require(sha(__file__)==checker_hash,'checker changed during verification')
    receipt={'status':'PASS_INDEPENDENT_GUARDED_STANDARD_UPPER_CAP',
             'source_pentachora':24,'target_pentachora':9,'moves':steps,
             'guarded_edge_collapses':steps,'source_tight_sha256':tight(source),
             'target_tight_sha256':tight(t),'target_isosig':t.isoSig(),
             'independent_coordinate_product_reconstruction':True,
             'embedded_boundary_sum_tetrahedron_checked':True,
             'exact_labelled_saved_target_reached':True,
             'all_move_guards_checked':True,'constructor_and_search_code_imported':False,
             'native_4d_simplify_used':False,
             'input_sha256':hashes,'input_hashes_stable_before_loading_after_loading_and_at_end':True,
             'regina_distribution':'7.4.1','checker_sha256':checker_hash,
             'trust_boundary':'Explicit staircase product and embedded boundary connected sum; Regina 7.4.1 guarded edge-collapse implementation preserves the PL manifold.',
             'conclusion':'The saved nine-pentachoron reference is PL homeomorphic to the boundary connected sum of two copies of Delta3 x S1.',
             'standardness_of_any_closed_4_manifold_asserted':False,
             'seconds':time.monotonic()-started}
    output.parent.mkdir(parents=True,exist_ok=True);tmp=output.with_suffix(output.suffix+'.tmp')
    tmp.write_text(json.dumps(receipt,indent=2)+'\n');tmp.replace(output)
    print(json.dumps(receipt),flush=True);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',required=True,type=Path);p.add_argument('--target',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);a=p.parse_args();verify(a.certificate,a.target,a.out)

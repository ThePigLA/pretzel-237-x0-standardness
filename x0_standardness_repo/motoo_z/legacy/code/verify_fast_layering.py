#!/usr/bin/env python3
"""Check every inverse4D shelling of an incrementally layered piece."""
from __future__ import annotations
import argparse,gzip,hashlib,json,time
from pathlib import Path
import regina

def require(ok,msg):
    if not ok: raise ValueError(msg)

def verify(source,folder,out):
    source,folder,out=Path(source),Path(folder),Path(out)
    start=time.monotonic()
    original=regina.Triangulation4(regina.open(str(source)))
    result=regina.Triangulation4(regina.open(str(folder/'layered_piece.rga')))
    with gzip.open(folder/'layering_events.jsonl.gz','rt') as f:
        events=[json.loads(line) for line in f]
    require(result.size()==original.size()+len(events),'wrong layered size')
    print(json.dumps({'stage':'begin_native_inverse_shelling_guards','source_pentachora':original.size(),
                      'layered_pentachora':result.size(),'moves':len(events)}),flush=True)
    for count,e in enumerate(reversed(events),1):
        i=e['added_pentachoron']
        require(i==result.size()-1,'layer append order differs')
        pent=result.simplex(i)
        for g in e['gluings']:
            adjacent=pent.adjacentSimplex(g['new_facet'])
            require(adjacent is not None and adjacent.index()==g['old_pentachoron'],'wrong attached neighbour')
            p=pent.adjacentGluing(g['new_facet'])
            require([p[j] for j in range(5)]==g['permutation'],'wrong attaching permutation')
        require(result.shellBoundary(pent),f'guarded inverse4D shelling failed at reversed step{count}')
        if count%100==0:
            print(json.dumps({'inverse_shellings_checked':count,'seconds':time.monotonic()-start}),flush=True)
    require(result==original,'guarded inverse shellings do not recover exact labelled source')
    report={'status':'PASS_ALL_INVERSE_4D_SHELLING_GUARDS','moves':len(events),
            'source_pentachora':original.size(),'source_recovered_exactly':True,
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'layered_sha256':hashlib.sha256((folder/'layered_piece.rga').read_bytes()).hexdigest(),
            'seconds':time.monotonic()-start}
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report),flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True)
    p.add_argument('--folder',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();verify(a.source,a.folder,a.out)

#!/usr/bin/env python3
"""Bounded guarded search from the explicit upper handlebody to a fixed cap.

This never calls four-dimensional simplify().  The construction is the
boundary sum of two staircase Delta3 x a three-edge circle products.
Only a successful finite path is exported; failure is inconclusive.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import random
import time
import regina
from complete_upper_handles import reference_cap
from reduce_finite_checked import legal,perform


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def tight(t):return hashlib.sha256(t.tightEncoding().encode()).hexdigest()
def save(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def search(target_path,out,seconds=120):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina 7.4.1')
    target_path,out=Path(target_path),Path(out);out.mkdir(parents=True,exist_ok=False)
    source,construction=reference_cap(simplify=False)
    target_sha=sha(target_path);target=regina.Triangulation4(regina.open(str(target_path)))
    require(target.size()==9 and target.isoSig()=='jafbGbafqajaadiLaNa','unexpected target')
    source_path=out/'constructed_upper24.rga'
    require(regina.PacketOfTriangulation4(source).save(str(source_path)),'could not save source')
    save(out/'SOURCE_CONSTRUCTION.json',{
        **construction,'source_sha256':sha(source_path),'source_tight_sha256':tight(source),
        'constructor_sha256':sha(Path(__file__).with_name('complete_upper_handles.py')),
        'target_sha256':target_sha,'target_tight_sha256':tight(target),
        'target_isosig':target.isoSig(),'regina_distribution':'7.4.1',
        'native_4d_simplify_used':False})
    started=time.monotonic();attempt=0;endpoints=Counter();best=source.size()
    priorities=[
        [('C',1),('S',4),('Z',2),('Z',1),('P',1),('P',0),('Z',0)],
        [('S',4),('C',1),('Z',2),('Z',1),('P',1),('P',0),('Z',0)],
        [('P',0),('P',1),('Z',2),('Z',1),('Z',0),('C',1),('S',4)],
    ]
    while time.monotonic()-started<seconds:
        rng=random.Random(attempt);tri=regina.Triangulation4(source);events=[]
        neutral=0;priority=priorities[attempt%len(priorities)]
        while time.monotonic()-started<seconds and len(events)<1000:
            iso=tri.isIsomorphicTo(target) if tri.size()==target.size() else None
            if iso is not None:
                require(sha(target_path)==target_sha,'target changed during search')
                endpoint=out/'guarded_upper9.rga'
                require(regina.PacketOfTriangulation4(tri).save(str(endpoint)),'could not save endpoint')
                log=out/'moves.jsonl';tmp=out/'moves.tmp'
                tmp.write_text(''.join(json.dumps(row,separators=(',',':'))+'\n' for row in events));tmp.replace(log)
                report={'status':'PASS_GUARDED_PATH_TO_SAVED_UPPER_CAP','attempt':attempt,
                    'moves':len(events),'source_pentachora':source.size(),'target_pentachora':tri.size(),
                    'source_sha256':sha(source_path),'source_tight_sha256':tight(source),
                    'target_file_sha256':target_sha,'endpoint_sha256':sha(endpoint),
                    'endpoint_tight_sha256':tight(tri),'endpoint_isosig':tri.isoSig(),
                    'moves_sha256':sha(log),'seconds':time.monotonic()-started,
                    'isomorphism_to_saved_target':{'simplices':[iso.simpImage(i) for i in range(tri.size())],
                         'permutations':[[iso.facetPerm(i)[j] for j in range(5)] for i in range(tri.size())]},
                    'move_counts':dict(Counter(str(row['operation'][:2]) for row in events)),
                    'all_moves_guarded':True,'native_4d_simplify_used':False,
                    'standardness_of_any_closed_4_manifold_asserted':False}
                save(out/'RESULT.json',report);print(json.dumps(report),flush=True);return report
            choices=[]
            for kind,dimension in priority:
                choices=legal(tri,kind,dimension)
                if choices:break
            reducing=bool(choices)
            if not choices and neutral<40:
                choices=legal(tri,'P',2)
            if not choices:break
            operation=choices[0] if attempt<3 else rng.choice(choices)
            before=tight(tri);perform(tri,operation)
            events.append({'index':len(events),'operation':operation,'before_sha256':before,
                           'after_sha256':tight(tri),'pentachora':tri.size()})
            neutral=0 if reducing else neutral+1
            if tri.size()<best:
                best=tri.size();print(json.dumps({'attempt':attempt,'best_pentachora':best,
                                    'moves':len(events),'seconds':time.monotonic()-started}),flush=True)
            if tri.size()<target.size():break
        endpoints[tri.isoSig()]+=1;attempt+=1
    report={'status':'INCONCLUSIVE_BOUNDED_SEARCH','attempts':attempt,'best_pentachora':best,
            'seconds':time.monotonic()-started,'endpoint_counts':dict(endpoints),
            'native_4d_simplify_used':False,'standardness_asserted':False}
    save(out/'RESULT.json',report);print(json.dumps(report),flush=True);return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--target',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);p.add_argument('--seconds',type=float,default=120)
    a=p.parse_args();search(a.target,a.out,a.seconds)

"""Equivariant ball-replacement compression, preserving the entire fixed circle.

Producer and independent replay use the accepted local ball replacement code.
This creates a new descendant; it never rewrites the Montesinos master.
"""
from pathlib import Path
from itertools import combinations
import json,gzip,sys,time,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from compact3 import Complex,digest
from verify_compact3 import replay_part

def read(p): return json.loads(gzip.decompress(Path(p).read_bytes()))
def save(p,x):
    p=Path(p);p.parent.mkdir(exist_ok=True,parents=True)
    p.write_bytes(gzip.compress((json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))

def build(source,output,seconds=120):
    st=time.monotonic(); raw=Path(source).read_bytes();x=json.loads(gzip.decompress(raw))
    tags=sorted({v for t in x['facets'] for v in t},key=repr)
    ix={v:i for i,v in enumerate(tags)};D=[ix[x['deck'][v]] for v in tags]
    rank=lambda i:min(tags[i],tags[D[i]])
    initial=[tuple(sorted((ix[v] for v in t),key=rank)) for t in x['facets']]
    K={tuple(sorted(t)) for t in initial}
    protected={i for i,d in enumerate(D) if i==d}
    worker=Complex(K,D,protected);steps=worker.run(seconds=seconds)
    ordered=[tuple(sorted(t,key=rank)) for t in sorted(worker.K)]
    p={'tags':tags,'deck':D,'protected':sorted(protected),'initial':sorted(K),
       'initial_ordered':initial,'initial_sha256':digest(sorted(K)),
       'endpoint':sorted(worker.K),'endpoint_sha256':digest(sorted(worker.K)),
       'ordered_endpoint':ordered,'steps':steps}
    print('COARSENED',len(K),'->',len(worker.K),'pairs',len(steps),flush=True)
    srcordered=[tuple(tags[v] for v in t) for t in initial]
    rawdeck={v:x['deck'][v] for v in tags}
    out,receipt=replay_part(p,srcordered,rawdeck,{tags[v] for v in protected})
    active=sorted({v for t in out for v in t});idx={v:i for i,v in enumerate(active)}
    fs=sorted(tuple(sorted(idx[v] for v in t)) for t in out)
    deck=[idx[x['deck'][v]] for v in active]
    y={'schema':'marked-pretzel-source-descendant-v1','facets':fs,'closed_cover_facets':fs,
       'deck':deck,'vertex_ranks':[min(i,deck[i]) for i in range(len(deck))],
       'parent_vertices':active,'parent_source_sha256':digest(x),
       'parent_file_sha256':hashlib.sha256(raw).hexdigest(),
       'ancestry':'equivariant pairs of independently replayed relative 3-ball replacements; fixed circle pointwise preserved',
       'ramification_edges':sorted({e for t in fs for e in combinations(t,2) if all(deck[v]==v for v in e)}),
       'source_interpretation':x['source_interpretation']}
    output=Path(output);output.mkdir(exist_ok=True,parents=True)
    (output/'master_source.json.gz').write_bytes(raw)
    save(output/'coarsening_certificate.json.gz',p)
    save(output/'source.json.gz',y)
    receipt.update({'passed':True,'source_sha256':digest(x),'endpoint_sha256':digest(y),'seconds':time.monotonic()-st})
    (output/'coarsening_verified.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)
    return y

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=ROOT/'seifert/data/source.json.gz')
    p.add_argument('--output',type=Path,default=ROOT/'coarsened');p.add_argument('--seconds',type=int,default=120)
    a=p.parse_args();build(a.source,a.output,a.seconds)

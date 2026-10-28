#!/usr/bin/env python3
"""Bounded native local-link admission on the exact X0 numbered exports.

No simplification, homology, isomorphism signature, or group call is made.
The original immutable source and the candidate master are read-only.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, resource, time
import regina

ROOT = Path(__file__).resolve().parents[1]


def main(which, output=None):
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 125))
    start, cpu = time.monotonic(), time.process_time()
    path = ROOT/'candidate'/f'{which}_regina.json.gz'
    raw = path.read_bytes()
    d = json.loads(gzip.decompress(raw))
    if d['simplex_count'] > 250000:
        raise ValueError('construction cell ceiling exceeded')
    t = regina.Triangulation4.fromGluings(
        d['simplex_count'],
        [(i,k,j,regina.Perm5(*p)) for i,k,j,p in d['gluings']])
    result = {
        'schema':'native-marked-X0-admission-v1',
        'source':path.name,
        'source_sha256':hashlib.sha256(raw).hexdigest(),
        'source_facets_sha256':d['source_facets_sha256'],
        'regina':regina.versionString(),
        'simplex_count':t.size(), 'valid':t.isValid(),
        'orientable':t.isOrientable(), 'connected':t.isConnected(),
        'closed':t.isClosed(), 'ideal':t.isIdeal(),
        'boundary_components':t.countBoundaryComponents(),
        'f_vector':[t.countVertices(),t.countEdges(),t.countTriangles(),
                    t.countTetrahedra(),t.size()],
        'chi':t.eulerCharTri(),
        'vertex_link_scope':'Regina validity, closedness and nonideal predicates recognize all vertex links as PL three-spheres; no homology substitution.',
        'wall_seconds':time.monotonic()-start,
        'cpu_seconds':time.process_time()-cpu,
        'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        'scope':'Closed PL-manifold admission only; no PL standardization.'}
    result['passed'] = (result['valid'] and result['orientable'] and
        result['connected'] and result['closed'] and not result['ideal'] and
        result['boundary_components']==0 and result['chi']==2)
    output = output or ROOT/'diagnostics'/f'current_{which}_native_admission.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:
        raise ValueError('Native admission failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('which',choices=('cover','quotient'))
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    main(args.which,args.output)

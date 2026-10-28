#!/usr/bin/env python3
"""Verify:  G6 (census 6-pentachoron triangulation of X0)  ~PL  standard S^4.            (v1.1)

Certificate = three forward move paths on canonical triangulations:
  (1) G6    -> W  (Pachner 3-3, 2-4 and guarded 2-0 edge/triangle moves)   [key G6_to_X]
  (2) std6  -> W  (same move types)                                        [key std6_to_X]
  (3) S^4_2 -> std6 (Pachner moves only)                                   [key sphere_to_std6]
where S^4_2 = cPkbbbbaaaaaaaa is two 4-simplices with corresponding facets glued by the identity
(the double of a 4-simplex, a standard PL 4-sphere), rebuilt here from scratch.

Each step (a, code, i, b): T = fromIsoSig(a), which must be canonical; face i of the named dimension
must exist; the guarded legality test must accept the move; the move must change the size by the
expected amount; and T.isoSig() must equal b exactly.  Consecutive steps must chain.

Pachner moves preserve PL type.  An admissible edge or triangle 2-0 move (degree-2 face in two
distinct pentachora, doubled-simplex link, no closed chain among the faces identified) preserves
the PL homeomorphism type of the ambient manifold.  (The flattening map itself is not a
homeomorphism; it is the local replacement that preserves PL type.)

v1.1: every check is an explicit test, so `python -O` cannot disable it (v1.0 used assert).
Exit status 0 = pass, 1 = fail.  The certificate file is unchanged from v1.0.
"""
import hashlib
import json
import sys
from collections import Counter

import regina

CERT_SHA256 = 'e9ba9f14831bb0da32043d7d93eb0148f0683e2a6e4f67d5156e8e9fc94adcdf'
G6 = 'gAMLMQaabcdeffff2a+a2avbmaoaxaOb0aUa'
SPHERE2 = 'cPkbbbbaaaaaaaa'
DIM = {'P0': 0, 'P1': 1, 'P2': 2, 'P3': 3, 'P4': 4, 'Z1': 1, 'Z2': 2}
COUNT = {0: 'countVertices', 1: 'countEdges', 2: 'countTriangles', 3: 'countTetrahedra', 4: 'size'}
GET = {0: 'vertex', 1: 'edge', 2: 'triangle', 3: 'tetrahedron', 4: 'pentachoron'}


class Fail(Exception):
    pass


def need(ok, msg):
    if not ok:
        raise Fail(msg)


def double_of_simplex():
    t = regina.Triangulation4()
    a, b = t.newPentachora(2)
    for facet in range(5):
        a.join(facet, b, regina.Perm5())
    return t.isoSig()


def check(steps, start, end, label, allowed):
    need(isinstance(steps, list) and len(steps) > 0, f'{label}: empty or malformed path')
    need(steps[0][0] == start, f'{label}: wrong start')
    need(steps[-1][3] == end, f'{label}: wrong end')
    prev = start
    for k, step in enumerate(steps):
        need(isinstance(step, list) and len(step) == 4, f'{label}: step {k} malformed')
        a, code, i, b = step
        need(a == prev, f'{label}: chain broken at step {k}')
        need(code in allowed, f'{label}: step {k} uses move {code!r}, not allowed on this path')
        T = regina.Triangulation4.fromIsoSig(a)
        need(T.isoSig() == a, f'{label}: step {k} source signature is not canonical')
        d = DIM[code]
        need(type(i) is int and 0 <= i < getattr(T, COUNT[d])(), f'{label}: step {k} face index out of range')
        f = getattr(T, GET[d])(i)
        n0 = T.size()
        if code[0] == 'P':
            need(T.hasPachner(f), f'{label}: step {k} Pachner move not legal')
            need(T.pachner(f), f'{label}: step {k} Pachner move failed')
            delta = 2 * d - 4
        else:
            need(T.has20(f), f'{label}: step {k} 2-0 move not legal')
            need(T.move20(f), f'{label}: step {k} 2-0 move failed')
            delta = -2
        need(T.size() == n0 + delta, f'{label}: step {k} wrong size change')
        need(T.isoSig() == b, f'{label}: step {k} result differs from the recorded next state')
        prev = b
    print(f'PASS {label}: {len(steps)} moves {dict(sorted(Counter(s[1] for s in steps).items()))}', flush=True)


def main():
    if len(sys.argv) != 2:
        print('usage: python3 verify_x0_standard.py x0_standard_certificate.json')
        return 2
    try:
        raw = open(sys.argv[1], 'rb').read()
        need(hashlib.sha256(raw).hexdigest() == CERT_SHA256, 'certificate SHA-256 differs from the published one')
        C = json.loads(raw)
        need(double_of_simplex() == SPHERE2, 'double of a 4-simplex has an unexpected signature')
        need(regina.Example4.sphere().isoSig() == SPHERE2, "Regina's Example4.sphere() disagrees")
        need(regina.Triangulation4.fromIsoSig(G6).isoSig() == G6, 'G6 signature is not canonical')
        mixed = {'P2', 'P3', 'Z1', 'Z2'}
        check(C['G6_to_X'], G6, C['X'], 'G6 -> W', mixed)
        check(C['std6_to_X'], C['std6'], C['X'], 'std6 -> W', mixed)
        check(C['sphere_to_std6'], SPHERE2, C['std6'], 'S^4_2 -> std6 (Pachner only)', {'P0', 'P1', 'P2', 'P3', 'P4'})
    except Fail as e:
        print('FAIL:', e)
        return 1
    print('ALL PASS: G6 is PL homeomorphic to the double of a 4-simplex (standard S^4). '
          'With the accepted X0 -> G6 certificate, X0 is PL standard, hence diffeomorphic to S^4.')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""Check the 594 closed 4D 2-0 moves without Regina's move implementation.

The face-identification guards are checked as forests, and the replacement
facet pairing is rebuilt by tracing across the flattened pillow.  Regina is
used for signature decoding/encoding, simplex construction, and face identity.
No move20, has20, or twoZeroMove call occurs in this program.

Reference: Regina 7.4.1 engine/triangulation/detail/moves-impl.h,
internal20; tag commit 219d658e5e51cc003807fc871fb3534a448a972f.
This checker is deliberately restricted to closed, unlocked 4D inputs.
"""
import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
import regina

CERTIFICATE_SHA256 = 'e9ba9f14831bb0da32043d7d93eb0148f0683e2a6e4f67d5156e8e9fc94adcdf'


def require(test, message):
    if not test:
        raise ValueError(message)


def compose(left, right):
    return tuple(left[right[i]] for i in range(5))


def inverse(p):
    q = [0]*5
    for i, j in enumerate(p):
        q[j] = i
    return tuple(q)


def perm_tuple(p):
    return tuple(p[i] for i in range(5))


def face_id(simplex, vertices):
    vertices = tuple(vertices)
    k = len(vertices)-1
    rest = tuple(x for x in range(5) if x not in vertices)
    p = regina.Perm5(*(vertices+rest))
    cls = (regina.Vertex4, regina.Edge4, regina.Triangle4, regina.Tetrahedron4)[k]
    return simplex.face(k, cls.faceNumber(p)).index()


def check_forest(pairs):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in pairs:
        ra, rb = find(a), find(b)
        require(ra != rb, 'closed chain in faces to be identified')
        parent[ra] = rb


def check_and_reconstruct(sig, code, index, expected):
    k = int(code[1])
    require(k in (1,2), 'unsupported face dimension')
    t = regina.Triangulation4.fromIsoSig(sig)
    require(t.isoSig() == sig, 'input signature not canonical')
    require(t.isValid() and t.isConnected() and not t.hasBoundaryFacets(),
            'requires valid connected closed facet pairing')
    require(not t.hasLocks(), 'locked input')
    f = t.face(k, index)
    require(f.isValid() and not f.isBoundary() and f.degree()==2,
            'face is not valid interior degree two')
    emb = list(f.embeddings())
    p, q = (e.simplex() for e in emb)
    pp, qp = (perm_tuple(e.vertices()) for e in emb)
    require(p.index()!=q.index(), 'both incidences belong to same simplex')
    cross = perm_tuple(p.adjacentGluing(pp[4]))
    require(all(p.adjacentSimplex(pp[j]) == q and
                perm_tuple(p.adjacentGluing(pp[j])) == cross
                for j in range(k+1,5)), 'link is not a double simplex')
    require(all(cross[pp[j]] == qp[j] for j in range(k+1)),
            'central face labelling is inconsistent')

    # For each j-face strictly containing the opposite (3-k)-face,
    # the flattening identifies one pair across the two pillow sides.
    # Requiring a forest excludes all prohibited identification cycles.
    # All boundary-chain exclusions are vacuous for these closed inputs.
    opposite = pp[k+1:]
    for j in range(3-k,4):
        pairs=[]
        for extra in combinations(pp[:k+1], j-(3-k)):
            vertices=opposite+extra
            pairs.append((face_id(p,vertices),
                          face_id(q,tuple(cross[v] for v in vertices))))
        check_forest(pairs)

    removed={p.index(),q.index()}
    pair={}
    for vertex in pp[:k+1]:
        pair[(p.index(),vertex)]=(q.index(),cross[vertex],cross)
        pair[(q.index(),cross[vertex])]=(p.index(),vertex,inverse(cross))
    out=regina.Triangulation4()
    survivors={i:out.newPentachoron() for i in range(t.size()) if i not in removed}
    for original_index, dest in survivors.items():
        original=t.pentachoron(original_index)
        for facet in range(5):
            if dest.adjacentSimplex(facet) is not None:
                continue
            nxt=original.adjacentSimplex(facet)
            require(nxt is not None, 'unexpected boundary')
            accumulated=perm_tuple(original.adjacentGluing(facet))
            nxt_index=nxt.index()
            nxt_facet=accumulated[facet]
            visited=set()
            while nxt_index in removed:
                half=(nxt_index,nxt_facet)
                require(half in pair and half not in visited,
                        'invalid or cyclic flattening traversal')
                visited.add(half)
                other_index, other_facet, across=pair[half]
                accumulated=compose(across,accumulated)
                other=t.pentachoron(other_index)
                nxt=other.adjacentSimplex(other_facet)
                require(nxt is not None, 'unexpected boundary after flattening')
                onward=perm_tuple(other.adjacentGluing(other_facet))
                accumulated=compose(onward,accumulated)
                nxt_index=nxt.index()
                nxt_facet=onward[other_facet]
            require(accumulated[facet]==nxt_facet,'facet-permutation mismatch')
            dest.join(facet,survivors[nxt_index],regina.Perm5(*accumulated))
    require(out.size()==t.size()-2,'wrong size change')
    require(out.isoSig()==expected,'rebuilt output signature mismatch')
    return t.size()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data=args.certificate.read_bytes()
    digest=hashlib.sha256(data).hexdigest()
    require(digest == CERTIFICATE_SHA256, 'certificate SHA-256 mismatch')
    c=json.loads(data)
    counts=Counter(); path_counts={}; max_size=0
    start=time.monotonic()
    for label in ('G6_to_X','std6_to_X','sphere_to_std6'):
        local=Counter()
        for n,(a,code,index,b) in enumerate(c[label]):
            if code in ('Z1','Z2'):
                try:
                    max_size=max(max_size,check_and_reconstruct(a,code,index,b))
                except Exception as e:
                    raise ValueError(f'{label} step {n}: {e}') from e
                counts[code]+=1; local[code]+=1
        path_counts[label]=dict(local)
    require(dict(counts)=={'Z1':498,'Z2':96},'unexpected number of 2-0 moves')
    report={'status':'PASS','certificate_sha256':digest,
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'regina_version':regina.versionString(),
            'two_zero_moves_checked':sum(counts.values()),'counts':dict(counts),
            'by_path':path_counts,'largest_source_size':max_size,
            'seconds':round(time.monotonic()-start,3),
            'regina_move_functions_called':False,
            'guard_method':'closed-input doubled-link check plus forests of identified faces',
            'output_method':'new facet pairing by explicit traversal through flattened pillow',
            'remaining_dependencies':['Regina signature codec','Regina face skeleton and validity','Regina simplex constructor']}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Check all inverse boundary shellings using raw face orbits.

This implements Regina7.4.1's documented shellBoundary conditions for a
previously admitted finite PL4-manifold. All face validity follows from
that admission and is preserved inductively. No4D skeleton is constructed.
An independent raw-facet/vertex-link admission report is REQUIRED and bound
to the exact initial layered RGA. No cached Regina property is changed.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import regina


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def orbit(tri, pentachoron, vertices, targets=None):
    """Find all embeddings of a local face, and whether it meets real boundary."""
    start = (pentachoron, tuple(sorted(vertices)))
    seen, stack = {start}, [start]
    boundary = False
    touched = set()
    while stack:
        index, face = stack.pop()
        if targets is not None and (index, face) in targets:
            touched.add((index, face))
        pent = tri.simplex(index)
        for facet in range(5):
            if facet in face:
                continue
            adjacent = pent.adjacentSimplex(facet)
            if adjacent is None:
                boundary = True
                continue
            p = pent.adjacentGluing(facet)
            next_face = (adjacent.index(), tuple(sorted(p[v] for v in face)))
            if next_face not in seen:
                seen.add(next_face)
                stack.append(next_face)
    return seen, boundary, touched


def shell_guard(tri, index):
    pent = tri.simplex(index)
    require(not pent.isLocked() and not any(pent.isFacetLocked(f) for f in range(5)),
            'locked simplex/facet in candidate shelling')
    boundary_facets = tuple(f for f in range(5) if pent.adjacentSimplex(f) is None)
    require(1 <= len(boundary_facets) <= 4, 'shell must have between one and four boundary facets')
    # f is opposite the intersection of all boundary facets. Its vertices
    # are precisely the indices of those boundary facets.
    middle = boundary_facets
    seen, boundary, _ = orbit(tri, index, middle)
    require(not boundary, 'the opposite face lies in the triangulation boundary')
    cofaces = {(index, tuple(sorted(middle + (v,))))
               for v in range(5) if v not in middle}
    unchecked = set(cofaces)
    traversed = len(seen)
    while unchecked:
        source = next(iter(unchecked))
        embeddings, _, touched = orbit(tri, source[0], source[1], cofaces)
        require(len(touched) == 1, 'two cofaces containing the opposite face are identified')
        unchecked -= touched
        traversed += len(embeddings)
    return {'boundary_facets': len(boundary_facets), 'local_face_embeddings_examined': traversed}


def verify(source, layered, admission, output):
    source, layered, admission, output = map(Path, (source, layered, admission, output))
    require(importlib.metadata.version('regina') == '7.4.1',
            'requires Regina distribution 7.4.1')
    source_hash, layered_hash, admission_hash = sha(source), sha(layered), sha(admission)
    proof = json.loads(admission.read_text())
    require(sha(admission) == admission_hash, 'admission file changed while it was loaded')
    require(proof['status'] == 'PASS_FINITE_ORIENTABLE_PL_4_MANIFOLD',
            'requires independent finite PL manifold admission')
    require(proof['source_sha256'] == layered_hash, 'admission is bound to another layered file')
    require(proof['all_face_stabilizers_trivial'] and all(x['recognised'] for x in proof['vertex_links']),
            'incomplete manifold admission')
    tri = regina.Triangulation4(regina.open(str(layered)))
    original = regina.Triangulation4(regina.open(str(source)))
    require(sha(source) == source_hash and sha(layered) == layered_hash,
            'a triangulation file changed while it was loaded')
    require(tri.size() >= original.size(), 'layered input is smaller than source')
    require(hashlib.sha256(tri.tightEncoding().encode()).hexdigest() == proof['source_tight_sha256'],
            'admitted labelled state differs')
    start, count, total_embeddings = time.monotonic(), 0, 0
    histogram = {}
    while tri.size() > original.size():
        index = tri.size() - 1
        check = shell_guard(tri, index)
        total_embeddings += check['local_face_embeddings_examined']
        k = check['boundary_facets']
        histogram[k] = histogram.get(k, 0) + 1
        tri.removeSimplexAt(index)
        count += 1
        if count % 1000 == 0:
            print(json.dumps({'checked_inverse_shellings': count, 'remaining_pentachora': tri.size(),
                              'seconds': time.monotonic() - start}), flush=True)
    require(tri.tightEncoding() == original.tightEncoding(),
            'checked inverse shellings did not recover exact source')
    require(sha(source) == source_hash and sha(layered) == layered_hash and
            sha(admission) == admission_hash, 'an input file changed during verification')
    result = {'status': 'PASS_ALL_RAW_INVERSE_SHELLING_GUARDS', 'moves': count,
              'source_sha256': source_hash, 'layered_sha256': layered_hash,
              'admission_sha256': admission_hash, 'source_recovered_exactly': True,
              'regina_distribution': '7.4.1', 'input_hashes_stable_during_verification': True,
              'boundary_facet_counts': histogram, 'face_embeddings_examined': total_embeddings,
              'initial_finite_PL_manifold_admitted': True,
              'all_face_validity_preserved_inductively': True,
              'four_dimensional_skeleton_used': False, 'cached_flags_modified': False,
              'guard_basis': 'Regina7.4.1 Triangulation4.shellBoundary documented conditions; validity supplied by the independent finite PL manifold admission.',
              'guard_documentation': 'https://regina-normal.github.io/engine-docs/classregina_1_1Triangulation_3_014_01_4.html',
              'seconds': time.monotonic() - start, 'checker_sha256': sha(__file__)}
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--layered', required=True)
    parser.add_argument('--admission', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    verify(args.source, args.layered, args.admission, args.out)

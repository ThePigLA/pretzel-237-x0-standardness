#!/usr/bin/env python3
"""Check a saved layering as a sequence of PL collar additions.

This is a conditional PL certificate: the source must already be a finite
PL manifold. It never recognizes a four-dimensional vertex link. Every
attaching abstract 3-ball is checked to have embedded open interior with
no boundary-to-interior identifications. Its boundary may self-identify.
The companion COLLAR_LEMMA.md gives the vanishing-height collar argument.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import itertools
import json
from pathlib import Path
import time

import regina

from audit_fast_closed_gluing import boundary_pairings


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tight_sha(t):
    return hashlib.sha256(t.tightEncoding().encode()).hexdigest()


def face(t, d, i):
    return t.simplex(i) if d == 3 else t.face(d, i)


def read_boundary(t):
    places, pairs, _ = boundary_pairings(t)
    b = regina.Triangulation3()
    b.newSimplices(len(places))
    for (i, f), (j, p) in pairs.items():
        if b.simplex(i).adjacentSimplex(f) is None:
            b.simplex(i).join(f, b.simplex(j), regina.Perm4(*p))
    return b, places


def check_support(t, labels, d, i):
    """Directly test injectivity on the interior of the abstract ball."""
    old = set(labels)
    centre = set.intersection(*(set(row) for row in labels.values()))
    missing = {j: next(iter(set(range(5)) - set(row))) for j, row in labels.items()}
    old_missing = set(missing.values())
    require(len(old) == 4 - d and len(old_missing) == len(old),
            'wrong number of distinct old tetrahedra/facets')
    require(len(centre) == d + 1 and centre == set(range(5)) - old_missing,
            'the old front is not the expected join ball')
    target = face(t, d, i)
    star = {i} if d == 3 else {e.simplex().index() for e in target.embeddings()}
    require(star == old, 'old tetrahedra do not exhaust the Pachner star')
    if d < 3:
        require(target.degree() == len(old), 'centre face has a repeated embedding')
    by_missing = {m: j for j, m in missing.items()}
    for j, row in labels.items():
        require(len(set(row)) == 4 and set(row) == set(range(5)) - {missing[j]},
                'old tetrahedron has incorrect abstract vertex labels')
        s = t.simplex(j)
        for f in range(4):
            if row[f] not in old_missing:
                continue
            other = s.adjacentSimplex(f)
            require(other is not None and other.index() == by_missing[row[f]],
                    'an intrinsic old-front triangle has the wrong neighbour')
            p = s.adjacentGluing(f)
            for k in range(4):
                expected = missing[j] if k == f else row[k]
                require(labels[other.index()][p[k]] == expected,
                        'an intrinsic old-front triangle has the wrong vertex map')

    # The relative interiors of precisely those faces containing the
    # entire centre simplex form the interior of centre * boundary(link).
    # Global face IDs and ordered vertex maps detect both unwanted face
    # identifications and nontrivial identifications within an open face.
    interior = {}
    interior_occurrences = {}
    interior_degrees = {}
    boundary = set()
    checked = 0
    for j, row in labels.items():
        s = t.simplex(j)
        for fd, count in [(0, 4), (1, 6), (2, 4)]:
            for local in range(count):
                p = s.faceMapping(fd, local)
                ordered = tuple(row[p[k]] for k in range(fd + 1))
                subset = frozenset(ordered)
                global_face = s.face(fd, local)
                require(global_face.isValid(), 'support contains an invalid global face')
                key = (fd, global_face.index())
                if centre.issubset(subset):
                    require(key not in boundary, 'support boundary meets its open interior')
                    prior = interior.setdefault(key, ordered)
                    require(prior == ordered, 'the abstract support interior is identified with itself')
                    interior_occurrences[key] = interior_occurrences.get(key, 0) + 1
                    interior_degrees[key] = global_face.degree()
                    checked += 1
                else:
                    require(key not in interior, 'support boundary meets its open interior')
                    boundary.add(key)
    require(interior_occurrences == interior_degrees,
            'an interior support face has an embedding outside the recorded support')
    # Top-dimensional interiors are embedded and mutually distinct by the
    # distinct-old-tetrahedra check. The dictionary also verifies that the
    # global canonical vertex order agrees on every repeated occurrence.
    return centre, missing, checked


def expected_front(t, labels, centre, missing):
    """Boundary of an attached abstract 4-simplex, computed locally.

    New tetrahedra are indexed by their missing abstract vertex. A local
    triangle is indexed by its omitted abstract vertex. Each complete map
    includes the omitted vertex, exactly as a Regina facet permutation.
    """
    old = set(labels)
    by_missing = {m: j for j, m in missing.items()}
    expected = {}
    for a in centre:
        for b in set(range(5)) - {a}:
            if b in centre:
                expected[a, b] = ('new', b,
                                  {x: a if x == b else x for x in set(range(5)) - {a}})
                continue
            j = by_missing[b]
            row = labels[j]
            omit = row.index(a)
            s = t.simplex(j)
            other = s.adjacentSimplex(omit)
            require(other is not None, 'the source 3D boundary is not closed')
            p = s.adjacentGluing(omit)
            if other.index() not in old:
                mapping = {x: p[omit] if x == b else p[row.index(x)]
                           for x in set(range(5)) - {a}}
                expected[a, b] = ('old', other.description(), mapping)
            else:
                target_row = labels[other.index()]
                target_a = target_row[p[omit]]
                require(target_a in centre,
                        'an exterior support triangle maps to an interior one')
                target_b = missing[other.index()]
                mapping = {x: target_b if x == b else target_row[p[row.index(x)]]
                           for x in set(range(5)) - {a}}
                expected[a, b] = ('new', target_a, mapping)
    return expected


def match_new_front(t, new_tets, centre, expected):
    """Solve the small front isomorphism, fixing all unaffected tetrahedra.

    This deliberately does not call the producer's label-transport code.
    At most four tetrahedra and 96 possible first-tetrahedron maps occur.
    """
    new_by_index = {s.index(): s for s in new_tets}
    first = new_tets[0].index()
    for missing in sorted(centre):
        for initial in itertools.permutations(sorted(set(range(5)) - {missing})):
            assigned = {first: (missing, initial)}
            todo = [first]
            okay = True
            while todo and okay:
                index = todo.pop()
                a, row = assigned[index]
                s = new_by_index[index]
                for f in range(4):
                    kind, destination, mapping = expected[a, row[f]]
                    other = s.adjacentSimplex(f)
                    if other is None:
                        okay = False
                        break
                    p = s.adjacentGluing(f)
                    if kind == 'old':
                        if other.index() in new_by_index or other.description() != destination:
                            okay = False
                            break
                        if any(p[k] != mapping[row[k]] for k in range(4)):
                            okay = False
                            break
                    else:
                        if other.index() not in new_by_index:
                            okay = False
                            break
                        target_row = [None] * 4
                        for k in range(4):
                            target_row[p[k]] = mapping[row[k]]
                        value = (destination, tuple(target_row))
                        if other.index() in assigned:
                            if assigned[other.index()] != value:
                                okay = False
                                break
                        else:
                            assigned[other.index()] = value
                            todo.append(other.index())
            if okay and len(assigned) == len(new_tets):
                if {v[0] for v in assigned.values()} == centre:
                    return assigned
    raise ValueError('the new 3D Pachner front does not match the exact 4D attachment')


def audit(source_path, path_file, layered_path, events_file, output, limit=None):
    now = time.monotonic()
    source_path, path_file, layered_path, events_file, output = map(
        Path, [source_path, path_file, layered_path, events_file, output])
    package_version = importlib.metadata.version('regina')
    require(package_version == '7.4.1', 'this labelled replay requires regina==7.4.1')
    input_paths = {'source': source_path, 'path': path_file,
                   'layered': layered_path, 'events': events_file}
    input_hashes = {key: sha(path) for key, path in input_paths.items()}
    checker_hash = sha(__file__)
    source = regina.Triangulation4(regina.open(str(source_path)))
    layered = regina.Triangulation4(regina.open(str(layered_path)))
    data = json.loads(path_file.read_text())
    require(all(sha(path) == input_hashes[key] for key, path in input_paths.items()),
            'an input changed while it was being loaded')
    require(layered.size() == source.size() + len(data['moves']), 'layer count differs')
    # Compare every original pairing, without invoking a 4D skeleton.
    for s in source.simplices():
        out_s = layered.simplex(s.index())
        for f in range(5):
            old = s.adjacentSimplex(f)
            actual = out_s.adjacentSimplex(f)
            if old is None:
                require(actual is None or actual.index() >= source.size(),
                        'an original boundary facet was paired inside the source')
            else:
                require(actual is not None and actual.index() == old.index() and
                        all(out_s.adjacentGluing(f)[k] == s.adjacentGluing(f)[k] for k in range(5)),
                        'an original internal gluing was changed')
    actual, places = read_boundary(source)
    t = regina.Triangulation3.tightDecoding(data['start_tight'])
    require(t.isValid() and t.isClosed() and t.isConnected(), 'invalid source boundary')
    iso = t.isIsomorphicTo(actual)
    require(iso is not None, 'path does not begin at the actual source boundary')
    locations = {}
    serial = 0
    for s in t.simplices():
        tag = str(serial); serial += 1
        s.setDescription(tag)
        pent, facet, vertices = places[iso.simpImage(s.index())]
        p = iso.facetPerm(s.index())
        locations[tag] = (pent, facet, tuple(vertices[p[k]] for k in range(4)))
    # This is an explicit initial boundary isomorphism. Every subsequent
    # saved attachment must agree with this map and its checked transport.
    interior_checks = 0
    move_counts = {str(d): 0 for d in range(4)}
    with gzip.open(events_file, 'rt') as events:
        for step, move in enumerate(data['moves']):
            if limit is not None and step >= limit:
                break
            raw = events.readline()
            require(raw, 'event log ends early')
            event = json.loads(raw)
            require(event['step'] == step and event['dimension'] == move['dimension'] and
                    event['face'] == move['face'], 'event and 3D path disagree')
            require(tight_sha(t) == move['before_sha256'], '3D source digest differs')
            require(event['source_boundary_sha256'] == move['before_sha256'] and
                    event['target_boundary_sha256'] == move['after_sha256'], 'event digest differs')
            d, i = move['dimension'], move['face']
            require(t.hasPachner(face(t, d, i)), '3D Pachner guard failed')
            source_n = source.size() + step
            require(event['added_pentachoron'] == source_n, 'wrong added pentachoron')
            by_place = {(p, f): (tag, vertices) for tag, (p, f, vertices) in locations.items()}
            tag_to_index = {s.description(): s.index() for s in t.simplices()}
            labels = {}
            old_tags = set()
            new_facets = set()
            for g in event['gluings']:
                key = (g['old_pentachoron'], g['old_facet'])
                require(key in by_place, 'attachment is not on the current boundary map')
                tag, vertices = by_place[key]
                p = g['permutation']
                require(sorted(p) == list(range(5)), 'invalid saved gluing permutation')
                missing = g['new_facet']
                require(p[missing] == key[1], 'saved facet does not map to the claimed facet')
                row = tuple(p.index(v) for v in vertices)
                require(set(row) == set(range(5)) - {missing}, 'invalid reconstructed labels')
                index = tag_to_index[tag]
                require(index not in labels and missing not in new_facets, 'repeated attaching facet')
                labels[index] = row
                old_tags.add(tag); new_facets.add(missing)
                s = layered.simplex(source_n)
                other = s.adjacentSimplex(missing)
                require(other is not None and other.index() == key[0] and
                        all(s.adjacentGluing(missing)[k] == p[k] for k in range(5)),
                        'saved event does not equal actual final 4D attachment')
            s = layered.simplex(source_n)
            observed_old = {f for f in range(5)
                            if s.adjacentSimplex(f) is not None and
                            s.adjacentSimplex(f).index() < source_n}
            require(new_facets == observed_old, 'an initial 4D attachment was omitted from the log')
            require(all(s.adjacentSimplex(f) is None or s.adjacentSimplex(f).index() != source_n
                        for f in range(5)), 'new pentachoron is self-glued')
            centre, missing, checks = check_support(t, labels, d, i)
            interior_checks += checks
            expected = expected_front(t, labels, centre, missing)
            require(t.pachner(face(t, d, i)), '3D Pachner execution failed')
            require(tight_sha(t) == move['after_sha256'], '3D target digest differs')
            new_tets = [s for s in t.simplices() if not s.description()]
            require(len(new_tets) == d + 1, 'wrong number of new boundary tetrahedra')
            assigned = match_new_front(t, new_tets, centre, expected)
            for tag in old_tags:
                del locations[tag]
            for s in new_tets:
                missing_new, row = assigned[s.index()]
                tag = str(serial); serial += 1
                s.setDescription(tag)
                locations[tag] = (source_n, missing_new, row)
            move_counts[str(d)] += 1
            if (step + 1) % 1000 == 0:
                print(json.dumps({'checked_collar_additions': step + 1,
                                  'boundary_tetrahedra': t.size(),
                                  'seconds': time.monotonic() - now}), flush=True)
        total = sum(move_counts.values())
        complete = total == len(data['moves'])
        if complete:
            require(not events.readline(), 'extra events after the 3D path ends')
            require(t.tightEncoding() == data['end_tight'], 'wrong final labelled boundary')
            final_boundary, final_places = read_boundary(layered)
            require(t.isIsomorphicTo(final_boundary) is not None, 'wrong actual final boundary')
            require({(p, f) for p, f, _ in locations.values()} ==
                    {(p, f) for p, f, _ in final_places}, 'final boundary facet coverage differs')
    require(all(sha(path) == input_hashes[key] for key, path in input_paths.items()),
            'an input changed during the replay')
    require(sha(__file__) == checker_hash, 'checker changed during the replay')
    report = {
        'status': 'PASS_ALL_LOCAL_COLLAR_ATTACHMENT_HYPOTHESES' if complete else 'PASS_PREFIX_LOCAL_COLLAR_HYPOTHESES',
        'source_sha256': sha(source_path), 'source_tight_sha256': tight_sha(source),
        'path_sha256': sha(path_file), 'layered_sha256': sha(layered_path),
        'layered_tight_sha256': tight_sha(layered), 'events_sha256': sha(events_file),
        'source_pentachora': source.size(), 'layered_pentachora': layered.size(),
        'checked_additions': total, 'move_dimension_counts': move_counts,
        'interior_face_occurrences_checked': interior_checks,
        'all_original_internal_gluings_unchanged': True,
        'each_saved_attachment_matches_actual_4d_facet_permutations': True,
        'support_open_interior_embedded_and_disjoint_from_boundary_image': True,
        'all_embeddings_of_each_interior_support_face_exhausted': True,
        'intrinsic_old_and_new_front_pairings_checked': True,
        'boundary_map_propagated_by_independent_local_front_solver': True,
        'final_raw_boundary_checked': complete,
        'source_finite_pl_manifold_premise_required': True,
        'four_dimensional_skeleton_used': False,
        'vertex_link_recognition_performed': False,
        'regina_package_version': package_version,
        'input_hashes_stable_before_loading_after_loading_and_at_end': True,
        'standardness_asserted': False,
        'seconds': time.monotonic() - now, 'checker_sha256': checker_hash}
    temporary = output.with_suffix(output.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(output)
    print(json.dumps(report), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--path', required=True)
    parser.add_argument('--layered', required=True)
    parser.add_argument('--events', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--limit', type=int)
    a = parser.parse_args()
    audit(a.source, a.path, a.layered, a.events, a.out, a.limit)

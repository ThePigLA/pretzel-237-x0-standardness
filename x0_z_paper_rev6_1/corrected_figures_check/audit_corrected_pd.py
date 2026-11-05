#!/usr/bin/env python3
"""Audit the exact two-crossing corrections and the empty slide strip.

This uses only standard Python. It checks the four supplied PD records,
not their transcription from images, and does not generate a Kirby move
certificate. It verifies the combinatorial premises of the geometric
handle-slide argument in Lemma 7.1.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_pd(path):
    data = json.loads(path.read_text())
    pd = data['pd_code']
    labels = data['component_labels']
    require(data['order'] == ['C' + str(i) for i in range(8)], 'unexpected component order')
    require(data['decorations'] == [0, 'x', 'x', -1, 'x', 0, 0, 0], 'unexpected framings or dotted roles')
    require(all(len(row) == 4 and all(type(x) is int for x in row) for row in pd), 'bad crossing row')
    owner = {}
    successor = {}
    for i, cycle in enumerate(labels):
        require(len(cycle) > 1, 'empty or trivial component label cycle')
        for j, x in enumerate(cycle):
            require(x not in owner, 'component label occurs in two cycles')
            owner[x] = i
            successor[x] = cycle[(j + 1) % len(cycle)]

    occurrences = defaultdict(list)
    sums = [[0] * len(labels) for _ in labels]
    for i, row in enumerate(pd):
        a, b, c, d = row
        require(owner[a] == owner[c] and owner[b] == owner[d], 'opposite half-edges do not follow a component')
        require((successor[a] == c) != (successor[c] == a), 'understrand orientation not determined')
        require((successor[b] == d) != (successor[d] == b), 'overstrand orientation not determined')
        under = 1 if successor[a] == c else -1
        over = 1 if successor[b] == d else -1
        if owner[a] != owner[b]:
            sums[owner[a]][owner[b]] += under * over
            sums[owner[b]][owner[a]] += under * over
        for j, edge in enumerate(row):
            occurrences[edge].append((i, j))
    require(set(occurrences) == set(owner), 'PD edge set differs from declared components')
    require(all(len(darts) == 2 for darts in occurrences.values()), 'edge does not have two ends')
    require(all(x % 2 == 0 for row in sums for x in row), 'non-integral linking number')
    linking = [[x // 2 for x in row] for row in sums]

    # A dart starts at a crossing and follows its edge. Traverse an edge,
    # then take the preceding half-edge in the recorded cyclic order.
    mate = {}
    for a, b in occurrences.values():
        mate[a] = b
        mate[b] = a
    todo = set(mate)
    faces = []
    while todo:
        start = min(todo)
        current = start
        cycle = []
        while True:
            require(current in todo, 'face cycle revisits a previously consumed dart')
            todo.remove(current)
            cycle.append(current)
            across = mate[current]
            current = (across[0], (across[1] - 1) % 4)
            if current == start:
                break
        faces.append(cycle)

    # Connected projection and Euler characteristic 2 certify that the
    # supplied rotation system is a cellular embedding in a sphere.
    adjacency = defaultdict(set)
    for darts in occurrences.values():
        i, j = darts[0][0], darts[1][0]
        adjacency[i].add(j)
        adjacency[j].add(i)
    reached = {0}
    stack = [0]
    while stack:
        for nxt in adjacency[stack.pop()]:
            if nxt not in reached:
                reached.add(nxt)
                stack.append(nxt)
    require(len(reached) == len(pd), 'projection graph disconnected')
    require(len(pd) - len(occurrences) + len(faces) == 2, 'rotation system is not a sphere')
    strip = [cycle for cycle in faces if len(cycle) == 4 and {i for i, _ in cycle} == {1, 3, 13, 29}]
    require(len(strip) == 1, 'expected empty quadrilateral face not found uniquely')
    strip = strip[0]
    strip_labels = [pd[i][j] for i, j in strip]
    require({owner[x] for x in strip_labels} == {0, 1, 2, 3}, 'strip boundary has wrong components')
    require(linking[0][3] == 0, 'mu has nonzero linking with the minus-one component')
    require(all(linking[0][j] == 0 for j in (3, 5, 6, 7)), 'mu links another framed component')
    return data, {
        'path': str(path),
        'sha256': sha(path),
        'crossings': len(pd),
        'projection_vertices': len(pd),
        'projection_edges': len(occurrences),
        'projection_faces': len(faces),
        'projection_euler_characteristic': 2,
        'empty_quadrilateral_darts': [list(d) for d in strip],
        'empty_quadrilateral_edge_labels': strip_labels,
        'empty_quadrilateral_components': [owner[x] for x in strip_labels],
        'signed_linking_matrix': linking,
        'mu_class_in_dotted_basis_C1_C2_C4': [linking[0][j] for j in (1, 2, 4)],
        'minus_one_class_in_dotted_basis_C1_C2_C4': [linking[3][j] for j in (1, 2, 4)],
        'mu_links_other_framed_components_zero': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    pairs = []
    for figure, expected_crossings, expected_edges in [(1, 64, [3, 45, 77, 17]), (2, 60, [3, 43, 73, 17])]:
        old, a = read_pd(args.directory / f'old_fig{figure}.json')
        new, b = read_pd(args.directory / f'new_fig{figure}.json')
        require(a['crossings'] == b['crossings'] == expected_crossings, 'unexpected crossing count')
        require(old['component_labels'] == new['component_labels'], 'component labels changed')
        require(old['decorations'] == new['decorations'], 'decorations changed')
        changed = []
        for i, (row, target) in enumerate(zip(old['pd_code'], new['pd_code'])):
            if row != target:
                shifts = [k for k in range(4) if row[k:] + row[:k] == target]
                require(len(shifts) == 1 and shifts[0] % 2 == 1, 'change is not exactly a crossing switch')
                changed.append({'row_zero_based': i, 'old': row, 'new': target, 'cyclic_shift': shifts[0]})
        require([x['row_zero_based'] for x in changed] == [13, 29], 'unexpected changed crossings')
        require(a['empty_quadrilateral_edge_labels'] == b['empty_quadrilateral_edge_labels'] == expected_edges, 'strip differs')
        mu = a['mu_class_in_dotted_basis_C1_C2_C4']
        before = a['minus_one_class_in_dotted_basis_C1_C2_C4']
        after = b['minus_one_class_in_dotted_basis_C1_C2_C4']
        require(mu == [-1, 1, 0] and before == [0, 1, 0] and after == [1, 0, 0], 'unexpected attaching homology classes')
        require([before[i] - mu[i] for i in range(3)] == after, 'homology does not match subtracting mu')
        for i in range(8):
            for j in range(8):
                if {i, j} not in ({1, 3}, {2, 3}):
                    require(a['signed_linking_matrix'][i][j] == b['signed_linking_matrix'][i][j], 'other linking number changed')
        pairs.append({'figure': figure, 'old': a, 'new': b, 'changed_crossings': changed,
                      'slide_homology': 'C3_new = C3_old - C0',
                      'framing_formula': '-1 + 0 - 2*0 = -1',
                      'other_linking_numbers_unchanged': True})
    result = {'status': 'PASS_TWO_CROSSING_CORRECTION_EMPTY_STRIP_AND_LINKING_CHECKS',
              'schema': 'x0_z.corrected_pd_local_audit.v1',
              'script_sha256': sha(Path(__file__)),
              'scope': 'Four supplied PD transcriptions only. Does not re-extract diagrams from PDF or certify a full standardization.',
              'method': 'Exact integer combinatorics of labelled PD records; standard Python only.',
              'sign_convention': 'Component directions follow declared cyclic arc labels. The chosen planar crossing sign may differ by a common global sign from another PD convention; all zero-linking and slide-consistency statements are unchanged.',
              'geometric_lemma': 'The empty quadrilateral verifies the band-strip premise of the geometric handle-slide proof. The handle-slide isotopy remains the mathematical local argument in Lemma 7.1.',
              'pairs': pairs}
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()

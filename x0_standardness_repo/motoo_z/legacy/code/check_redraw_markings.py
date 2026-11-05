#!/usr/bin/env python3
"""Compare framed link matrices up to component permutation and orientation.

This tests a necessary condition for a marked-link isotopy. A matching matrix
does not prove link isotopy or identify a four-manifold.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

from spherogram import Link


def load(path):
    path = Path(path)
    data = json.loads(path.read_text())
    link = Link(data['pd_code'], check_planarity=True)
    matrix = [[int(x) for x in row] for row in link.linking_matrix()]
    decorations = data['decorations']
    if len(decorations) != len(matrix):
        raise ValueError('Decoration count differs from PD component count')
    identities = list(range(len(matrix)))
    if 'pd_component_to_source' in data and 'components' in data:
        order = data['pd_component_to_source']
        identities = [data['components'][i]['id'] for i in order]
        if 'linking_matrix' in data:
            saved = data['linking_matrix']
            if matrix != [[saved[i][j] for j in order] for i in order]:
                raise ValueError('Saved source linking matrix differs from fresh PD calculation')
    elif 'original_source_ids_in_pd_order' in data:
        identities = data['original_source_ids_in_pd_order']
    return data, matrix, decorations, identities


def witnesses(old_matrix, old_decorations, new_matrix, new_decorations):
    n = len(old_matrix)
    if n != len(new_matrix):
        return []
    if n > 9:
        raise ValueError('This focused exhaustive matcher supports at most nine components')
    result = []
    for perm in itertools.permutations(range(n)):
        if any(new_decorations[i] != old_decorations[perm[i]] for i in range(n)):
            continue
        if any(abs(new_matrix[i][j]) != abs(old_matrix[perm[i]][perm[j]])
               for i in range(n) for j in range(n)):
            continue
        signs = [None] * n
        valid = True
        for start in range(n):
            if signs[start] is not None:
                continue
            signs[start] = 1
            todo = [start]
            while todo and valid:
                i = todo.pop()
                for j in range(n):
                    a, b = new_matrix[i][j], old_matrix[perm[i]][perm[j]]
                    if not b:
                        continue
                    needed = (a // b) * signs[i]
                    if signs[j] is None:
                        signs[j] = needed
                        todo.append(j)
                    elif signs[j] != needed:
                        valid = False
                        break
            if not valid:
                break
        if valid:
            if any(new_matrix[i][j] != signs[i] * signs[j] * old_matrix[perm[i]][perm[j]]
                   for i in range(n) for j in range(n)):
                raise RuntimeError('Internal witness verification failed')
            result.append({'new_pd_component_to_old_pd_component': list(perm),
                           'orientation_signs': signs})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--redraw', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    _, old, old_dec, old_ids = load(args.original)
    _, new, new_dec, _ = load(args.redraw)
    matches = witnesses(old, old_dec, new, new_dec)
    for row in matches:
        row['new_pd_component_to_old_source_id'] = [
            old_ids[i] for i in row['new_pd_component_to_old_pd_component']]
    report = {
        'status': 'PASS_LINKING_NECESSARY_CONDITION' if matches else 'LINKING_MISMATCH',
        'original_sha256': hashlib.sha256(args.original.read_bytes()).hexdigest(),
        'redraw_sha256': hashlib.sha256(args.redraw.read_bytes()).hexdigest(),
        'old_matrix_in_pd_order': old,
        'new_matrix_in_pd_order': new,
        'old_decorations_in_pd_order': old_dec,
        'new_decorations_in_pd_order': new_dec,
        'witnesses': matches,
        'orientation_witnesses_normalized_on_each_linking_graph_component': True,
        'link_isotopy_proved': False,
        'four_manifold_identified': False,
        'scope': 'Necessary framed-link invariant only; source geometry and crossing transcription require separate checks.',
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise ValueError('Choose a new output file')
    args.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'witnesses': len(matches), 'out': str(args.out)}))


if __name__ == '__main__':
    main()

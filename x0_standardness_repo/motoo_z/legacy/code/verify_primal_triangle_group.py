#!/usr/bin/env python3
"""Replay short primal-relator eliminations without importing the solver.

Triangle relators are regenerated through their first tetrahedron embedding
and tetrahedron edge mappings, independently of the solver's Triangle3
edgeMapping calls. Substitutions are verified one by one; no search is done.
"""
from __future__ import annotations
import argparse
from array import array
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


def reconstruct(tri, tree_edges):
    parent = list(range(tri.countVertices()))
    def root(x):
        while x != parent[x]:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree = set(tree_edges)
    require(len(tree) == len(tree_edges) == tri.countVertices() - 1,
            'wrong number of tree edges')
    for index in tree_edges:
        require(type(index) is int and 0 <= index < tri.countEdges(), 'invalid tree edge')
        edge = tri.edge(index)
        a, b = root(edge.vertex(0).index()), root(edge.vertex(1).index())
        require(a != b, 'tree has a cycle')
        parent[a] = b
    require(len({root(x) for x in range(tri.countVertices())}) == 1,
            'tree is not spanning')
    variables = [0] * tri.countEdges()
    generators = 0
    for i in range(tri.countEdges()):
        if i not in tree:
            generators += 1
            variables[i] = generators
    # For each tetrahedron, lookup[a,b] is the signed global edge variable
    # corresponding to the oriented local edge a->b.
    lookup = array('i', [0]) * (16 * tri.size())
    for tetrahedron in tri.tetrahedra():
        offset = 16 * tetrahedron.index()
        pairs = set()
        for local_edge in range(6):
            perm = tetrahedron.edgeMapping(local_edge)
            a, b = perm[0], perm[1]
            pairs.add(tuple(sorted((a, b))))
            g = variables[tetrahedron.edge(local_edge).index()]
            lookup[offset + 4 * a + b] = g
            lookup[offset + 4 * b + a] = -g
        require(len(pairs) == 6, 'inconsistent tetrahedron edge mappings')
    relators = []
    digest = hashlib.sha256()
    for triangle in tri.triangles():
        emb = triangle.embedding(0)
        perm = emb.vertices()
        offset = 16 * emb.simplex().index()
        a, b, c = perm[0], perm[1], perm[2]
        row = tuple(x for x in (lookup[offset + 4*a+b],
                               lookup[offset + 4*b+c],
                               lookup[offset + 4*c+a]) if x)
        relators.append(row)
        digest.update((json.dumps(row, separators=(',', ':')) + '\n').encode())
    return generators, relators, digest.hexdigest()


def verify(source, directory, output):
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina7.4.1')
    source, directory, output = map(Path, (source, directory, output))
    started = time.monotonic()
    meta = json.loads((directory / 'SOURCE.json').read_text())
    result = json.loads((directory / 'RESULT.json').read_text())
    source_hash = sha(source)
    require(source_hash == meta['source_sha256'] == result['source_sha256'],
            'source checksum mismatch')
    packet = regina.open(str(source))
    require(isinstance(packet, regina.Triangulation3), 'not a3D triangulation')
    tri = regina.Triangulation3(packet)
    require(tri.isValid() and not tri.isIdeal() and tri.isConnected(),
            'source is not a valid finite connected3-manifold')
    generators, relators, digest = reconstruct(tri, meta['tree_edges'])
    require(generators == meta['generators'] == result['generators'],
            'generator count differs')
    require(digest == meta['presentation_sha256'] == result['presentation_sha256'],
            'independently regenerated triangle presentation differs')
    substitution = {}
    def image(g):
        current, sign = g, 1
        path = []
        while current in substitution:
            path.append((current, sign))
            current, step_sign = substitution[current]
            sign *= step_sign
        for old, prefix_sign in path:
            substitution[old] = (current, sign * prefix_sign)
        return current, sign
    def reduced(raw):
        answer = []
        for literal in raw:
            g, sign = image(abs(literal))
            if not g:
                continue
            token = g * sign * (1 if literal > 0 else -1)
            if answer and answer[-1] == -token:
                answer.pop()
            else:
                answer.append(token)
        while len(answer) > 1 and answer[0] == -answer[-1]:
            answer = answer[1:-1]
        return answer
    log = directory / 'eliminations.jsonl'
    require(sha(log) == result['elimination_log_sha256'], 'elimination log checksum mismatch')
    steps = 0
    for line in log.read_text().splitlines():
        row = json.loads(line)
        require(row['step'] == steps, 'nonsequential proof step')
        index = row['triangle']
        require(type(index) is int and 0 <= index < len(relators), 'invalid triangle')
        word = reduced(relators[index])
        require(word == row['reduced_relator'], 'recorded reduced relator differs')
        child, ancestor, sign = row['child'], row['parent'], row['sign']
        require(type(child) is int and 1 <= child <= generators and
                type(ancestor) is int and 0 <= ancestor <= generators and
                sign in (-1, 1), 'invalid substitution')
        require(child not in substitution and ancestor not in substitution and
                child != ancestor, 'substitution is not between distinct current roots')
        if len(word) == 1:
            require(abs(word[0]) == child and ancestor == 0,
                    'length-one relation does not justify the asserted elimination')
        else:
            require(len(word) == 2 and ancestor != 0 and
                    {abs(x) for x in word} == {child, ancestor},
                    'length-two elimination must involve two distinct roots')
        substitution[child] = (ancestor, sign)
        require(not reduced(word), 'asserted signed substitution does not satisfy its relation')
        steps += 1
    remaining = sorted({image(g)[0] for g in range(1, generators + 1)} - {0})
    require(steps == result['eliminations'] and
            remaining == result['remaining_roots'] and
            len(remaining) == result['remaining_generator_classes'],
            'final root accounting differs')
    trivial = not remaining
    require(trivial == result['all_original_generators_proved_identity'],
            'claimed triviality differs')
    kind = None
    boundary = {'components': tri.countBoundaryComponents()}
    if trivial and tri.isOrientable():
        if tri.isClosed():
            kind = 'sphere'
        elif tri.countBoundaryComponents() == 1:
            bc = tri.boundaryComponent(0)
            boundary.update(orientable=bc.isOrientable(), euler_characteristic=bc.eulerChar())
            if bc.isOrientable() and bc.eulerChar() == 2:
                kind = 'ball'
    report = {'status': 'PASS_INDEPENDENT_TRIVIAL_GROUP_PROOF' if trivial
                        else 'PASS_ELIMINATION_PREFIX_INCONCLUSIVE',
              'source_sha256': source_hash, 'presentation_sha256': digest,
              'proof_result_sha256': sha(directory / 'RESULT.json'),
              'elimination_log_sha256': sha(log), 'eliminations': steps,
              'original_generators': generators, 'remaining_generator_classes': len(remaining),
              'triangulation_valid_finite_connected': True,
              'triangulation_orientable': tri.isOrientable(),
              'triangulation_closed': tri.isClosed(), 'boundary': boundary,
              'recognised_link_kind_by_poincare_and_3d_schoenflies': kind,
              'triangle_relators_regenerated_via_tetrahedron_embeddings': True,
              'producer_code_imported': False,
              'regina_group_computation_used': False,
              'geometric_simplification_used': False,
              'standardness_of_any_4_manifold_asserted': False,
              'seconds': time.monotonic() - started, 'checker_sha256': sha(__file__)}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(output)
    print(json.dumps(report), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    verify(args.source, args.directory, args.out)

#!/usr/bin/env python3
"""Prove triviality of a 3D triangulation group using short primal relators.

Only length-one and length-two relators eliminate generators. A torsion
relation such as x*x is retained. Failure to eliminate all generators is
inconclusive. Every elimination is justified by an original triangle under
the substitutions already proved. No geometric simplification is used.
"""
from __future__ import annotations
import argparse
from collections import deque
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


def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


class UnsignedUF:
    def __init__(self, size):
        self.parent = list(range(size))
        self.weight = [1] * size

    def find(self, a):
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        if self.weight[a] > self.weight[b]:
            a, b = b, a
        self.parent[a] = b
        self.weight[b] += self.weight[a]
        return True


class SignedUF:
    """generator[a] = generator[parent[a]] ** sign[a]; root zero is identity."""
    def __init__(self, generators):
        self.parent = list(range(generators + 1))
        self.sign = [1] * (generators + 1)
        self.weight = [1] * (generators + 1)

    def find(self, a):
        parent = self.parent[a]
        if parent == a:
            return a, 1
        root, sign = self.find(parent)
        self.parent[a] = root
        self.sign[a] *= sign
        return root, self.sign[a]

    def reduce(self, raw):
        word = []
        for literal in raw:
            root, sign = self.find(abs(literal))
            if not root:
                continue
            token = root * sign * (1 if literal > 0 else -1)
            if word and word[-1] == -token:
                word.pop()
            else:
                word.append(token)
        while len(word) > 1 and word[0] == -word[-1]:
            word = word[1:-1]
        return word

    def merge_roots(self, child, parent, sign):
        require(child != parent and child != 0 and
                self.parent[child] == child and self.parent[parent] == parent,
                'merge requires two distinct current roots, with nonzero child')
        require(sign in (-1, 1), 'invalid sign')
        self.parent[child] = parent
        self.sign[child] = sign
        self.weight[parent] += self.weight[child]


def presentation(tri):
    require(tri.isValid() and not tri.isIdeal() and tri.isConnected(),
            'requires a valid finite connected 3D triangulation')
    vertices = tri.countVertices()
    tree = UnsignedUF(vertices)
    edge_variables = [0] * tri.countEdges()
    tree_edges = []
    generators = 0
    for edge in tri.edges():
        if tree.join(edge.vertex(0).index(), edge.vertex(1).index()):
            tree_edges.append(edge.index())
        else:
            generators += 1
            edge_variables[edge.index()] = generators
    require(len(tree_edges) == vertices - 1, 'primal spanning tree incomplete')
    # Oriented boundary of a triangle: 0->1, 1->2, 2->0.
    orientation = {(0, 1): (0, 1), (1, 0): (0, -1),
                   (1, 2): (1, 1), (2, 1): (1, -1),
                   (2, 0): (2, 1), (0, 2): (2, -1)}
    relators = []
    digest = hashlib.sha256()
    for triangle in tri.triangles():
        word = [0, 0, 0]
        seen = 0
        for local_edge in range(3):
            mapping = triangle.edgeMapping(local_edge)
            place, sign = orientation[mapping[0], mapping[1]]
            require(not seen & (1 << place), 'duplicate edge of local triangle')
            seen |= 1 << place
            word[place] = sign * edge_variables[triangle.edge(local_edge).index()]
        require(seen == 7, 'incomplete local triangle boundary')
        row = tuple(x for x in word if x)
        relators.append(row)
        digest.update((json.dumps(row, separators=(',', ':')) + '\n').encode())
    info = {'vertices': vertices, 'edges': tri.countEdges(),
            'triangles': tri.countTriangles(), 'tetrahedra': tri.size(),
            'generators': generators, 'tree_edges': tree_edges,
            'presentation_sha256': digest.hexdigest()}
    return relators, info


def solve(source, out, seconds=120):
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina7.4.1')
    source, out = Path(source), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    source_hash = sha(source)
    packet = regina.open(str(source))
    require(isinstance(packet, regina.Triangulation3), 'requires a 3D RGA')
    tri = regina.Triangulation3(packet)
    relators, info = presentation(tri)
    require(sha(source) == source_hash, 'source changed while reading')
    atomic_json(out / 'SOURCE.json', {'source_sha256': source_hash,
                'source_file': str(source.resolve()), **info,
                'triangulation_valid_finite_connected': True,
                'presentation_basis': 'primal1-skeleton spanning tree and oriented triangle boundaries',
                'checker_sha256': sha(__file__)})
    n = info['generators']
    uf = SignedUF(n)
    incidence = [set() for _ in range(n + 1)]
    for index, relator in enumerate(relators):
        for literal in relator:
            incidence[abs(literal)].add(index)
    queue = deque(range(len(relators)))
    queued = bytearray([1]) * len(relators)
    active = bytearray([1]) * len(relators)
    events = []
    remaining_roots = n
    examined = 0
    last_report = started
    while queue and remaining_roots and time.monotonic() - started < seconds:
        index = queue.popleft()
        queued[index] = 0
        if not active[index]:
            continue
        examined += 1
        word = uf.reduce(relators[index])
        if not word:
            active[index] = 0
            continue
        if len(word) == 1:
            child, parent, sign = abs(word[0]), 0, 1
        elif len(word) == 2 and abs(word[0]) != abs(word[1]):
            a, b = word
            child, parent = abs(a), abs(b)
            if uf.weight[child] > uf.weight[parent]:
                child, parent = parent, child
            sign = -1 if (a > 0) == (b > 0) else 1
        else:
            # Includes x*x=1: this is a genuine torsion relation, not x=1.
            continue
        events.append({'step': len(events), 'triangle': index,
                       'reduced_relator': word, 'child': child,
                       'parent': parent, 'sign': sign})
        affected = incidence[child]
        uf.merge_roots(child, parent, sign)
        remaining_roots -= 1
        active[index] = 0
        if parent:
            incidence[parent].update(affected)
        for relation in affected:
            if active[relation] and not queued[relation]:
                queue.append(relation)
                queued[relation] = 1
        incidence[child] = set()
        if time.monotonic() - last_report > 5:
            progress = {'stage': 'SIGNED_UNION_FIND_ELIMINATION',
                        'eliminations': len(events), 'remaining_roots': remaining_roots,
                        'examined_relators': examined, 'queue': len(queue),
                        'seconds': time.monotonic() - started}
            atomic_json(out / 'PROGRESS.json', progress)
            print(json.dumps(progress), flush=True)
            last_report = time.monotonic()
    unresolved = sorted({uf.find(x)[0] for x in range(1, n + 1)} - {0})
    require(len(unresolved) == remaining_roots, 'root accounting differs')
    log = out / 'eliminations.jsonl'
    temporary = log.with_suffix('.tmp')
    temporary.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n'
                                 for row in events))
    temporary.replace(log)
    report = {'status': 'PASS_TRIVIAL_PRIMAL_FUNDAMENTAL_GROUP' if not unresolved
                        else 'INCONCLUSIVE_PRIMAL_ELIMINATION',
              'source_sha256': source_hash, 'presentation_sha256': info['presentation_sha256'],
              'generators': n, 'eliminations': len(events),
              'remaining_generator_classes': len(unresolved),
              'remaining_roots': unresolved, 'elimination_log_sha256': sha(log),
              'all_original_generators_proved_identity': not unresolved,
              'examined_relators': examined, 'seconds': time.monotonic() - started,
              'bound_seconds_including_extraction': seconds,
              'geometric_simplification_used': False,
              'regina_group_computation_used': False,
              'torsion_relations_dropped': False,
              'nontriviality_claimed': False, 'standardness_of_source_claimed': False}
    atomic_json(out / 'RESULT.json', report)
    print(json.dumps({k: v for k, v in report.items() if k != 'remaining_roots'}), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=120)
    args = parser.parse_args()
    solve(args.source, args.out, args.seconds)

#!/usr/bin/env python3
"""Check a closed PL 4-manifold directly from its raw facet gluings.

No Regina 4D skeleton, vertex classification, isValid(), or isIdeal() call is
made.  Face stabilizers are checked by independent union-find structures;
vertex links are explicitly assembled as 3D triangulations.  Each link is
simplified and recognized using Regina 7.4.1.  This is manifold admission,
not a proof of simple connectedness or standard-sphere recognition.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import itertools
import json
from pathlib import Path
import time
import regina


def require(ok, message):
    if not ok:
        raise ValueError(message)


class DSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = bytearray(n)

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def sign(permutation):
    return -1 if sum(permutation[i] > permutation[j]
                     for i in range(len(permutation))
                     for j in range(i + 1, len(permutation))) % 2 else 1


def raw_pairings(tri, closed=True):
    n = tri.size()
    require(n > 0, 'empty input')
    pairs = []
    orientation = [None] * n
    adjacency = [[] for _ in range(n)]
    boundary_facets = 0
    for simplex in tri.simplices():
        i = simplex.index()
        for facet in range(5):
            neighbor = simplex.adjacentSimplex(facet)
            if neighbor is None:
                require(not closed, 'input has a real boundary facet')
                boundary_facets += 1
                continue
            j = neighbor.index()
            p = [simplex.adjacentGluing(facet)[v] for v in range(5)]
            require(sorted(p) == list(range(5)), 'invalid facet permutation')
            target = p[facet]
            require((i, facet) != (j, target), 'facet glued to itself')
            require(neighbor.adjacentSimplex(target).index() == i,
                    'nonreciprocal facet neighbor')
            q = neighbor.adjacentGluing(target)
            require(all(q[p[v]] == v for v in range(5)), 'noninverse facet pairings')
            adjacency[i].append((j, -sign(p)))
            if (i, facet) < (j, target):
                pairs.append((i, facet, j, target, p))
    require(len(pairs) * 2 + boundary_facets == n * 5, 'not all facets accounted for')
    orientation[0] = 1
    stack = [0]
    while stack:
        i = stack.pop()
        for j, factor in adjacency[i]:
            wanted = factor * orientation[i]
            if orientation[j] is None:
                orientation[j] = wanted
                stack.append(j)
            else:
                require(orientation[j] == wanted, 'nonorientable facet gluings')
    require(all(x is not None for x in orientation), 'disconnected dual graph')
    return pairs, orientation


def face_stabilizers(n, pairs):
    """Check that no local positive-dimensional face identifies two of its vertices."""
    counts = {}
    for k in (2, 3, 4):
        templates = list(itertools.combinations(range(5), k))
        indexes = {face: index for index, face in enumerate(templates)}
        width = len(templates)
        faces = DSU(n * width)
        corners = DSU(n * width * k)
        for i, omitted, j, opposite, permutation in pairs:
            for face in itertools.combinations([v for v in range(5) if v != omitted], k):
                mapped = tuple(sorted(permutation[v] for v in face))
                a, b = i * width + indexes[face], j * width + indexes[mapped]
                faces.union(a, b)
                positions = {v: index for index, v in enumerate(mapped)}
                for position, vertex in enumerate(face):
                    corners.union(a * k + position,
                                  b * k + positions[permutation[vertex]])
        for face in range(n * width):
            require(len({corners.find(face * k + position) for position in range(k)}) == k,
                    f'nontrivial self-identification in dimension {k - 1}, local face {face}')
        counts[k - 1] = len({faces.find(face) for face in range(n * width)})
    return counts


def explicit_vertex_links(tri, pairs):
    n = tri.size()
    corners = DSU(5 * n)
    for i, facet, j, target, permutation in pairs:
        for vertex in range(5):
            if vertex != facet:
                corners.union(5 * i + vertex, 5 * j + permutation[vertex])
    groups = {}
    for corner in range(5 * n):
        groups.setdefault(corners.find(corner), []).append(divmod(corner, 5))
    for vertex_index, group in enumerate(sorted(groups.values(), key=lambda g: g[0])):
        link = regina.Triangulation3()
        link.newSimplices(len(group))
        local = {corner: index for index, corner in enumerate(group)}
        orders = [[v for v in range(5) if v != apex] for i, apex in group]
        for source, (i, apex) in enumerate(group):
            pent = tri.simplex(i)
            for facet, old_omitted in enumerate(orders[source]):
                if link.simplex(source).adjacentSimplex(facet) is not None:
                    continue
                neighbor = pent.adjacentSimplex(old_omitted)
                if neighbor is None:
                    continue
                permutation = pent.adjacentGluing(old_omitted)
                target = local[neighbor.index(), permutation[apex]]
                q = [orders[target].index(permutation[v]) for v in orders[source]]
                link.simplex(source).join(facet, link.simplex(target), regina.Perm4(*q))
        yield vertex_index, group, link


def audit(path, out, save_links=False, closed=True):
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina distribution 7.4.1')
    path, out = Path(path), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    source_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    packet = regina.open(str(path))
    require(isinstance(packet, regina.Triangulation4), 'not a 4D input')
    require(hashlib.sha256(path.read_bytes()).hexdigest() == source_sha,
            'source file changed while it was loaded')
    tri = regina.Triangulation4(packet)
    pairs, orientation = raw_pairings(tri, closed=closed)
    counts = face_stabilizers(tri.size(), pairs)
    print(json.dumps({'stage': 'facet_and_face_checks_passed', 'pentachora': tri.size(),
                      'seconds': time.monotonic() - start}), flush=True)
    vertices = []
    for index, embeddings, link in explicit_vertex_links(tri, pairs):
        initial = link.size()
        require(link.isValid() and link.isConnected() and link.isOrientable(),
                f'vertex {index} link is not a valid orientable 3-manifold')
        boundary_vertex = link.hasBoundaryFacets()
        require(not closed or not boundary_vertex, 'closed input has a vertex with boundary link')
        if save_links:
            require(regina.PacketOfTriangulation3(link).save(str(out / f'link_{index:04d}.rga')),
                    'link save failed')
        then = time.monotonic()
        link.simplify()
        recognised = link.isBall() if boundary_vertex else link.isSphere()
        row = {'vertex': index, 'original_tetrahedra': initial,
               'simplified_tetrahedra': link.size(), 'simplified_isosig': link.isoSig(),
               'recognised_kind': 'ball' if boundary_vertex else 'sphere',
               'recognised': bool(recognised),
               'is_sphere': bool(recognised) if not boundary_vertex else False,
               'seconds': time.monotonic() - then,
               'native_3d_simplification_used': True,
               'first_local_corner': list(embeddings[0])}
        vertices.append(row)
        (out / 'VERTEX_LINKS.json').write_text(json.dumps(vertices, indent=2) + '\n')
        print(json.dumps({'link': row}), flush=True)
        require(recognised, f'vertex {index} has a link that is not the required ' +
                ('ball' if boundary_vertex else 'sphere') + ' (nonspherical link for an internal vertex)')
    vector = [len(vertices), counts[1], counts[2], counts[3], tri.size()]
    result = {'status': 'PASS_CLOSED_ORIENTABLE_PL_4_MANIFOLD' if closed else 'PASS_FINITE_ORIENTABLE_PL_4_MANIFOLD',
              'closed_required': closed,
              'real_boundary_facets': sum(s.adjacentSimplex(f) is None for s in tri.simplices() for f in range(5)),
              'source_file': str(path.resolve()), 'source_sha256': source_sha,
              'source_tight_sha256': hashlib.sha256(tri.tightEncoding().encode()).hexdigest(),
              'f_vector': vector, 'euler_characteristic': sum((-1) ** i * n for i, n in enumerate(vector)),
              'all_face_stabilizers_trivial': True, 'facet_pairings_reciprocal': True,
              'simplex_orientations': orientation, 'vertex_links': vertices,
              'four_dimensional_skeleton_used': False,
              'native_3d_recognition_used': True, 'native_4d_recognition_used': False,
              'simply_connectedness_asserted': False, 'standardness_asserted': False,
              'seconds': time.monotonic() - start,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out / ('CLOSED_MANIFOLD_ADMISSION.json' if closed else 'FINITE_MANIFOLD_ADMISSION.json')).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'f_vector', 'euler_characteristic', 'seconds')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--save-links', action='store_true')
    parser.add_argument('--finite', action='store_true', help='allow real boundary and recognise its vertex links as balls')
    args = parser.parse_args()
    audit(args.source, args.out, args.save_links, closed=not args.finite)

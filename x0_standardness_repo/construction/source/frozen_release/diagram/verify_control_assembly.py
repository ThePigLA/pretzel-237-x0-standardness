#!/usr/bin/env python3
"""Independent stored-data check of the complete marked control assembly.

No producer, including assemble_control.py, is imported.  Regenerates every
product pentachoron from the ordered input pieces and checks native covering
and orbit-quotient gluing squares, the full fixed subcomplex, both marked
surfaces, source seam carriers and the source-oriented product anchor.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip
import json
import time

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def req(x, why):
    if not x:
        raise ValueError(why)


def load(name):
    p = ROOT / name
    b = p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix == '.gz' else b)


def digest(x):
    return sha256(json.dumps(x, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def filehash(name):
    return sha256((ROOT / name).read_bytes()).hexdigest()


def tup(x):
    return tuple(map(tup, x)) if isinstance(x, list) else x


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def determinant(rows):
    a = [list(map(Q, row)) for row in rows]
    out = Q(1)
    for k in range(len(a)):
        pivot_row = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot_row is None:
            return Q(0)
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            out = -out
        pivot = a[k][k]
        out *= pivot
        for i in range(k + 1, len(a)):
            scale = a[i][k] / pivot
            for j in range(k + 1, len(a)):
                a[i][j] -= scale * a[k][j]
    return out


def product_cells(base, a, b):
    # The sole vertical step can follow any base vertex.  In particular
    # the positively oriented last circle edge is (3,0), never (0,3).
    for k in range(len(base)):
        yield [(v, a) for v in base[:k + 1]] + [(v, b) for v in base[k:]]


def boundary(fs):
    counts = Counter(face for f in fs for face in combinations(f, len(f) - 1))
    req(all(n in (1, 2) for n in counts.values()), 'Invalid facet incidence in a piece')
    return {face for face, n in counts.items() if n == 1}


def closure(fs):
    ans = set()
    for f in fs:
        for n in range(1, len(f) + 1):
            ans.update(combinations(f, n))
    return ans


def connected(adj):
    req(adj, 'Empty adjacency graph')
    seen = {next(iter(adj))}
    queue = deque(seen)
    while queue:
        v = queue.popleft()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return len(seen) == len(adj)


def surface(triangles):
    triangles = set(triangles)
    verts = {v for f in triangles for v in f}
    edge_count = Counter(e for f in triangles for e in combinations(f, 2))
    req(set(edge_count.values()) == {2}, 'Surface has noncircular edge links')
    vertex_links = {v: defaultdict(set) for v in verts}
    graph = defaultdict(set)
    for a, b in edge_count:
        graph[a].add(b)
        graph[b].add(a)
    req(connected(graph), 'Surface disconnected')
    for f in triangles:
        for v in f:
            a, b = [w for w in f if w != v]
            vertex_links[v][a].add(b)
            vertex_links[v][b].add(a)
    for link in vertex_links.values():
        req(all(len(n) == 2 for n in link.values()) and connected(link), 'Surface vertex link not S1')
    fvec = [len(verts), len(edge_count), len(triangles)]
    req(fvec[0] - fvec[1] + fvec[2] == 2, 'Marked surface is not a sphere')
    return fvec


def half_facet_maps(data):
    n = data['simplex_count']
    maps = {}
    for i, k, j, p in data['gluings']:
        req(0 <= i < n and 0 <= j < n and 0 <= k < 5 and sorted(p) == list(range(5)), 'Invalid native gluing row')
        l = p[k]
        req((i, k) != (j, l) and (i, k) not in maps and (j, l) not in maps, 'Repeated/self half-facet gluing')
        maps[i, k] = (j, tuple(p))
        maps[j, l] = (i, tuple(p.index(v) for v in range(5)))
    req(len(maps) == 5 * n, 'Native gluing has an unpaired face')
    return maps


class UnionFind:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, v):
        while self.p[v] != v:
            self.p[v] = self.p[self.p[v]]
            v = self.p[v]
        return v

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.p[max(a, b)] = min(a, b)


def verify():
    start = time.monotonic()
    m = load('control/marked_control.json.gz')
    cov, quo = load('control/cover_regina.json.gz'), load('control/quotient_regina.json.gz')
    adaptor = load('adaptor/model.json.gz')
    central = load('central_roll/model.json.gz')
    central_cert = load('central_roll/certificate.json.gz')
    bridge, cap = load('cap/data/source_cap_bridge.json.gz'), load('cap/data/cap_m2.json.gz')
    cut, cut_ori = load('peripheral/cut_source.json.gz'), load('peripheral/cut_orientation.json')
    bindings = {'adaptor': digest(adaptor), 'central': digest(central), 'bridge': digest(bridge), 'cap': digest(cap)}
    req(m['input_bindings'] == bindings, 'Assembly input bindings differ')
    req(central['central_certificate_sha256'] == digest(central_cert), 'Central certificate binding')
    central_gate = load('central_roll/independent_verification.json')
    req(central_gate['passed'] and central_gate['model_sha256'] == digest(central)
        and central_gate['certificate_sha256'] == digest(central_cert), 'Independent central geometry gate stale')
    req(central_gate['straightening']['full_integer_winding'] == 1,
        'Independent central gate did not certify full roll winding')
    req(m['twist'] == 2 and m['roll'] == 1 and m['grid'] == [3, 4, 4], 'Wrong control parameters')
    req(cap['twist'] == 2 and cap['beta']['cap_full_turns'] == -1, 'Full m2 cap winding absent')
    fs = list(map(tuple, m['facets']))
    req(fs == sorted(set(fs)) and all(len(f) == 5 and tuple(sorted(set(f))) == f for f in fs), 'Folded or repeated pentachora')
    n = len(fs)
    req(n <= m['construction_limits']['max_pentachora'], 'Declared construction cap exceeded')
    req(digest(fs) == m['facets_sha256'] == cov['source_facets_sha256'] == quo['source_facets_sha256'], 'Native source facets differ')
    tags = list(map(tup, m['vertices']))
    idx = {t: i for i, t in enumerate(tags)}
    req(len(idx) == len(tags), 'Repeated global vertex tag')
    deck = m['deck']
    req(len(deck) == len(tags) and sorted(deck) == list(range(len(tags))) and all(deck[deck[v]] == v for v in range(len(tags))), 'Invalid deck permutation')
    f_index = {f: i for i, f in enumerate(fs)}
    ids = lambda vertices: tuple(sorted(idx[v] for v in vertices))
    parts = {}
    deck_expected = {}

    def assign(v, w):
        req(v not in deck_expected or deck_expected[v] == w, 'Conflicting deck on attachment')
        deck_expected[v] = w

    def prod(v, j):
        return ('B', v, j % 4)

    stages = {s: list(map(tup, ts)) for s, ts in adaptor['stages'].items()}
    d3 = {tup(a): tup(b) for a, b in adaptor['deck']}
    for side in ('N', 'E'):
        part = []
        for tet in stages[side]:
            for j in range(4):
                part.extend(ids(prod(v, t) for v, t in cell) for cell in product_cells(tet, j, (j + 1) % 4))
        for v in {v for tet in stages[side] for v in tet}:
            for j in range(4):
                assign(prod(v, j), prod(d3[v], j))
        parts[side] = part

    def cent(v, j):
        if v[0] in ('a', 'b'):
            return prod(('N' if v[0] == 'a' else 'E', 'g', v[1], j % 4), v[2])
        return ('C', v, j % 4)

    central_tets = list(map(tup, central['tetrahedra']))
    parts['central'] = [ids(cent(v, t) for v, t in cell) for tet in central_tets for j in range(4)
                        for cell in product_cells(tet, j, (j + 1) % 4)]
    for v in {v for tet in central_tets for v in tet}:
        for j in range(4):
            assign(cent(v, j), cent(v, j + 2))

    source_port = bridge['ports']['source']
    cap_port = bridge['ports']['cap']
    req(len(source_port['vertices']) == len(bridge['source_boundary']['vertices']), 'Source port size')
    smap = {v: prod(('N', 'o', a), j) for v, (a, j) in zip(source_port['vertices'], bridge['source_boundary']['vertices'])}
    cmap = {v: ('K', j) for j, v in enumerate(cap_port['vertices'])}
    req(not (set(smap) & set(cmap)), 'Bridge ports intersect')
    bmap = {v: smap.get(v, cmap.get(v, ('H', v))) for v in range(len(bridge['vertices']))}
    parts['bridge'] = [ids(bmap[v] for v in f) for f in bridge['facets']]
    parts['cap'] = [ids(('K', v) for v in f) for f in cap['facets']]
    for v, w in enumerate(bridge['deck']):
        assign(bmap[v], bmap[w])
    for v, w in enumerate(cap['deck']):
        assign(('K', v), ('K', w))
    req(set(deck_expected) == set(tags), 'Unused/missing assembly vertices')
    req([idx[deck_expected[v]] for v in tags] == deck, 'Saved deck differs from source pieces')
    regenerated = [f for pp in parts.values() for f in pp]
    req(sorted(regenerated) == fs, 'Full facet regeneration differs')
    for name, pp in parts.items():
        req(sorted(pp) == list(map(tuple, m['parts'][name])), 'Saved part differs: ' + name)
    print('Regenerated every product, bridge and cap pentachoron', flush=True)

    bulk = parts['N'] + parts['E'] + parts['central']
    bulk_boundary = boundary(bulk)
    expected_source_port = {ids(smap[v] for v in f) for f in source_port['facets']}
    iface = bridge['source_interface']
    literal_source_product = {ids(prod(('N', 'o', v), t) for v, t in cell)
                              for tri in iface['ordered_triangles'] for j in range(4)
                              for cell in product_cells(tri, j, (j + 1) % 4)}
    req(bulk_boundary == expected_source_port == literal_source_product == set(map(tuple, m['bulk_boundary'])), 'Bulk puncture/source bridge port mismatch')
    expected_cap_port = {ids(cmap[v] for v in f) for f in cap_port['facets']}
    req(expected_cap_port == {ids(('K', v) for v in f) for f in cap['boundary']}, 'Cap port mismatch')
    req(m['source_port_map'] == [[v, idx[smap[v]]] for v in source_port['vertices']], 'Saved source map differs')
    req(m['cap_port_map'] == [[v, idx[('K', v)]] for v in range(len(cap_port['vertices']))], 'Saved cap map differs')

    # Cover face maps and orientation are checked against literal abstract
    # simplex faces.  The incidence count separately rules out extra sheets.
    req(cov['simplex_count'] == n, 'Cover simplex count')
    cfaces = half_facet_maps(cov)
    face_counts = Counter(face for f in fs for face in combinations(f, 4))
    req(set(face_counts.values()) == {2}, 'Assembled cover is not closed with degree-two facets')
    signs = m['orientation']
    req(len(signs) == n and set(signs) <= {-1, 1}, 'Invalid orientation chain')
    graph = {i: [] for i in range(n)}
    for i, k, j, p in cov['gluings']:
        l = p[k]
        req(all(fs[i][a] == fs[j][p[a]] for a in range(5) if a != k), 'Cover gluing changes face vertex')
        req(fs[i][:k] + fs[i][k + 1:] == fs[j][:l] + fs[j][l + 1:], 'Cover gluing mismatches abstract face')
        req(signs[j] == -signs[i] * parity(p), 'Orientation does not cancel on a face')
        graph[i].append(j)
        graph[j].append(i)
    req(connected(graph), 'Assembled cover disconnected')
    cell_deck = []
    for i, f in enumerate(fs):
        req(len({min(v, deck[v]) for v in f}) == 5, 'Simplex has two vertices in the same deck orbit')
        image = tuple(deck[v] for v in f)
        j = f_index.get(tuple(sorted(image)))
        req(j is not None and j != i, 'Deck top-cell action failed')
        req(signs[j] == signs[i] * parity(image), 'Deck reverses the anchored ambient orientation')
        cell_deck.append(j)

    # The global source anchor uses a genuinely unchanged N tetrahedron.
    anchor = m['orientation_anchor']
    req(anchor['cut_orientation_sha256'] == digest(cut_ori), 'Orientation anchor cut binding')
    k = anchor['source_N_tetrahedron_index']
    original = anchor['source_cut_vertices_in_product_order']
    sorted_original = tuple(sorted(original))
    req(tuple(cut['N'][k]) == sorted_original, 'Anchor original tetrahedron differs')
    base = tuple(('N', 'o', v) for v in original)
    req(base in stages['N'], 'Anchor is not an unchanged ordered source tetrahedron')
    source_sign = cut_ori['N_tetrahedron_orientations'][k]
    base_sign = source_sign * parity(original)
    ordered = [idx[prod(v, 0)] for v in base] + [idx[prod(base[-1], 1)]]
    facet = tuple(sorted(ordered))
    canonical_sign = base_sign * parity(ordered)
    req(anchor['source_sorted_chain_sign'] == source_sign and anchor['branch_order_sign'] == base_sign,
        'Anchor source/branch orientation differs')
    req(anchor['ordered_product_vertices'] == ordered and tuple(anchor['canonical_pentachoron']) == facet,
        'Anchor product vertices differ')
    req(anchor['time_edge'] == [0, 1] and anchor['pentachoron_index'] == f_index[facet], 'Anchor time/index differs')
    req(anchor['canonical_orientation'] == canonical_sign == signs[f_index[facet]], 'Global orientation not source-anchored')
    standard_base = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    coordinates = {}
    for v in facet:
        tag = tags[v]
        coordinates[v] = (*standard_base[sorted_original.index(tag[1][2])], Q(tag[2], 4))
    points = [coordinates[v] for v in facet]
    actual_determinant = determinant([[b - a for a, b in zip(points[0], point)] for point in points[1:]])
    req(abs(actual_determinant) == Q(1, 4) and source_sign * (1 if actual_determinant > 0 else -1) == canonical_sign,
        'Canonical anchor affine determinant disagrees with source orientation')
    # Expressing the same simplex in the given branch order gives the
    # last-vertical product (0,0),(e1,0),(e2,0),(e3,0),(e3,1/4).
    local_points = [(*p, Q(0)) for p in standard_base] + [(*standard_base[-1], Q(1, 4))]
    local_determinant = determinant([[b - a for a, b in zip(local_points[0], point)] for point in local_points[1:]])
    req(local_determinant == Q(1, 4), 'Last-vertical staircase determinant failed')
    anchor_proof = {'source_N_tetrahedron_index': k, 'source_sorted_chain_sign': source_sign,
                    'source_vertices_in_product_order': original, 'branch_order_sign': base_sign,
                    'last_vertical_product_affine_determinant': str(local_determinant),
                    'canonical_vertices_affine_determinant_in_sorted_source_then_time': str(actual_determinant),
                    'canonical_pentachoron_index': f_index[facet], 'canonical_orientation': canonical_sign,
                    'meaning': 'canonical Seifert master = signed-diagram ambient, followed by positive time'}

    # Full fixed locus, including isolated vertices or edges that a scan of
    # fixed triangles alone would miss.  The orbit condition above rules
    # out inversions, so all fixed points lie in these fixed-vertex faces.
    fixed_vertices = {v for v, w in enumerate(deck) if v == w}
    fixed_faces = set()
    normal_edges = defaultdict(set)
    for f in fs:
        fixed = tuple(v for v in f if v in fixed_vertices)
        req(len(fixed) <= 3, 'Fixed locus has dimension above two')
        fixed_faces.update(closure([fixed]))
        if len(fixed) == 3:
            normal_edges[fixed].add(tuple(v for v in f if v not in fixed_vertices))
    ram = set(map(tuple, m['ramification_triangles']))
    req(fixed_faces == closure(ram), 'Full fixed subcomplex contains omitted lower-dimensional components')
    rf = surface(ram)
    for tri, es in normal_edges.items():
        link = defaultdict(set)
        for a, b in es:
            link[a].add(b)
            link[b].add(a)
        req(all(len(z) == 2 for z in link.values()) and connected(link), 'Fixed triangle normal link is not S1')
        req(all(deck[v] != v and deck[v] in link for v in link), 'Normal-link deck has a fixed vertex')
        req({tuple(sorted((deck[a], deck[b]))) for a, b in es} == es, 'Normal-link deck changes an edge')
        req(all(tuple(sorted((deck[a], deck[b]))) != (a, b) for a, b in es), 'Normal-link deck fixes an edge')
    req(set(normal_edges) == ram, 'Missing fixed-triangle normal link')
    twin = set(map(tuple, m['twin_triangles']))
    req(twin == {ids(('K', v) for v in f) for f in cap['twin_triangles']}, 'Twin changed under cap embedding')
    tf = surface(twin)
    all_triangles = {face for f in fs for face in combinations(f, 3)}
    req(twin <= all_triangles and ram <= all_triangles, 'Marked surface missing from complex')
    intersection = closure(twin) & fixed_faces
    req(len(intersection) == 2 and all(len(f) == 1 for f in intersection), 'Twin/fixed intersection not exactly two vertices')
    intersection_proof = []
    for (v,) in sorted(intersection):
        tag = tags[v]
        req(tag[0] == 'K', 'Twin intersection outside cap')
        local = tag[1]
        req({f for f in fs if v in f} == {ids(('K', w) for w in f) for f in cap['facets'] if local in f},
            'Twin intersection star changed outside cap')
        intersection_proof.append({'global_vertex': v, 'cap_vertex': local, 'full_ambient_star_unchanged_from_cap': True})

    # Exact orbit-cell quotient, not a quotient formed by merging vertices.
    reps = [i for i, j in enumerate(cell_deck) if i < j]
    rep_index = {r: k for k, r in enumerate(reps)}
    req(quo['cell_deck'] == cell_deck and quo['cover_cell_representatives'] == reps, 'Quotient cell orbits differ')
    req(quo['simplex_count'] == len(reps) and 2 * len(reps) == n, 'Quotient size')
    qfaces = half_facet_maps(quo)
    covering = []
    for i, f in enumerate(fs):
        r = min(i, cell_deck[i])
        p = [fs[r].index(deck[v] if i != r else v) for v in f]
        covering.append([rep_index[r], p])
    req(quo['covering_simplex_maps'] == covering, 'Covering simplex maps differ')
    orbit_labels = [min(v, deck[v]) for v in range(len(deck))]
    req(quo['vertex_orbit_labels'] == orbit_labels, 'Quotient vertex orbit labels differ')
    for i, k, j, p in cov['gluings']:
        qi, a = covering[i]
        qj, b = covering[j]
        dest, qp = qfaces[qi, a[k]]
        req(dest == qj and [qp[a[v]] for v in range(5)] == [b[p[v]] for v in range(5)],
            'Orbit quotient face square does not commute')
    uf = UnionFind(5 * len(reps))
    for i, k, j, p in quo['gluings']:
        req(signs[reps[j]] == -signs[reps[i]] * parity(p), 'Quotient orientation does not descend')
        for v in range(5):
            if v != k:
                uf.join(5 * i + v, 5 * j + p[v])
    class_label, label_class = {}, {}
    for i, r in enumerate(reps):
        for v, source_v in enumerate(fs[r]):
            cl, lab = uf.find(5 * i + v), orbit_labels[source_v]
            req(cl not in class_label or class_label[cl] == lab, 'Quotient gluing merges distinct vertex orbits')
            req(lab not in label_class or label_class[lab] == cl, 'Quotient splits a vertex orbit')
            class_label[cl], label_class[lab] = lab, cl
    req(set(label_class) == set(orbit_labels), 'Quotient omits a vertex orbit')
    print('Full fixed subcomplex and every orbit-quotient gluing square verified', flush=True)

    seam_proofs = []
    for row in central_cert['source_seam_witnesses']:
        req(Q(row['full_unwrapped_source_winding']) == Q(row['input_s']), 'Source winding witness changed')
        for cc in row['carriers']:
            tagged = list(map(tup, cc['expected_product_carrier']))
            ordered = [idx[v] for v in tagged]
            f = tuple(sorted(ordered))
            req(f in f_index and f in set(map(tuple, m['parts'][cc['part']])), 'Source seam carrier is absent from actual part')
            weights = list(map(Q, cc['product_barycentric_weights']))
            req(len(weights) == 5 and min(weights) >= 0 and sum(weights) == 1, 'Source seam barycentric point invalid')
            if cc['part'] in ('N', 'E'):
                base = tup(cc['base_tetrahedron'])
                a, b = (3, 0) if cc['time'] == 1 else (0, 1)
                possible = [list(prod(v, t) for v, t in cell) for cell in product_cells(base, a, b)]
                req(tagged in possible, 'Source seam carrier ignored cyclic edge branching')
                if cc['time'] == 1:
                    req(weights[0] == 0 and tagged[0][-1] == 3 and all(z[-1] == 0 for z in tagged[1:]),
                        'Time-one seam carrier does not use the outgoing 3->0 edge')
            seam_proofs.append({'input_s': row['input_s'], 'time': cc['time'], 'part': cc['part'],
                                'pentachoron_index': f_index[f], 'canonical_vertices': f,
                                'canonical_barycentric_weights': [str(weights[ordered.index(v)]) for v in f]})
    req(len(seam_proofs) == 10, 'Incomplete source seam witness list')

    fvec = [len(tags), len({e for f in fs for e in combinations(f, 2)}), len(all_triangles), len(face_counts), n]
    req(sum((-1) ** i * z for i, z in enumerate(fvec)) == 2, 'Cover Euler characteristic differs from two')
    native = {}
    for which, data in [('cover', cov), ('quotient', quo)]:
        receipt = load(f'control/{which}_native_admission.json')
        req(receipt['source_sha256'] == filehash(f'control/{which}_regina.json.gz'), 'Native receipt source changed')
        req(all(receipt[k] for k in ('passed', 'valid', 'closed', 'connected', 'orientable')) and not receipt['ideal'],
            'Native closed-manifold admission failed')
        req(receipt['simplex_count'] == data['simplex_count'] and receipt['chi'] == 2, 'Native dimensions mismatch')
        native[which] = {'bound_to_current_file': True, 'simplex_count': data['simplex_count'], 'closed_nonideal_valid': True}
    req(load('control/cover_native_admission.json')['f_vector'] == fvec, 'Independent/native cover f-vectors disagree')
    source_files = ['control/marked_control.json.gz', 'control/cover_regina.json.gz', 'control/quotient_regina.json.gz',
                    'control/cover_native_admission.json', 'control/quotient_native_admission.json',
                    'adaptor/model.json.gz', 'central_roll/model.json.gz', 'central_roll/certificate.json.gz',
                    'central_roll/independent_verification.json',
                    'cap/data/source_cap_bridge.json.gz', 'cap/data/cap_m2.json.gz',
                    'peripheral/cut_orientation.json', 'diagram/verify_control_assembly.py']
    return {'schema': 'independent-marked-control-assembly-verification-v1', 'passed': True,
            'source_model_canonical_sha256': digest(m), 'source_facets_sha256': digest(fs),
            'parts_regenerated': {k: len(v) for k, v in parts.items()}, 'f_vector': fvec,
            'source_port_tetrahedra': len(bulk_boundary), 'cap_port_tetrahedra': len(expected_cap_port),
            'cover_face_pairings': len(cov['gluings']), 'quotient_face_pairings': len(quo['gluings']),
            'quotient_cells': len(reps), 'quotient_vertices': len(label_class),
            'all_covering_face_squares_commute': True, 'deck_orientation_preserving': True,
            'no_simplex_inversions': True, 'full_fixed_subcomplex_checked': True,
            'omitted_fixed_vertices_or_edges': 0, 'ramification_f_vector': rf, 'twin_f_vector': tf,
            'normal_branch_circles_checked': len(normal_edges), 'twin_intersection': intersection_proof,
            'orientation_anchor': anchor_proof, 'source_seam_carriers': seam_proofs,
            'full_source_roll_winding': '1', 'quarter_source_displacement': ['1/4', '0'],
            'bound_native_admission': native, 'input_bindings': bindings,
            'file_sha256': {name: filehash(name) for name in source_files},
            'scope': 'Actual assembled finite cells, source orientation and markings, and orbit quotient. Global source-twist localization and standardization are separate gates.',
            'seconds': time.monotonic() - start}


if __name__ == '__main__':
    result = verify()
    (HERE / 'control_assembly_verification.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('file_sha256', 'source_seam_carriers')}, indent=2))

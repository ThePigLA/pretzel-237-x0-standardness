#!/usr/bin/env python3
"""Exact marked-cusp certificate; verification needs only Python's stdlib.

The input .tri files are read literally, without a SnapPy import that could
renormalize their peripheral curves.  Isomorphisms are enumerated by choosing
the image of tetrahedron 0 and its four vertices and then propagating over
every face.  A peripheral comparison is certified by an integer potential
on the vertices of the actual triangulated cusp, not a numerical isometry.

Run with --write-certificate after regenerating the source files.  Without
that flag, the full deterministic certificate is recomputed and compared.
"""

from collections import deque
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import argparse
import gzip
import json


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent

# Right-sheet corners, anticlockwise as viewed from the cusp.  These are the
# other endpoints of t3m.TruncatedSimplexCorners at the selected tet vertex.
CORNERS = ((1, 2, 3), (0, 3, 2), (3, 0, 1), (2, 1, 0))
PINNED_COMMIT = "42217b3867ce3a00d1125a1b6860db335684c333"
PRIMARY_SOURCES = {
    "kernel/kernel_code/peripheral_curves.c":
        "2725d13a8ac5030c707b21a181212f52f6accd10a8bbd55c5e3ffae5972a0a2c",
    "kernel/kernel_code/intersection_numbers.c":
        "7b67cc7122372acb5db318c1570bad5568eae443aceb4b3734e3678dcf3ca972",
    "kernel/kernel_code/orient.c":
        "9a70ee73ae6c8ae5b8c01d925c2a106ec51bc6d5691a8bde4d87a4f1f0c61374",
    "kernel/headers/kernel_typedefs.h":
        "80017c52afb8dab93fd48a9585518c808c10966d1a4cfecf8b11ec3d5cfb1590",
    "cython/core/triangulation.pyx":
        "eacdf67d5a09e786e92e70d884bf800ea23fffe9d4d3b87c910d1d8eeba5fae9",
}


def require(test, message):
    if not test:
        raise ValueError(message)


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def canonical_hash(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_json(path):
    data = Path(path).read_bytes()
    return json.loads(gzip.decompress(data) if Path(path).suffix == ".gz" else data)


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))


def inverse(p):
    return [p.index(i) for i in range(4)]


def read_tri(path):
    """Read the oriented, unfilled torus-cusp .tri format used in this bundle."""
    lines = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    require(lines[0] == "% Triangulation", "Not a SnapPea triangulation")
    require(lines[3] == "oriented_manifold", "Expected explicitly oriented input")
    cusp_counts = list(map(int, lines[5].split()))
    require(cusp_counts == [1, 0], "This certificate expects exactly one torus cusp")
    cusp = lines[6].split()
    require(cusp[0] == "torus" and all(Q(z) == 0 for z in cusp[1:]), "Unexpected filling")
    n = int(lines[7])
    require(len(lines) == 8 + 8 * n, "Unexpected .tri record length")
    data, cusp_indices, curves = [], [], []
    for t in range(n):
        block = lines[8 + 8 * t:16 + 8 * t]
        neighbors = list(map(int, block[0].split()))
        glue = [list(map(int, s)) for s in block[1].split()]
        ci = list(map(int, block[2].split()))
        rows = [list(map(int, line.split())) for line in block[3:7]]
        require(len(neighbors) == len(glue) == len(ci) == 4, "Bad tet record")
        require(all(sorted(p) == [0, 1, 2, 3] for p in glue), "Bad permutation")
        require(all(len(r) == 16 for r in rows), "Bad peripheral row")
        require(all(0 <= j < n for j in neighbors), "Unpaired face")
        require(all(z in (-1, 0) for z in ci), "Unexpected cusp index")
        data.append([neighbors, glue])
        cusp_indices.append(ci)
        curves.extend(rows)
    result = {"data": data, "cusp_indices": cusp_indices, "curves": curves}
    check_gluings(result)
    return result


def check_gluings(tri):
    for t, (neighbors, glue) in enumerate(tri["data"]):
        for f, j in enumerate(neighbors):
            p = glue[f]
            require(parity(p) == -1, "Input tetrahedra do not carry a consistent positive orientation")
            require(tri["data"][j][0][p[f]] == t, "Face pairing not reciprocal")
            require(tri["data"][j][1][p[f]] == inverse(p), "Face maps not inverse")
            for v in range(4):
                if v != f:
                    require(tri["cusp_indices"][t][v] == tri["cusp_indices"][j][p[v]], "Cusp index not glued")


class UnionFind:
    def __init__(self, values):
        self.parent = {v: v for v in values}

    def find(self, v):
        p = self.parent[v]
        if p != v:
            self.parent[v] = self.find(p)
        return self.parent[v]

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def side_corners(v, f):
    cs = CORNERS[v]
    k = cs.index(f)
    return cs[(k + 1) % 3], cs[(k + 2) % 3]


def cusp_complex(tri):
    triangles = [(t, v) for t, row in enumerate(tri["cusp_indices"])
                 for v, cusp in enumerate(row) if cusp == 0]
    germs = [(t, v, w) for t, v in triangles for w in range(4) if w != v]
    uf = UnionFind(germs)
    for t, v, w in germs:
        for f in range(4):
            if f not in (v, w):
                j, p = tri["data"][t][0][f], tri["data"][t][1][f]
                uf.join((t, v, w), (j, p[v], p[w]))
    groups = {}
    for g in germs:
        groups.setdefault(uf.find(g), []).append(g)
    classes = [groups[k] for k in sorted(groups)]
    vertex = {g: i for i, group in enumerate(classes) for g in group}
    edges = []
    side_to_edge = {}
    for t, v in triangles:
        for f in range(4):
            if f == v:
                continue
            side = (t, v, f)
            j, p = tri["data"][t][0][f], tri["data"][t][1][f]
            other = (j, p[v], p[f])
            require(side != other, "Self-paired cusp side")
            if side > other:
                continue
            a, b = side_corners(v, f)
            aa, bb = side_corners(p[v], p[f])
            va, vb = vertex[(t, v, a)], vertex[(t, v, b)]
            require(va == vertex[(j, p[v], bb)] and vb == vertex[(j, p[v], aa)], "Cusp edge orientations do not oppose")
            index = len(edges)
            edges.append({"side": list(side), "opposite_side": list(other), "vertices": [va, vb]})
            side_to_edge[side] = (index, 1)
            side_to_edge[other] = (index, -1)
    require(len(classes) - len(edges) + len(triangles) == 0, "Cusp Euler characteristic is not zero")
    require(len(side_to_edge) == 3 * len(triangles), "Unpaired cusp edge")
    # The tree computation also proves connectedness of the cusp.
    obj = {"vertices": classes, "edges": edges, "triangles": triangles, "side_to_edge": side_to_edge}
    tree_potential(obj, [0] * len(edges))
    return obj


def get_flows(tri, curve):
    # Only the right orientation-cover sheet is populated in these oriented
    # inputs.  Do not silently add a possibly inconsistent left-sheet curve.
    flows = {}
    for t, ci in enumerate(tri["cusp_indices"]):
        row = tri["curves"][4 * t + 2 * curve]
        left = tri["curves"][4 * t + 2 * curve + 1]
        require(not any(left), "Unexpected nonzero left-sheet curve")
        for v in range(4):
            require(row[4 * v + v] == 0, "Diagonal curve entry must vanish")
            for f in range(4):
                if ci[v] < 0 or v == f:
                    require(row[4 * v + f] == 0, "Curve at finite vertex")
                else:
                    flows[(t, v, f)] = row[4 * v + f]
    return flows


def cochain(complex_, flows):
    # Kernel F > 0 enters the triangle.  On the right-sheet oriented side
    # a->b the corresponding primal cochain is -F.
    result = []
    for e in complex_["edges"]:
        s, o = tuple(e["side"]), tuple(e["opposite_side"])
        require(flows[s] == -flows[o], "Normal flow does not match across edge")
        result.append(-flows[s])
    for t, v in complex_["triangles"]:
        require(sum(flows[(t, v, f)] for f in range(4) if f != v) == 0,
                "Normal flow is not closed in triangle")
    return result


def tree_potential(complex_, values):
    n = len(complex_["vertices"])
    adjacent = [[] for _ in range(n)]
    for i, e in enumerate(complex_["edges"]):
        a, b = e["vertices"]
        adjacent[a].append((b, values[i]))
        adjacent[b].append((a, -values[i]))
    pot = [None] * n
    pot[0] = 0
    todo = deque([0])
    while todo:
        a = todo.popleft()
        for b, value in adjacent[a]:
            if pot[b] is None:
                pot[b] = pot[a] + value
                todo.append(b)
    require(all(x is not None for x in pot), "Disconnected cusp")
    return pot


def residual(complex_, values):
    h = tree_potential(complex_, values)
    return [z - h[e["vertices"][1]] + h[e["vertices"][0]]
            for z, e in zip(values, complex_["edges"])]


def wedge_integral(complex_, first, second):
    """Exact integral of affine representatives; total is the cup pairing."""
    def val(c, t, v, f):
        i, sign = complex_["side_to_edge"][(t, v, f)]
        return sign * c[i]
    total = Q(0)
    for t, v in complex_["triangles"]:
        a, b, c = CORNERS[v]
        # side missing c is a->b; side missing a is b->c.
        total += Q(val(first, t, v, c) * val(second, t, v, a)
                   - val(second, t, v, c) * val(first, t, v, a), 2)
    require(total.denominator == 1, "Cusp pairing is not integral")
    return int(total)


def exact_coordinates(complex_, value, meridian, longitude):
    rm, rl, rv = (residual(complex_, c) for c in (meridian, longitude, value))
    pair = None
    for i in range(len(rm)):
        for j in range(i + 1, len(rm)):
            det = rm[i] * rl[j] - rm[j] * rl[i]
            if det:
                pair = (i, j, det)
                break
        if pair is not None:
            break
    require(pair is not None, "Target peripheral cochains do not span H1")
    i, j, det = pair
    a = Q(rv[i] * rl[j] - rv[j] * rl[i], det)
    b = Q(rm[i] * rv[j] - rm[j] * rv[i], det)
    require(a.denominator == b.denominator == 1, "Nonintegral cusp map")
    a, b = int(a), int(b)
    delta = [z - a * m - b * l for z, m, l in zip(value, meridian, longitude)]
    h = tree_potential(complex_, delta)
    for e, z in zip(complex_["edges"], delta):
        va, vb = e["vertices"]
        require(h[vb] - h[va] == z, "Peripheral difference is not an exact integer cochain")
    return {"target_M_L_coefficients": [a, b], "integer_vertex_potential": h,
            "basis_detection_edges": [i, j], "basis_detection_determinant": det}


def verify_iso(source, target, image, perms):
    n = len(source["data"])
    require(len(target["data"]) == n and sorted(image) == list(range(n)), "Tet images not a bijection")
    require(all(sorted(p) == list(range(4)) for p in perms), "Invalid vertex map")
    signs = [parity(p) for p in perms]
    require(len(set(signs)) == 1, "Isomorphism has inconsistent orientation")
    for t, (neighbors, glue) in enumerate(source["data"]):
        tt, p = image[t], perms[t]
        for v in range(4):
            require(source["cusp_indices"][t][v] == target["cusp_indices"][tt][p[v]], "Isomorphism changes cusp label")
        for f, j in enumerate(neighbors):
            jj = target["data"][tt][0][p[f]]
            q = target["data"][tt][1][p[f]]
            require(jj == image[j], "Isomorphism changes adjacency")
            require([q[p[v]] for v in range(4)] == [perms[j][glue[f][v]] for v in range(4)], "Isomorphism gluing square fails")
    return signs[0]


def enumerate_isos(source, target):
    """At most 24*n seeds; all connected triangulation isomorphisms."""
    n = len(source["data"])
    require(len(target["data"]) == n, "Different tetrahedron counts")
    answers = []
    for first_tet in range(n):
        for first_perm in permutations(range(4)):
            images, perms = {0: first_tet}, {0: list(first_perm)}
            used = {first_tet: 0}
            todo = deque([0])
            good = True
            while todo and good:
                t = todo.popleft()
                tt, p = images[t], perms[t]
                for f, j in enumerate(source["data"][t][0]):
                    s = source["data"][t][1][f]
                    invs = inverse(s)
                    jj = target["data"][tt][0][p[f]]
                    q = target["data"][tt][1][p[f]]
                    pj = [q[p[invs[w]]] for w in range(4)]
                    if j in images:
                        if images[j] != jj or perms[j] != pj:
                            good = False
                            break
                    elif jj in used:
                        good = False
                        break
                    else:
                        images[j], perms[j], used[jj] = jj, pj, j
                        todo.append(j)
            if good:
                require(len(images) == n, "Disconnected source")
                image = [images[t] for t in range(n)]
                perm = [perms[t] for t in range(n)]
                sign = verify_iso(source, target, image, perm)
                answers.append({"tetrahedron_images": image, "vertex_permutations": perm,
                                "orientation_sign": sign})
    return answers


def map_flows(source_flows, image, perms):
    return {(image[t], perms[t][v], perms[t][f]): z for (t, v, f), z in source_flows.items()}


def source_periods(raw, chart, deck, rows):
    """Integrate the given normal dual chains through exact source charts.

    In each lifted triangle a signed influx F_e follows its side midpoint
    towards the barycentre.  Sum F_e*(barycentre-midpoint_e); the barycentre
    term cancels because total influx is zero.  No orientation assumption
    or numerical approximation is needed for this calculation.
    """
    periods = chart["downstairs_periods"]
    edge_index = {tuple(e): i for i, e in enumerate(chart["torus_downstairs"]["edge_representatives"])}
    xy = {v: (Q(u), Q(w)) for v, u, w in chart["downstairs_vertex_coordinates"]}

    def period(a, b):
        e = tuple(sorted((a, b)))
        rep = min(e, tuple(sorted(deck[v] for v in e)))
        sign = 1 if (a, b) == rep or (deck[a], deck[b]) == rep else -1
        return tuple(sign * z for z in periods[edge_index[rep]])

    answer = [[Q(0), Q(0)] for _ in range(2)]
    for t in raw["cone_tetrahedra"]:
        labels = raw["tetrahedron_source_vertices"][t]
        v = labels.index(-1)
        corners = [w for w in range(4) if w != v]
        a = labels[corners[0]]
        lifted = {}
        for w in corners:
            b = labels[w]
            disp = (0, 0) if b == a else period(a, b)
            lifted[w] = tuple(z + d for z, d in zip(xy[min(b, deck[b])], disp))
        # Check that all three side displacements agree with these lifts.
        for w in corners:
            for z in corners:
                if w != z:
                    expected = tuple(b - a + d for a, b, d in
                                     zip(xy[min(labels[w], deck[labels[w]])],
                                         xy[min(labels[z], deck[labels[z]])],
                                         period(labels[w], labels[z])))
                    require(tuple(b - a for a, b in zip(lifted[w], lifted[z])) == expected,
                            "Source chart does not lift consistently over triangle")
        for curve in range(2):
            require(sum(rows[4 * t + 2 * curve][4 * v + f] for f in corners) == 0,
                    "Source curve not closed")
            for f in corners:
                flow = rows[4 * t + 2 * curve][4 * v + f]
                a, b = [w for w in corners if w != f]
                for k in range(2):
                    answer[curve][k] -= Q(flow, 2) * (lifted[a][k] + lifted[b][k])
    require(all(q.denominator == 1 for row in answer for q in row), "Nonintegral source period")
    return [[int(q) for q in row] for row in answer]


def audit_native_import(peripheral):
    raw_path = peripheral / "marked_coned_exterior.json.gz"
    native_path = peripheral / "native_import_peripheral_arrays.json.gz"
    chart_path = peripheral / "peripheral_chart.json.gz"
    cut_path = peripheral / "cut_source.json.gz"
    raw, native, chart, cut = [read_json(p) for p in (raw_path, native_path, chart_path, cut_path)]
    require(native["raw_source_sha256"] == canonical_hash(raw), "Native arrays refer to another raw source")
    require(raw["chart_sha256"] == canonical_hash(chart), "Raw chart hash mismatch")
    require(raw["cut_sha256"] == canonical_hash(cut), "Raw cut hash mismatch")
    check_gluings(raw)
    require(native["cusp_indices"] == raw["cusp_indices"], "Importer changed cusp labels")
    require(len(native["curves"]) == len(raw["curves"]), "Importer changed curve dimensions")
    changes = [0, 0, 0, 0]
    for i, (r, n) in enumerate(zip(raw["curves"], native["curves"])):
        expected = [-z for z in r] if i % 4 == 0 else r
        require(n == expected, "Importer performed more than the audited right-sheet M negation")
        changes[i % 4] += sum(a != b for a, b in zip(r, n))
    require(native["normalization"] == -1, "Unexpected normalization metadata")
    # The saved .tri must literally have these same arrays and gluings.
    imported = read_tri(peripheral / "marked_coned_exterior.tri")
    require(imported["data"] == raw["data"], "Saved import changed gluing arrays")
    require(imported["cusp_indices"] == native["cusp_indices"] and imported["curves"] == native["curves"],
            "Saved import does not contain the audited native arrays")
    complex_ = cusp_complex(raw)
    raw_curves = [cochain(complex_, get_flows(raw, k)) for k in range(2)]
    native_tri = {**raw, "curves": native["curves"]}
    native_curves = [cochain(complex_, get_flows(native_tri, k)) for k in range(2)]
    raw_pairing = wedge_integral(complex_, raw_curves[1], raw_curves[0])
    native_pairing = wedge_integral(complex_, native_curves[1], native_curves[0])
    require((raw_pairing, native_pairing) == (-1, 1), "Unexpected peripheral handedness")
    raw_periods = source_periods(raw, chart, cut["deck"], raw["curves"])
    native_periods = source_periods(raw, chart, cut["deck"], native["curves"])
    require(raw_periods == [[0, 1], [1, 0]], "Unexpected raw M,L source periods")
    require(native_periods == [[0, -1], [1, 0]], "Unexpected native M,L source periods")
    return {
        "raw_source_canonical_sha256": canonical_hash(raw),
        "native_arrays_file_sha256": file_hash(native_path),
        "raw_to_native_M_L_matrix": [[-1, 0], [0, 1]],
        "changed_entries_by_row_M_right_M_left_L_right_L_left": changes,
        "raw_M_L_periods_in_source_u_v": raw_periods,
        "native_M_L_periods_in_source_u_v": native_periods,
        "raw_pairing_L_M": raw_pairing,
        "native_pairing_L_M": native_pairing,
        "source_cusp_counts_vertices_edges_triangles": [len(complex_["vertices"]), len(complex_["edges"]), len(complex_["triangles"])],
        "diagram_coordinates_from_source_u_v": [[1, 0], [0, -1]],
        "coordinate_statement": "u_diagram=u_source; v_diagram=-v_source. Native M has diagram periods (0,1), native L has (1,0).",
        "trust_boundary": "Exact source/import audit and endpoint isomorphism; finite-vertex removal and intermediate simplification remain trusted native PL move transport.",
    }


def build_certificate(peripheral):
    source_path = peripheral / "matched_marked_exterior.tri"
    target_path = peripheral / "diagram_exterior.tri"
    source, target = read_tri(source_path), read_tri(target_path)
    require(all(c == 0 for row in source["cusp_indices"] for c in row), "Endpoint source not all ideal")
    require(all(c == 0 for row in target["cusp_indices"] for c in row), "Endpoint target not all ideal")
    source_complex, target_complex = cusp_complex(source), cusp_complex(target)
    source_flows = [get_flows(source, k) for k in range(2)]
    for f in source_flows:
        cochain(source_complex, f)
    target_curves = [cochain(target_complex, get_flows(target, k)) for k in range(2)]
    require(wedge_integral(target_complex, target_curves[1], target_curves[0]) == 1,
            "Target curves are not a positively normalized primitive peripheral basis")
    isos = enumerate_isos(source, target)
    require(isos, "No exact tetrahedron isomorphism")
    selected = None
    for index, iso in enumerate(isos):
        values = [cochain(target_complex, map_flows(f, iso["tetrahedron_images"], iso["vertex_permutations"]))
                  for f in source_flows]
        proofs = [exact_coordinates(target_complex, value, *target_curves) for value in values]
        cols = [p["target_M_L_coefficients"] for p in proofs]
        matrix = [[cols[j][i] for j in range(2)] for i in range(2)]
        determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        require(determinant == iso["orientation_sign"], "Cusp map determinant disagrees with tetrahedron orientation")
        iso.update({"cusp_matrix_columns_images_of_source_M_L": matrix,
                    "meridian_homology_certificate": proofs[0],
                    "longitude_homology_certificate": proofs[1]})
        if matrix == [[1, 0], [0, 1]] and iso["orientation_sign"] == 1:
            selected = index
    require(selected is not None, "No orientation-preserving isomorphism carrying M,L with their signs")
    audit = audit_native_import(peripheral)
    return {
        "schema": "exact-source-native-marked-cusp-isomorphism-v1",
        "source_file": str(source_path.relative_to(ROOT)),
        "source_file_sha256": file_hash(source_path),
        "target_file": str(target_path.relative_to(ROOT)),
        "target_file_sha256": file_hash(target_path),
        "tetrahedra": len(source["data"]),
        "conventions": {
            "normal_entries": "Positive enters the cusp triangle across the specified tetrahedron face.",
            "row_order": ["M/right", "M/left", "L/right", "L/left"],
            "right_sheet_oriented_corner_orders": CORNERS,
            "edge_cochain": "For each oriented right-sheet side a->b, c(a,b)=-normal_entry.",
            "potential_equation": "For source image minus a*M_target minus b*L_target, h(b)-h(a) equals its edge cochain. Equivalently the normal-flow difference is h(a)-h(b).",
            "matrix": "Column 0 is image(source M); column 1 is image(source L), in target (M,L).",
            "kernel_normalization": "Intersection(L,M)=+1; orient.c negates M alone when the initial value is -1.",
        },
        "target_cusp": {
            "vertices_as_germ_classes_tet_cuspvertex_othervertex": target_complex["vertices"],
            "oriented_edges": target_complex["edges"],
            "triangles": target_complex["triangles"],
            "pairing_L_M": 1,
        },
        "enumeration": "Exhaust every image of tetrahedron 0 and all 24 vertex permutations; propagate across every face and verify all gluing squares.",
        "all_isomorphisms": isos,
        "selected_identity_peripheral_isomorphism": selected,
        "native_import_audit": audit,
        "primary_convention_sources": {
            "snappy_release": "3.3.2_as_released", "commit": PINNED_COMMIT,
            "files": [{"url": f"https://raw.githubusercontent.com/3-manifolds/SnapPy/{PINNED_COMMIT}/{p}", "sha256": h}
                      for p, h in PRIMARY_SOURCES.items()],
        },
        "verified": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--peripheral-dir", type=Path, default=ROOT / "peripheral")
    parser.add_argument("--certificate", type=Path, default=HERE / "marked_isomorphism_certificate.json")
    args = parser.parse_args()
    result = build_certificate(args.peripheral_dir)
    # JSON conversion normalizes the tuples used internally.
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write_certificate:
        args.certificate.write_text(encoded)
    else:
        require(json.loads(args.certificate.read_text()) == json.loads(encoded), "Saved certificate differs from recomputation")
    selected = result["all_isomorphisms"][result["selected_identity_peripheral_isomorphism"]]
    print(json.dumps({"verified": True, "isomorphisms": len(result["all_isomorphisms"]),
                      "selected": selected,
                      "normalization": result["native_import_audit"]}, indent=2))


if __name__ == "__main__":
    main()

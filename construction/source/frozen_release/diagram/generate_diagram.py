#!/usr/bin/env python3
"""Build the explicitly marked signed pretzel diagram, without a knot-name lookup.

The source of the diagram is the indicated sum of rational tangles. Spherogram
constructs it independently of the direct twist-box wiring below. All outputs
are confined to this directory. No 3D cover or 4D object is claimed here.
"""
from __future__ import annotations
import hashlib
import inspect
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PARAMETERS = [-2, 3, 7]
SLOPES = [[-1, 2], [1, 3], [1, 7]]


def free(word):
    ans = []
    for a in word:
        if ans and ans[-1] == -a:
            ans.pop()
        else:
            ans.append(a)
    return ans


def inv(word):
    return [-a for a in reversed(word)]


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def direct_wiring():
    crossing_ids = [f"R{r}C{i}" for r, q in enumerate(PARAMETERS) for i in range(abs(q))]
    adjacency = {c: [None] * 4 for c in crossing_ids}
    ports, centres, boundaries, connections = {}, {}, [], []

    def connect(a, b, points, kind):
        a, b = tuple(a), tuple(b)
        if adjacency[a[0]][a[1]] is not None or adjacency[b[0]][b[1]] is not None:
            raise ValueError("port used twice")
        adjacency[a[0]][a[1]] = b
        adjacency[b[0]][b[1]] = a
        connections.append({"ends": [a, b], "points": points, "kind": kind})

    for r, q in enumerate(PARAMETERS):
        ids = [f"R{r}C{i}" for i in range(abs(q))]
        tl, tr, bl, br = (3, 2, 0, 1) if q > 0 else (2, 1, 3, 0)
        x = 6 * r
        for i, c in enumerate(ids):
            y = -2 - 3 * i
            centres[c] = [x, y]
            ports[c] = {tl: [x-1, y+1, 0], tr: [x+1, y+1, 0],
                        bl: [x-1, y-1, 0], br: [x+1, y-1, 0]}
        for a, b in zip(ids, ids[1:]):
            connect((a, bl), (b, tl), [ports[a][bl], ports[b][tl]], "within_twist_box")
            connect((a, br), (b, tr), [ports[a][br], ports[b][tr]], "within_twist_box")
        bd = {"TL": (ids[0], tl), "TR": (ids[0], tr),
              "BL": (ids[-1], bl), "BR": (ids[-1], br)}
        boundaries.append(bd)

    def terminal(r, name):
        return [6*r + (-1 if name.endswith("L") else 1), 0 if name.startswith("T") else -22, 0]

    for r in range(2):
        for side in ("T", "B"):
            a, b = boundaries[r][side+"R"], boundaries[r+1][side+"L"]
            connect(a, b, [ports[a[0]][a[1]], terminal(r, side+"R"),
                           terminal(r+1, side+"L"), ports[b[0]][b[1]]], "horizontal_tangle_sum")
    for side, height in (("T", 2), ("B", -24)):
        a, b = boundaries[0][side+"L"], boundaries[2][side+"R"]
        connect(a, b, [ports[a[0]][a[1]], terminal(0, side+"L"), [-1, height, 0],
                       [13, height, 0], terminal(2, side+"R"), ports[b[0]][b[1]]],
                "numerator_" + side.lower() + "_closure")
    return crossing_ids, adjacency, ports, centres, boundaries, connections


def make_record():
    import spherogram
    from spherogram import RationalTangle, Crossing
    from spherogram.links import tangles, links_base
    from sympy import Matrix, symbols, Poly

    ts = []
    for r, (a, b) in enumerate(SLOPES):
        t = RationalTangle(a, b)
        for i, c in enumerate(x for x in t.crossings if isinstance(x, Crossing)):
            c.label = f"R{r}C{i}"
        ts.append(t)
    link = (ts[0] + ts[1] + ts[2]).numerator_closure()
    ids, adj, ports, centres, boundaries, connections = direct_wiring()
    norm = {c.label: [(d.label, j) for d, j in c.adjacent] for c in link.crossings}
    rotations = []
    for ds in itertools.product((0, 2), repeat=len(ids)):
        d = dict(zip(ids, ds))
        if all(norm[c][(i+d[c]) % 4] == (v, (j+d[v]) % 4)
               for c in ids for i, (v, j) in enumerate(adj[c])):
            rotations.append(d)
    if len(rotations) != 1:
        raise ValueError("direct twist boxes did not identify uniquely with the rational-tangle source")
    rotation = rotations[0]

    # The chosen direction agrees with the source library, but is specified here
    # by an actual port: from right to left along the outside top closure.
    start = boundaries[0]["TL"]
    entries, seen, cur = [], set(), start
    while cur not in seen:
        seen.add(cur)
        entries.append(cur)
        c, i = cur
        cur = tuple(adj[c][(i+2) % 4])
    if cur != start or len(entries) != 2 * len(ids):
        raise ValueError("not one knot component")
    library_entries = {(e.crossing.label, (e.strand_index-rotation[e.crossing.label]) % 4)
                       for e in link.link_components[0]}
    if len(link.link_components) != 1 or set(entries) != library_entries:
        raise ValueError("source orientation mismatch")
    local = {}
    for c in ids:
        under = next(i for v, i in entries if v == c and i % 2 == 0)
        over = next(i for v, i in entries if v == c and i % 2 == 1)
        sign = 1 if (under-over) % 4 == 1 else -1
        local[c] = {"under_in": under, "over_in": over, "sign": sign}
    by_id = {c.label: c for c in link.crossings}
    if any(local[c]["sign"] != by_id[c].sign for c in ids):
        raise ValueError("crossing-sign mismatch")

    arc = 1
    traversal = []
    under_order = []
    for c, i in entries:
        old = arc
        if i % 2 == 0:
            arc = arc % len(ids) + 1
            local[c].update({"incoming_arc": old, "outgoing_arc": arc})
            under_order.append(c)
        else:
            local[c]["over_arc"] = arc
        traversal.append({"crossing": c, "entry_port": i, "exit_port": (i+2) % 4,
                          "kind": "under" if i % 2 == 0 else "over",
                          "incoming_arc": old, "outgoing_arc": arc})
    if arc != 1:
        raise ValueError("arc labelling did not close")

    relators, blackboard, proof = [], [], []
    for c in under_order:
        z = local[c]
        letter = z["sign"] * z["over_arc"]
        rel = [-letter, z["incoming_arc"], letter, -z["outgoing_arc"]]
        relators.append({"crossing": c, "word": rel})
        blackboard.append(letter)
        proof.append({"relator_index": len(relators)-1, "conjugator": list(blackboard)})
    writhe = sum(local[c]["sign"] for c in ids)
    correction = [1 if -writhe >= 0 else -1] * abs(writhe)
    longitude = free(blackboard + correction)
    certificate = []
    for p in proof:
        w = p["conjugator"]
        certificate += w + relators[p["relator_index"]]["word"] + inv(w)
    commutator = [1] + longitude + [-1] + inv(longitude)
    if free(certificate) != free(commutator):
        raise ValueError("peripheral commutator certificate failure")

    lookup = {}
    for e in connections:
        a, b = map(tuple, e["ends"])
        lookup[a] = e["points"]
        lookup[b] = list(reversed(e["points"]))
    polygon = [ports[start[0]][start[1]]]
    for c, i in entries:
        centre = centres[c] + [F(-1, 4) if i % 2 == 0 else F(1, 4)]
        out = (c, (i+2) % 4)
        polygon.extend([centre, ports[c][out[1]]])
        polygon.extend(lookup[out][1:])
    if polygon[0] != polygon[-1]:
        raise ValueError("geometric polygon did not close")
    polygon.pop()
    marked = [6, 2, 0]
    hit = next(i for i in range(len(polygon))
               if polygon[i] == [13, 2, 0] and polygon[(i+1) % len(polygon)] == [-1, 2, 0])
    polygon.insert(hit+1, marked)
    hit += 1
    polygon = polygon[hit:] + polygon[:hit]
    polygon = [[str(F(q)) for q in p] for p in polygon]

    colour = [[0]*len(ids) for _ in ids]
    for row, c in enumerate(under_order):
        z = local[c]
        colour[row][z["over_arc"]-1] += 2
        colour[row][z["incoming_arc"]-1] -= 1
        colour[row][z["outgoing_arc"]-1] -= 1
    colour_minor_det = int(Matrix([r[1:] for r in colour[1:]]).det())
    V = link.seifert_matrix()
    t = symbols("t")
    alex = Poly((t*Matrix(V)-Matrix(V).transpose()).det(), t)
    seifert_fill_matrix = [[2, 0, 0, -1], [0, 3, 0, 1], [0, 0, 7, 1], [1, 1, 1, 0]]

    record = {
        "schema": "marked-pretzel-diagram-v1",
        "classification": "EXPLICIT_ORIENTED_DIAGRAM_AND_PERIPHERAL_WORDS",
        "knot_label": "J=P(-2,3,7)",
        "source": {
            "construction": "numerator closure of T(-1/2)+T(1/3)+T(1/7), in that left-to-right order",
            "signed_pretzel_parameters": PARAMETERS,
            "rational_tangle_slopes": SLOPES,
            "library_expression": "(RationalTangle(-1,2)+RationalTangle(1,3)+RationalTangle(1,7)).numerator_closure()",
            "spherogram_version": spherogram.__version__,
            "library_source_hashes": {Path(inspect.getfile(m)).name: hashlib.sha256(Path(inspect.getfile(m)).read_bytes()).hexdigest()
                                      for m in (tangles, links_base)},
            "uses_knot_name_lookup": False,
            "conventions_url": "https://snappy.computop.org/spherogram.html",
            "convention_note": "The official tutorial explicitly names this rational-tangle sum the (-2,3,7) pretzel. Its oriented crossings here all have Spherogram sign -1. No silent mirroring is performed."
        },
        "raw_crossing_convention": {
            "slots": "0,1,2,3 in counterclockwise order; 0--2 is the underpass and 1--3 the overpass",
            "positive_box_physical_ports": {"SW": 0, "SE": 1, "NE": 2, "NW": 3},
            "negative_box_physical_ports": {"SW": 3, "SE": 0, "NE": 1, "NW": 2},
            "coordinate_orientation": "standard (x,y,z); diagram viewed from positive z; overpass has z=+1/4",
            "oriented_sign": "+1 when (under_entry-over_entry) mod 4 is 1, otherwise -1",
            "source_to_normalized_slot": "normalized_slot = (raw_slot + recorded_rotation[crossing]) mod 4",
        },
        "crossing_ids": ids,
        "raw_adjacency": adj,
        "tangle_boundary_ports": boundaries,
        "library_normalization_rotation": rotation,
        "pd_code": [list(row) for row in link.PD_code()],
        "normalized_adjacency": norm,
        "oriented_crossings": local,
        "orientation": {"start_entry": start, "description": "right to left along the outside top closure, then into the top-left port of the -2 box",
                        "traversal": traversal, "component_count": 1, "writhe": writhe},
        "spatial_diagram": {
            "coordinate_field": "Q", "units": "arbitrary consistent diagram units",
            "crossing_centres_xy": centres,
            "crossing_ports_xyz": {c: [ports[c][i] for i in range(4)] for c in ids},
            "connection_polylines": connections,
            "oriented_polygon_vertices": polygon,
            "closed_by": "last vertex to first vertex",
            "geometry_purpose": "an explicit PL embedding realizing the marked diagram; not a triangulation of its complement"
        },
        "marking": {
            "knot_basepoint_xyz": ["6", "2", "0"],
            "oriented_tangent": [-1, 0, 0],
            "diagram_arc_generator": 1,
            "pd_strand_label": int(by_id[start[0]].strand_labels[(start[1]+rotation[start[0]]) % 4]),
            "raw_diagram_edge_endpoints": [boundaries[2]["TR"], boundaries[0]["TL"]],
            "puncture_ball": {"shape": "closed axis-aligned PL cube", "centre": ["6", "2", "0"],
                              "half_width": "1/4", "branch_arc_endpoints_in_oriented_order": [["25/4", "2", "0"], ["23/4", "2", "0"]],
                              "statement": "The cube meets the knot in its straight central x-directed arc only."},
            "positive_normal_basis": [[0, 1, 0], [0, 0, -1]],
            "positive_meridian_polygon": [["6", "17/8", "0"], ["6", "2", "-1/8"], ["6", "15/8", "0"], ["6", "2", "1/8"]],
            "meridian_orientation": "normal basis follows tangent to give the ambient orientation; this meridian has linking number +1",
            "peripheral_basepoint_xyz": ["6", "17/8", "0"],
            "whisker_convention": "Wirtinger stems to the upper half-space, starting at the marked top-closure arc; meridians are oriented with positive linking number",
        },
        "wirtinger": {
            "generator_count": len(ids), "generator_meaning": "positive meridian of each successive under-to-under arc; arc 1 contains the marked point",
            "under_crossing_order": under_order,
            "relation_convention": "o^(-epsilon) * incoming * o^(epsilon) * outgoing^(-1) = 1",
            "relators": relators,
            "meridian_word": [1], "blackboard_longitude_word": blackboard,
            "blackboard_linking_number": writhe,
            "zero_longitude_raw_word": blackboard + correction,
            "zero_longitude_word": longitude,
            "zero_longitude_geometric_prescription": "Follow the oriented right-side blackboard push-off once from the peripheral basepoint, then append -writhe positive meridian turns; this corrects the push-off linking number to zero.",
            "commutator_convention": "[mu,lambda]=mu*lambda*mu^-1*lambda^-1",
            "commutator_relator_product": proof,
            "certificate_meaning": "product of conjugator*relator*conjugator^-1 freely equals [mu,lambda]",
        },
        "invariant_checks": {
            "fox_colouring_matrix": colour,
            "deleted_row_column": [0, 0], "fox_minor_determinant": colour_minor_det,
            "knot_determinant": abs(colour_minor_det),
            "spherogram_seifert_matrix": V,
            "alexander_polynomial_coefficients_ascending": [int(alex.nth(i)) for i in range(alex.degree()+1)],
            "alexander_at_minus_one": int(alex.eval(-1)),
            "role": "cross-checks of the explicit construction, not knot-identification substitutes"
        },
        "seifert_cover_interface": {
            "status": "SOURCE_THEOREM_AND_EXACT_FILLING_PRESCRIPTION; NO_CELLWISE_PAIR_MAP_IN_THIS_FILE",
            "base": "oriented pair of pants", "fibre": "oriented S1 with generator h",
            "boundary_section_convention": "c1+c2+c3=0 in H1; group relation c1*c2*c3=1 before filling",
            "filling_pairs_alpha_beta": [[2, -1], [3, 1], [7, 1]],
            "filling_meridians": ["2*c1-h", "3*c2+h", "7*c3+h"],
            "filling_relation_convention": "ci^alpha_i * h^beta_i = 1, with h central",
            "abelianization_matrix_columns": ["c1", "c2", "c3", "h"],
            "abelianization_matrix": seifert_fill_matrix,
            "abelianization_matrix_determinant": int(Matrix(seifert_fill_matrix).det()),
            "sum_beta_over_alpha": str(sum(F(a, b) for a, b in SLOPES)),
            "normalized_description": {"b": -1, "pairs": [[2, 1], [3, 1], [7, 1]], "convention": "q1*q2*q3=h^b"},
            "basis_transport_requirement": "A source-bound branched-cover implementation must transport these pillowcase/fibre bases and the diagram marks; homology or a named Seifert manifold is insufficient.",
            "reference": "https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf"
        },
        "scope": {
            "completed_here": ["explicit signed diagram", "source-to-direct-wiring equivalence", "rational PL knot embedding", "orientation and trivial-arc puncture", "diagrammatic meridian and zero longitude with group certificate", "determinant cross-check"],
            "not_completed_here": ["triangulated branched double cover", "cellwise source-pair and peripheral transport", "relative roll map", "even-twist bulk", "beta0/beta2 caps", "X2 or X0", "standardization"],
            "handedness": "This is the exact rational-tangle convention documented above. Matching to any paper using another pretzel convention requires an explicit convention comparison."
        }
    }
    return record


def write_svg(record, path):
    def xy(p):
        return (70 + 32*float(F(p[0])), 126 - 24*float(F(p[1])))
    def pts(ps):
        return " ".join(f"{x:.3f},{y:.3f}" for x, y in map(xy, ps))
    body = ['<svg xmlns="http://www.w3.org/2000/svg" width="570" height="760" viewBox="0 0 570 760">',
            '<rect width="570" height="760" fill="white"/>',
            '<g font-family="sans-serif" fill="#172033"><text x="24" y="24" font-size="18">Marked J = P(−2,3,7)</text>',
            '<text x="24" y="46" font-size="11">T(−1/2) + T(1/3) + T(1/7), numerator closure</text></g>',
            '<g fill="none" stroke="#172033" stroke-width="2.5" stroke-linejoin="round">']
    for e in record["spatial_diagram"]["connection_polylines"]:
        body.append(f'<polyline points="{pts(e["points"])}"/>')
    for c in record["crossing_ids"]:
        p = record["spatial_diagram"]["crossing_ports_xyz"][c]
        body.append(f'<polyline points="{pts([p[0], p[2]])}"/>')
        body.append(f'<polyline points="{pts([p[1], p[3]])}" stroke="white" stroke-width="8"/>')
        body.append(f'<polyline points="{pts([p[1], p[3]])}"/>')
    body.append('</g>')
    for r, slope in enumerate(("−1/2", "1/3", "1/7")):
        x, _ = xy([6*r, 0])
        body.append(f'<text x="{x}" y="729" text-anchor="middle" font-family="sans-serif" font-size="14">T({slope})</text>')
    x, y = xy([6, 2])
    body.extend([f'<rect x="{x-8}" y="{y-6}" width="16" height="12" fill="none" stroke="#14866d" stroke-dasharray="3 2"/>',
                 f'<circle cx="{x}" cy="{y}" r="3.5" fill="#14866d"/>',
                 f'<path d="M {x-12} {y} l 8 -4 l 0 8 Z" fill="#14866d"/>',
                 f'<text x="{x+16}" y="{y-8}" font-family="sans-serif" font-size="12" fill="#14866d">marked trivial arc</text>',
                 '<text x="24" y="752" font-family="sans-serif" font-size="10" fill="#475569">Exact source wiring; all 12 oriented crossing signs are −1 in Spherogram conventions.</text>', '</svg>'])
    path.write_text("\n".join(body)+"\n")


def main():
    record = make_record()
    target = ROOT / "pretzel_237_diagram.json"
    target.write_text(json.dumps(record, indent=2)+"\n")
    write_svg(record, ROOT / "pretzel_237_diagram.svg")
    print(json.dumps({"output": str(target), "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                      "crossings": len(record["crossing_ids"]), "components": 1,
                      "writhe": record["orientation"]["writhe"],
                      "determinant": record["invariant_checks"]["knot_determinant"],
                      "classification": record["classification"]}, indent=2))


if __name__ == "__main__":
    main()

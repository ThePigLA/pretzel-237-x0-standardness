#!/usr/bin/env python3
"""Independent stdlib checks of the signed diagram, PL embedding and markings.

This does not import the generator. Optional --with-library regenerates the
documented rational-tangle expression and compares the oriented PD record.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def reduce_word(word):
    result = []
    for v in word:
        require(isinstance(v, int) and v != 0, "invalid free-group letter")
        if result and result[-1] == -v:
            result.pop()
        else:
            result.append(v)
    return result


def reverse_word(word):
    return [-v for v in word[::-1]]


def determinant(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    require(all(len(row) == n for row in a), "not square")
    sign, previous = 1, 1
    if not n:
        return 1
    for k in range(n-1):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[k], a[p] = a[p], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                z, rem = divmod(a[i][j]*pivot-a[i][k]*a[k][j], previous)
                require(rem == 0, "nonintegral Bareiss division")
                a[i][j] = z
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def expected_wiring():
    # Explicit local rational-tangle convention, independent of Spherogram and
    # of the producer: each 1/q tangle is a stack of q crossings in a rectangle.
    labels = [f"R{r}C{i}" for r, n in enumerate((2, 3, 7)) for i in range(n)]
    pair = {}
    boundary = []
    def glue(a, b):
        require(a not in pair and b not in pair, "duplicate expected endpoint")
        pair[a], pair[b] = b, a
    for r, n in enumerate((2, 3, 7)):
        tl, tr, bl, br = (2, 1, 3, 0) if r == 0 else (3, 2, 0, 1)
        names = [f"R{r}C{i}" for i in range(n)]
        for i in range(n-1):
            glue((names[i], bl), (names[i+1], tl))
            glue((names[i], br), (names[i+1], tr))
        boundary.append({"TL": [names[0], tl], "TR": [names[0], tr],
                         "BL": [names[-1], bl], "BR": [names[-1], br]})
    for a, b in ((0, 1), (1, 2)):
        for k in ("T", "B"):
            glue(tuple(boundary[a][k+"R"]), tuple(boundary[b][k+"L"]))
    glue(tuple(boundary[0]["TL"]), tuple(boundary[2]["TR"]))
    glue(tuple(boundary[0]["BL"]), tuple(boundary[2]["BR"]))
    return labels, pair, boundary


def cross2(u, v):
    return u[0]*v[1]-u[1]*v[0]


def on_segment(p, a, b):
    return cross2([p[0]-a[0], p[1]-a[1]], [b[0]-a[0], b[1]-a[1]]) == 0 and all(
        min(a[k], b[k]) <= p[k] <= max(a[k], b[k]) for k in range(2))


def projected_intersection(a, b, c, d):
    u = [b[k]-a[k] for k in range(2)]
    v = [d[k]-c[k] for k in range(2)]
    w = [c[k]-a[k] for k in range(2)]
    det = cross2(u, v)
    if det:
        t, s = cross2(w, v)/det, cross2(w, u)/det
        if 0 <= t <= 1 and 0 <= s <= 1:
            return [(a[0]+t*u[0], a[1]+t*u[1])]
        return []
    if cross2(w, u):
        return []
    return sorted(set(tuple(p[:2]) for p in (a, b, c, d)
                      if on_segment(p, a, b) and on_segment(p, c, d)))


def height_at(p, a, b):
    k = next(k for k in range(2) if a[k] != b[k])
    t = (p[k]-a[k])/(b[k]-a[k])
    return a[2]+t*(b[2]-a[2])


def cube_intersection(a, b, centre, radius):
    lo, hi = Q(0), Q(1)
    for k in range(3):
        v = b[k]-a[k]
        if v == 0:
            if not centre[k]-radius <= a[k] <= centre[k]+radius:
                return None
        else:
            p, q = sorted(((centre[k]-radius-a[k])/v, (centre[k]+radius-a[k])/v))
            lo, hi = max(lo, p), min(hi, q)
            if lo > hi:
                return None
    return [[a[k]+t*(b[k]-a[k]) for k in range(3)] for t in (lo, hi)]


def verify(data, with_library=False):
    require(data["schema"] == "marked-pretzel-diagram-v1", "wrong schema")
    require(data["source"]["signed_pretzel_parameters"] == [-2, 3, 7], "wrong signed pretzel")
    require(data["source"]["rational_tangle_slopes"] == [[-1, 2], [1, 3], [1, 7]], "wrong rational slopes")
    labels, expected, bd = expected_wiring()
    require(data["crossing_ids"] == labels, "crossing labels changed")
    require(data["tangle_boundary_ports"] == bd, "boundary positions changed")
    pair = {(c, i): tuple(v) for c, row in data["raw_adjacency"].items() for i, v in enumerate(row)}
    require(pair == expected, "wiring does not realize the prescribed rational-tangle sum")
    require(all(pair[v] == k and k != v for k, v in pair.items()), "not a fixed-point-free dart pairing")
    require(len(pair) == 48, "wrong half-edge count")

    # Rotation-system check proves that this connected diagram graph is on S2.
    unseen, faces = set(pair), []
    while unseen:
        first = cur = min(unseen)
        f = []
        while cur in unseen:
            unseen.remove(cur)
            f.append(cur)
            c, i = pair[cur]
            cur = (c, (i+1) % 4)
        require(cur == first, "face walk failed")
        faces.append(f)
    require(len(labels)-len(pair)//2+len(faces) == 2, "rotation system is not spherical")

    rot = data["library_normalization_rotation"]
    require(set(rot) == set(labels) and set(rot.values()) <= {0, 2}, "normalization must preserve underpass")
    normalized = {(c, (i+rot[c]) % 4): (d, (j+rot[d]) % 4) for (c, i), (d, j) in pair.items()}
    supplied_normal = {(c, i): tuple(v) for c, row in data["normalized_adjacency"].items() for i, v in enumerate(row)}
    require(normalized == supplied_normal, "normalization is not a graph isomorphism")
    occurrences = {}
    for c, row in zip(labels, data["pd_code"]):
        require(len(row) == 4, "bad PD crossing")
        for i, e in enumerate(row):
            occurrences.setdefault(e, []).append((c, i))
    require(sorted(occurrences) == list(range(24)), "bad PD strand labels")
    require(all(len(z) == 2 and normalized[z[0]] == z[1] for z in occurrences.values()), "PD disagrees with wiring")

    start = tuple(bd[0]["TL"])
    require(tuple(data["orientation"]["start_entry"]) == start, "wrong oriented start")
    loop, cur = [], start
    while cur not in loop:
        loop.append(cur)
        cur = pair[cur[0], (cur[1]+2) % 4]
    require(cur == start and len(loop) == 24, "not a single knot traversing every strand")
    require(len({(c, i % 2) for c, i in loop}) == 24, "duplicate strand pass")
    require([(x["crossing"], x["entry_port"]) for x in data["orientation"]["traversal"]] == loop, "wrong recorded orientation")
    local, arc, unders = {}, 1, []
    for c in labels:
        u = next(i for x, i in loop if x == c and i % 2 == 0)
        o = next(i for x, i in loop if x == c and i % 2 == 1)
        local[c] = {"under_in": u, "over_in": o, "sign": 1 if (u-o) % 4 == 1 else -1}
    for (c, i), event in zip(loop, data["orientation"]["traversal"]):
        old = arc
        if i % 2 == 0:
            arc = arc % 12 + 1
            local[c].update(incoming_arc=old, outgoing_arc=arc)
            unders.append(c)
        else:
            local[c]["over_arc"] = arc
        require(event["incoming_arc"] == old and event["outgoing_arc"] == arc and event["exit_port"] == (i+2) % 4,
                "incorrect oriented arc event")
    require(arc == 1 and local == data["oriented_crossings"], "incorrect Wirtinger arcs/signs")
    writhe = sum(c["sign"] for c in local.values())
    require(writhe == -12 == data["orientation"]["writhe"], "handedness or writhe mismatch")
    require(data["marking"]["pd_strand_label"] == data["pd_code"][0][(start[1]+rot[start[0]]) % 4], "marked PD strand mismatch")
    require(data["marking"]["pd_strand_label"] == 22, "changed marked source strand")

    # Check the geometry against the signed port convention before checking its
    # embeddedness: a reflected picture would otherwise pass intersection tests.
    geo = data["spatial_diagram"]
    ports = {c: [list(map(Q, p)) for p in row] for c, row in geo["crossing_ports_xyz"].items()}
    centres = {}
    for r, n in enumerate((2, 3, 7)):
        loc = {"SW": 3, "SE": 0, "NE": 1, "NW": 2} if r == 0 else {"SW": 0, "SE": 1, "NE": 2, "NW": 3}
        for i in range(n):
            c, x, y = f"R{r}C{i}", 6*r, -2-3*i
            centres[c] = (Q(x), Q(y))
            require(tuple(map(Q, geo["crossing_centres_xy"][c])) == centres[c], "wrong crossing centre")
            for corner, slot in loc.items():
                expect = [Q(x + (-1 if corner[1] == "W" else 1)), Q(y + (1 if corner[0] == "N" else -1)), Q(0)]
                require(ports[c][slot] == expect, "reflected or altered crossing chart")
    connectors = {}
    for e in geo["connection_polylines"]:
        a, b = map(tuple, e["ends"])
        path = [list(map(Q, p)) for p in e["points"]]
        require(pair[a] == b and a not in connectors and b not in connectors, "connection graph mismatch")
        require(path[0] == ports[a[0]][a[1]] and path[-1] == ports[b[0]][b[1]], "connection endpoint mismatch")
        connectors[a], connectors[b] = path, path[::-1]
    require(set(connectors) == set(pair), "missing connector")
    polygon = [ports[start[0]][start[1]]]
    for c, i in loop:
        out = (c, (i+2) % 4)
        polygon += [list(centres[c])+[Q(-1, 4) if i % 2 == 0 else Q(1, 4)], ports[c][out[1]]]
        polygon += connectors[out][1:]
    require(polygon[-1] == polygon[0], "polygon is open")
    polygon.pop()
    mark = [Q(6), Q(2), Q(0)]
    i = next(i for i in range(len(polygon)) if polygon[i] == [13, 2, 0] and polygon[(i+1) % len(polygon)] == [-1, 2, 0])
    polygon.insert(i+1, mark)
    polygon = polygon[i+1:]+polygon[:i+1]
    stored = [list(map(Q, p)) for p in geo["oriented_polygon_vertices"]]
    require(polygon == stored and len(set(map(tuple, polygon))) == len(polygon), "wrong or repeated polygon vertices")
    segments = list(zip(polygon, polygon[1:]+polygon[:1]))
    expected_crossings = set(centres.values())
    found = set()
    intersection_pairs = 0
    for i, (a, b) in enumerate(segments):
        require(a[:2] != b[:2], "vertical or zero projected segment")
        for j in range(i+1, len(segments)):
            if j == i+1 or (i == 0 and j == len(segments)-1):
                continue
            c, d = segments[j]
            hits = projected_intersection(a, b, c, d)
            require(len(hits) <= 1, "overlapping projected nonadjacent segments")
            for p in hits:
                require(p in expected_crossings, "unintended projected crossing")
                h1, h2 = height_at(p, a, b), height_at(p, c, d)
                require({h1, h2} == {Q(-1, 4), Q(1, 4)}, "spatial crossing or wrong over-under height")
                found.add(p)
                intersection_pairs += 1
    require(found == expected_crossings, "missing projected crossing")
    require(intersection_pairs == 48, "unexpected PL crossing segment incidences")
    cube = data["marking"]["puncture_ball"]
    require(list(map(Q, cube["centre"])) == mark and Q(cube["half_width"]) == Q(1, 4), "changed puncture cube")
    clips = [hit for a, b in segments if (hit := cube_intersection(a, b, mark, Q(1, 4))) is not None]
    require(len(clips) == 2 and all(p[1:] == [Q(2), Q(0)] for hit in clips for p in hit), "puncture contains other knot pieces")
    require(sorted(set(p[0] for hit in clips for p in hit)) == [Q(23, 4), Q(6), Q(25, 4)], "wrong puncture arc endpoints")
    require(polygon[1][0] < polygon[0][0], "wrong marked orientation")
    require(data["marking"]["oriented_tangent"] == [-1, 0, 0] and data["marking"]["positive_normal_basis"] == [[0, 1, 0], [0, 0, -1]], "normal orientation mismatch")
    meridian = [list(map(Q, p)) for p in data["marking"]["positive_meridian_polygon"]]
    require(meridian == [[Q(6), Q(17, 8), Q(0)], [Q(6), Q(2), Q(-1, 8)], [Q(6), Q(15, 8), Q(0)], [Q(6), Q(2), Q(1, 8)]], "changed positive meridian")
    meridian_normal = [(p[1]-2, -p[2]) for p in meridian]
    require(sum(cross2(a, b) for a, b in zip(meridian_normal, meridian_normal[1:]+meridian_normal[:1])) > 0,
            "meridian does not wind positively in the declared normal disk")

    w = data["wirtinger"]
    require(w["generator_count"] == 12 and w["under_crossing_order"] == unders and w["meridian_word"] == [1], "wrong marked group presentation")
    rels, bb = [], []
    for c in unders:
        z = local[c]
        k = z["sign"]*z["over_arc"]
        rels.append({"crossing": c, "word": [-k, z["incoming_arc"], k, -z["outgoing_arc"]]})
        bb.append(k)
    require(w["relators"] == rels and w["blackboard_longitude_word"] == bb, "peripheral words not derived from this diagram")
    require(w["blackboard_linking_number"] == writhe, "wrong push-off framing")
    longitude = reduce_word(bb+[1]*12)
    require(w["zero_longitude_raw_word"] == bb+[1]*12 and w["zero_longitude_word"] == longitude, "missing zero-framing correction")
    require(sum(1 if a > 0 else -1 for a in longitude) == 0, "longitude nonzero in H1")
    product = []
    for item in w["commutator_relator_product"]:
        c = item["conjugator"]
        product += c + rels[item["relator_index"]]["word"] + reverse_word(c)
    comm = [1]+longitude+[-1]+reverse_word(longitude)
    require(reduce_word(product) == reduce_word(comm), "commutator proof does not freely telescope")

    colour = [[0]*12 for _ in range(12)]
    for i, c in enumerate(unders):
        z = local[c]
        colour[i][z["over_arc"]-1] += 2
        colour[i][z["incoming_arc"]-1] -= 1
        colour[i][z["outgoing_arc"]-1] -= 1
    facts = data["invariant_checks"]
    require(colour == facts["fox_colouring_matrix"], "wrong colouring matrix")
    det = determinant([row[1:] for row in colour[1:]])
    require(det == facts["fox_minor_determinant"] == -1 and facts["knot_determinant"] == 1, "determinant mismatch")
    V = facts["spherogram_seifert_matrix"]
    require(len(V) == 10 and all(len(row) == 10 for row in V), "unexpected Seifert-matrix dimensions")
    require(determinant([[V[i][j]+V[j][i] for j in range(10)] for i in range(10)]) == -1, "Seifert determinant cross-check failed")
    coefficients = facts["alexander_polynomial_coefficients_ascending"]
    require(len(coefficients) == 11, "wrong Alexander degree bound")
    # Both polynomials have degree at most ten, so eleven exact evaluations
    # verify their equality, without SymPy or the producer's determinant code.
    for value in range(-5, 6):
        evaluated = determinant([[value*V[i][j]-V[j][i] for j in range(10)] for i in range(10)])
        require(evaluated == sum(a*value**i for i, a in enumerate(coefficients)), "Alexander polynomial mismatch")
    require(sum(a*(-1)**i for i, a in enumerate(coefficients)) == facts["alexander_at_minus_one"] == -1,
            "Alexander value mismatch")
    s = data["seifert_cover_interface"]
    require(s["filling_pairs_alpha_beta"] == [[2, -1], [3, 1], [7, 1]], "wrong source filling slopes")
    mat = [[2, 0, 0, -1], [0, 3, 0, 1], [0, 0, 7, 1], [1, 1, 1, 0]]
    require(s["abelianization_matrix"] == mat and determinant(mat) == s["abelianization_matrix_determinant"] == 1,
            "wrong Seifert H1 check")
    require(sum(Q(b, a) for a, b in s["filling_pairs_alpha_beta"]) == Q(s["sum_beta_over_alpha"]) == Q(-1, 42), "wrong rational Euler sum")

    library = "NOT_REQUESTED"
    if with_library:
        from spherogram import RationalTangle
        knot = (RationalTangle(-1, 2)+RationalTangle(1, 3)+RationalTangle(1, 7)).numerator_closure()
        require([list(row) for row in knot.PD_code()] == data["pd_code"], "regenerated rational-tangle PD differs")
        require(len(knot.link_components) == 1 and [c.sign for c in knot.crossings] == [-1]*12, "library handedness mismatch")
        require(knot.seifert_matrix() == V, "library Seifert-matrix mismatch")
        library = "PASS"
    return {"classification": "DIAGRAM_AND_MARKINGS_VERIFIED", "crossings": 12, "components": 1,
            "diagram_faces": len(faces), "writhe": writhe, "knot_determinant": 1,
            "polygon_vertices": len(polygon), "unwanted_projected_crossings": 0,
            "spatial_self_intersections": 0, "puncture_intersected_segments": len(clips),
            "meridian_pd_strand": 22, "longitude_exponent_sum": 0,
            "peripheral_commutator_relator_product": "PASS", "seifert_filling_matrix_determinant": 1,
            "optional_rational_tangle_regeneration": library,
            "not_certified": ["triangulated branched double cover", "cellwise source-pair map", "peripheral transport to a 3D triangulation", "X2", "X0"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=Path(__file__).with_name("pretzel_237_diagram.json"))
    parser.add_argument("--with-library", action="store_true")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    result = verify(data, args.with_library)
    # These are rejection controls for this exact source, not old control reruns.
    mutations = []
    for name in ("mirror_slope", "wrong_zero_longitude", "wrong_puncture", "over_under_chart"):
        z = copy.deepcopy(data)
        if name == "mirror_slope":
            z["source"]["rational_tangle_slopes"][0][0] = 1
        elif name == "wrong_zero_longitude":
            z["wirtinger"]["zero_longitude_word"] = z["wirtinger"]["blackboard_longitude_word"]
        elif name == "wrong_puncture":
            z["marking"]["puncture_ball"]["half_width"] = "10"
        else:
            z["spatial_diagram"]["crossing_ports_xyz"]["R0C0"][0], z["spatial_diagram"]["crossing_ports_xyz"]["R0C0"][1] = z["spatial_diagram"]["crossing_ports_xyz"]["R0C0"][1], z["spatial_diagram"]["crossing_ports_xyz"]["R0C0"][0]
        try:
            verify(z)
        except (ValueError, KeyError, StopIteration) as exc:
            mutations.append({"name": name, "result": "REJECTED", "reason": str(exc)})
        else:
            raise ValueError("mutation was admitted: "+name)
    result["new_source_negative_controls"] = mutations
    result["input_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["result"] = "PASS"
    if args.receipt:
        args.receipt.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

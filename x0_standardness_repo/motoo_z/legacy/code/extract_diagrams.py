#!/usr/bin/env python3
"""Extract the two framed dotted link diagrams from Tange's vector PDF.

This is a diagram transcription, not a proof that either handlebody closes.
Bezier control points and every artificial continuation across a drawing gap
are retained so that the transcription can be independently checked.

PD convention: zero-based arc labels, four incident arcs in counterclockwise
order (mathematical, y-up coordinates), starting with the incoming underarc.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import fitz

SOURCE_SHA256 = "493c2637a41dc6144a45d3a10361bb1d077da96775d973e287330edc5a378573"
DECORATIONS = [
    ("top_zero_meridian", 2, 0),
    ("left_dotted", 1, None),
    ("upper_dotted", 1, None),
    ("minus_one", 2, -1),
    ("bottom_dotted", 1, None),
    ("lower_inner_zero", 2, 0),
    ("lower_outer_zero", 2, 0),
    ("central_zero", 2, 0),
]
# The actual component cycles have been visually checked after automatic gap
# matching. These checks catch any upstream change in vector ordering or
# accidental change in the chosen gap matching.
EXPECTED_CYCLES_1 = [
    [0, 1], [2, 9, 8, 7, 6, 5, 4, 3],
    [10, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11],
    [24, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25],
    [41, 45, 44, 43, 42], [46, 50, 49, 48, 47], [51, 55, 54, 53, 52],
    [56, 62, 61, 60, 59, 58, 57],
]
EXPECTED_CYCLES_2 = EXPECTED_CYCLES_1[:-1] + [[56, 60, 59, 58, 57]]


def vadd(a, b): return (a[0] + b[0], a[1] + b[1])
def vsub(a, b): return (a[0] - b[0], a[1] - b[1])
def vmul(a, t): return (a[0] * t, a[1] * t)
def dot(a, b): return a[0] * b[0] + a[1] * b[1]
def cross(a, b): return a[0] * b[1] - a[1] * b[0]
def norm(a): return math.hypot(*a)
def unit(a): return vmul(a, 1.0 / norm(a))
def midpoint(a, b): return vmul(vadd(a, b), 0.5)


def cubic_flatten(p, tolerance, t0=0.0, t1=1.0, depth=0):
    """Polyline with source-parameter values and a control-polygon error bound."""
    chord = vsub(p[3], p[0])
    length = norm(chord)
    if length:
        err = max(abs(cross(chord, vsub(q, p[0]))) / length for q in p[1:3])
    else:
        err = max(norm(vsub(q, p[0])) for q in p[1:])
    if err <= tolerance or depth >= 26:
        return [(p[0], t0), (p[3], t1)]
    a, b, c = midpoint(p[0], p[1]), midpoint(p[1], p[2]), midpoint(p[2], p[3])
    d, e = midpoint(a, b), midpoint(b, c)
    f = midpoint(d, e)
    tm = (t0 + t1) / 2
    left = cubic_flatten([p[0], a, d, f], tolerance, t0, tm, depth + 1)
    right = cubic_flatten([f, e, c, p[3]], tolerance, tm, t1, depth + 1)
    return left[:-1] + right


def segment_intersection(p, q, r, s, include_ends=False):
    a, b = vsub(q, p), vsub(s, r)
    den = cross(a, b)
    if abs(den) <= 1e-14:
        return None
    t, u = cross(vsub(r, p), b) / den, cross(vsub(r, p), a) / den
    eps = 1e-10
    if (eps < t < 1-eps and -eps <= u <= 1+eps) if not include_ends else (-eps <= t <= 1+eps and -eps <= u <= 1+eps):
        return t, u, vadd(p, vmul(a, t))
    return None


def read_strokes(page, figure):
    drawings = page.get_drawings()
    selected = [(i, d) for i, d in enumerate(drawings)
                if d["type"] == "s" and ((d["rect"].y0 < 350) == (figure == 1))]
    result = []
    for local, (pdf_id, d) in enumerate(selected):
        items = []
        for item in d["items"]:
            if item[0] not in ("l", "c"):
                raise ValueError(f"Unexpected vector command: {item[0]}")
            items.append({"type": item[0], "points": [[float(p.x), float(p.y)] for p in item[1:]]})
        result.append({"id": local, "pdf_drawing_index": pdf_id, "items": items})
    return result


def endpoints(strokes):
    result = {}
    for s in strokes:
        first, last = s["items"][0]["points"], s["items"][-1]["points"]
        result[s["id"], 0] = (first[0], unit(vsub(first[0], first[1])))
        result[s["id"], 1] = (last[-1], unit(vsub(last[-1], last[-2])))
    return result


def match_gaps(strokes):
    ep = endpoints(strokes)
    choices = []
    keys = sorted(ep)
    for ii, x in enumerate(keys):
        p, v = ep[x]
        for y in keys[ii+1:]:
            if x[0] == y[0]:
                continue
            q, w = ep[y]
            distance = norm(vsub(q, p))
            if distance > 15:
                continue
            direction = unit(vsub(q, p))
            score = distance + 40 * (1-dot(direction, v)) + 40 * (1+dot(direction, w))
            choices.append((score, x, y, distance))
    used, mates, selected = set(), {}, []
    for score, x, y, distance in sorted(choices):
        if x in used or y in used:
            continue
        used.update([x, y])
        mates[x], mates[y] = y, x
        selected.append({"endpoints": [list(x), list(y)], "distance": distance,
                         "tangent_pairing_score": score})
    if len(used) != len(ep):
        raise ValueError("Some vector endpoints were not paired")
    if max(g["distance"] for g in selected) >= 10:
        raise ValueError("Unusually long drawing-gap continuation")
    return ep, mates, selected


def make_cycles(strokes, mates):
    visited, cycles, lookup = set(), [], {}
    for i in range(len(strokes)):
        if i in visited:
            continue
        current, start = (i, 0), (i, 0)
        cycle = []
        while current[0] not in visited:
            k, end = current
            visited.add(k)
            cycle.append({"stroke": k, "direction": 1 if end == 0 else -1})
            lookup[k] = (len(cycles), len(cycle)-1, 1 if end == 0 else -1)
            current = mates[k, 1-end]
        if current != start:
            raise ValueError("Gap matching does not give disjoint oriented circles")
        cycles.append(cycle)
    return cycles, lookup


def flatten_strokes(strokes, tolerance):
    result = []
    for s in strokes:
        for item_i, item in enumerate(s["items"]):
            pts = item["points"]
            poly = [(pts[0], 0.0), (pts[1], 1.0)] if item["type"] == "l" else cubic_flatten(pts, tolerance)
            for (a, ta), (b, tb) in zip(poly, poly[1:]):
                result.append({"stroke": s["id"], "item": item_i, "a": a, "b": b,
                               "t0": ta, "t1": tb})
    return result


def extract_crossings(strokes, ep, mates, cycles, lookup, tolerance):
    segments = flatten_strokes(strokes, tolerance)
    crossings, gap_records = [], []
    # Parameterisation is along each oriented stroke followed by its gap.
    for ci, cycle in enumerate(cycles):
        for pos, entry in enumerate(cycle):
            si, ori = entry["stroke"], entry["direction"]
            end = (si, 1 if ori == 1 else 0)
            nxt = mates[end]
            p, q = ep[end][0], ep[nxt][0]
            hits = []
            for seg in segments:
                if seg["stroke"] in (si, nxt[0]):
                    # Continuations cannot cross either adjacent stroke in these
                    # drawings; a separate complete intersection audit checks it.
                    continue
                inter = segment_intersection(p, q, seg["a"], seg["b"])
                if inter is None:
                    continue
                t, u, point = inter
                if any(norm(vsub(point, h["position_pdf"])) < 1e-5 for h in hits):
                    continue
                oi, opos, oori = lookup[seg["stroke"]]
                parameter = seg["item"] + seg["t0"] + u * (seg["t1"]-seg["t0"])
                if oori == -1:
                    parameter = len(strokes[seg["stroke"]]["items"]) - parameter
                over_vector = vmul(vsub(seg["b"], seg["a"]), oori)
                under_vector = vsub(q, p)
                # In y-up coordinates the determinant has the opposite sign.
                det_y_up = -cross(over_vector, under_vector)
                sine = abs(cross(unit(over_vector), unit(under_vector)))
                if sine < 0.05:
                    raise ValueError("Nearly tangential crossing requires manual review")
                hits.append({
                    "position_pdf": list(point),
                    "under_component": ci, "over_component": oi,
                    "under_order": [pos, len(strokes[si]["items"])+t],
                    "over_order": [opos, parameter],
                    "under_tangent_pdf": list(unit(under_vector)),
                    "over_tangent_pdf": list(unit(over_vector)),
                    "sign": 1 if det_y_up > 0 else -1,
                    "crossing_angle_sine": sine,
                    "under_gap": [list(end), list(nxt)],
                    "over_stroke": seg["stroke"],
                    "over_item": seg["item"],
                })
            # Some drawing gaps span two nearby overstrands. Retain both in
            # their geometric order instead of imposing one crossing per gap.
            gap_records.append({"from": list(end), "to": list(nxt), "component": ci,
                                "start_pdf": list(p), "end_pdf": list(q),
                                "number_of_crossings": len(hits)})
            crossings.extend(hits)
    crossings.sort(key=lambda c: (c["position_pdf"][1], c["position_pdf"][0]))
    for i, c in enumerate(crossings):
        c["id"] = i
    return crossings, gap_records, segments


def make_pd(crossings, number_components):
    visits = [[] for _ in range(number_components)]
    for c in crossings:
        visits[c["under_component"]].append((c["under_order"], c["id"], "under"))
        visits[c["over_component"]].append((c["over_order"], c["id"], "over"))
    arcs, component_arcs, ordered_visits = {}, [], []
    next_label = 0
    for comp_visits in visits:
        comp_visits.sort()
        if not comp_visits:
            raise ValueError("Crossing-free component needs explicit encoding")
        labels = list(range(next_label, next_label+len(comp_visits)))
        next_label += len(comp_visits)
        for j, (_, cid, role) in enumerate(comp_visits):
            arcs[cid, role] = {"incoming": labels[(j-1) % len(labels)], "outgoing": labels[j]}
        component_arcs.append(labels)
        ordered_visits.append([{"crossing": cid, "role": role} for _, cid, role in comp_visits])
    pd = []
    for c in crossings:
        cid = c["id"]
        u, o = arcs[cid, "under"], arcs[cid, "over"]
        uv, ov = c["under_tangent_pdf"], c["over_tangent_pdf"]
        rays = [(u["incoming"], vmul(uv, -1)), (u["outgoing"], uv),
                (o["incoming"], vmul(ov, -1)), (o["outgoing"], ov)]
        base = math.atan2(-rays[0][1][1], rays[0][1][0])
        rays.sort(key=lambda z: (math.atan2(-z[1][1], z[1][0])-base) % (2*math.pi))
        labels = [z[0] for z in rays]
        # Avoid a possible atan2 wrap-round of the reference ray itself.
        at = labels.index(u["incoming"])
        labels = labels[at:] + labels[:at]
        if labels[2] != u["outgoing"]:
            raise ValueError("Crossing rays do not alternate transversely")
        c["under_arcs"], c["over_arcs"], c["pd_tuple"] = u, o, labels
        pd.append(labels)
    occurrences = [0] * next_label
    for row in pd:
        for label in row:
            occurrences[label] += 1
    if set(occurrences) != {2}:
        raise ValueError("A PD arc does not have exactly two ends")
    return pd, component_arcs, ordered_visits


def linking_matrix(crossings, n):
    twice = [[0]*n for _ in range(n)]
    writhe = [0]*n
    for c in crossings:
        i, j = c["under_component"], c["over_component"]
        if i == j:
            writhe[i] += c["sign"]
        else:
            twice[i][j] += c["sign"]
            twice[j][i] += c["sign"]
    if any(x % 2 for row in twice for x in row):
        raise ValueError("An intercomponent crossing-sign sum is odd")
    return [[x//2 for x in row] for row in twice], writhe


def audit_all_intersections(result):
    """Check that no unrecorded solid/solid or gap/gap crossing was ignored.

    The count is repeated from a spatially indexed union of all flattened
    solid strokes and all explicit gap continuations, including continuations
    adjacent to their own stroke (which the main extraction skips).
    """
    strokes = result["geometry"]["strokes"]
    segments = flatten_strokes(strokes, result["geometry"]["bezier_flattening_tolerance_points"])
    for s in segments:
        s["kind"] = "solid"
    for i, g in enumerate(result["geometry"]["oriented_gaps"]):
        segments.append({"kind": "gap", "gap_id": i,
                         "a": g["start_pdf"], "b": g["end_pdf"]})
    cells = {}
    side = 3.0
    for i, s in enumerate(segments):
        x0, x1 = sorted([s["a"][0], s["b"][0]])
        y0, y1 = sorted([s["a"][1], s["b"][1]])
        s["bbox"] = (x0, y0, x1, y1)
        for x in range(math.floor(x0/side), math.floor(x1/side)+1):
            for y in range(math.floor(y0/side), math.floor(y1/side)+1):
                cells.setdefault((x,y), []).append(i)
    seen_pairs, hits = set(), []
    for bucket in cells.values():
        for ii, i in enumerate(bucket):
            a = segments[i]
            ax0, ay0, ax1, ay1 = a["bbox"]
            for j in bucket[ii+1:]:
                b = segments[j]
                bx0, by0, bx1, by1 = b["bbox"]
                if ax0 > bx1 or bx0 > ax1 or ay0 > by1 or by0 > ay1:
                    continue
                key = (i,j) if i<j else (j,i)
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                hit = segment_intersection(a["a"], a["b"], b["a"], b["b"])
                if hit is None:
                    continue
                point = hit[2]
                if any(norm(vsub(point, old["point"])) < 1e-4 for old in hits):
                    continue
                hits.append({"kind": sorted([a["kind"], b["kind"]]), "point": point,
                             "segment_a": i, "segment_b": j})
    unexplained = [h for h in hits if h["kind"] != ["gap", "solid"]]
    if unexplained:
        raise ValueError(f"Ambiguous unbroken-stroke crossings: {unexplained}")
    if len(hits) != result["crossing_count"]:
        raise ValueError(f"Full intersection audit counted {len(hits)} instead of {result['crossing_count']}")
    max_error = 0.0
    for c in result["crossings"]:
        error = min(norm(vsub(c["position_pdf"], h["point"])) for h in hits)
        max_error = max(max_error, error)
        if error > 1e-4:
            raise ValueError("Recorded crossing absent from complete geometric audit")
    return {"figure": result["figure"], "full_intersections": len(hits),
            "solid_solid_crossings": 0, "gap_gap_crossings": 0,
            "all_crossings_explained_by_underpass_gaps": True,
            "maximum_crossing_coordinate_difference_points": max_error,
            "segment_count_including_gap_continuations": len(segments)}


def annotate(pdf, result, out):
    doc = fitz.open(pdf)
    p = doc[0]
    colours = [(0.6,0.2,0.8),(0.05,0.5,0.1),(0.1,0.2,0.85),(0.9,0.4,0),
               (0,0.65,0.65),(0.55,0.35,0.05),(0.65,0.1,0.4),(0.9,0,0)]
    for comp in result["components"]:
        for item in comp["oriented_stroke_cycle"]:
            s = result["geometry"]["strokes"][item["stroke"]]
            shape = p.new_shape()
            for seg in s["items"]:
                pts = [fitz.Point(*q) for q in seg["points"]]
                if seg["type"] == "c": shape.draw_bezier(*pts)
                else: shape.draw_line(*pts)
            shape.finish(color=colours[comp["index"]], width=1.3, closePath=False)
            shape.commit()
    for c in result["crossings"]:
        x,y = c["position_pdf"]
        # Small labels are intended for large-resolution inspection.
        p.insert_text((x+1.4,y-1.4), str(c["id"]), fontsize=3.4, color=(0,0,0), overlay=True)
    clip = fitz.Rect(195,104,410,335) if result["figure"] == 1 else fitz.Rect(195,356,410,586)
    p.get_pixmap(matrix=fitz.Matrix(7,7), clip=clip).save(out)


def extract(pdf, figure, tolerance):
    page = fitz.open(pdf)[0]
    strokes = read_strokes(page, figure)
    ep, mates, pairing = match_gaps(strokes)
    cycles, lookup = make_cycles(strokes, mates)
    expected = EXPECTED_CYCLES_1 if figure == 1 else EXPECTED_CYCLES_2
    if [[x["stroke"] for x in c] for c in cycles] != expected:
        raise ValueError("Unexpected component cycles; inspect the vector matching")
    crossings, gaps, segments = extract_crossings(strokes, ep, mates, cycles, lookup, tolerance)
    pd, component_arcs, visits = make_pd(crossings, len(cycles))
    linking, writhe = linking_matrix(crossings, len(cycles))
    components = []
    for i, (name, handle_index, framing) in enumerate(DECORATIONS):
        components.append({"index": i, "id": f"C{i}", "name": name,
                           "handle_index": handle_index, "dotted": handle_index == 1,
                           "framing": framing, "arc_labels": component_arcs[i],
                           "oriented_stroke_cycle": cycles[i], "crossing_visits": visits[i],
                           "self_writhe": writhe[i]})
    empty = [g for g in gaps if g["number_of_crossings"] == 0]
    return {
        "schema": "motoo.vector_link.v1", "figure": figure,
        "source_sha256": SOURCE_SHA256, "source_page": 1,
        "status": "TRANSCRIBED_DIAGRAM_REQUIRES_INDEPENDENT_AUDIT",
        "pd_convention": "Zero-based arcs; counterclockwise in mathematical y-up coordinates, beginning at the incoming underarc.",
        "sign_convention": "sign(det(over_tangent, under_tangent)) in mathematical y-up coordinates",
        "crossing_count": len(crossings), "component_count": len(components),
        "pd_code": pd, "components": components, "crossings": crossings,
        "linking_matrix": linking,
        "geometry": {"coordinate_system": "PDF page points, x right and y down",
                     "strokes": strokes, "gap_pairing": pairing,
                     "oriented_gaps": gaps, "empty_gaps": empty,
                     "bezier_flattening_tolerance_points": tolerance,
                     "flattened_segment_count": len(segments)},
        "qualification": "Framing labels are interpreted as the displayed integral Kirby framings. This encoding alone does not certify dotted-unlink status, handlebody closure, or a homotopy sphere."
    }


def validate_with_libraries(result, link_python):
    """Independently import the PD and match each crossing's marked components."""
    script = r'''
import importlib.metadata, json, sys
import regina
from spherogram import Link
d = json.load(sys.stdin)
link = Link(d['pd_code'], check_planarity=True)
mapping = {}
for i, c in enumerate(link.crossings):
    if c.label != i:
        raise ValueError('Spherogram crossing numbering differs from input')
    expected = d['crossings'][i]
    if c.sign != expected['sign']:
        raise ValueError('Crossing orientation/sign mismatch')
    for imported, source in [(c.strand_components[0], expected['under_component']),
                             (c.strand_components[1], expected['over_component'])]:
        if imported in mapping and mapping[imported] != source:
            raise ValueError('Inconsistent source-to-PD component mapping')
        mapping[imported] = source
order = [mapping[i] for i in range(len(link.link_components))]
if sorted(order) != list(range(d['component_count'])):
    raise ValueError('Component mapping is not a bijection')
matrix = link.linking_matrix()
for i, row in enumerate(matrix):
    for j, x in enumerate(row):
        if x != d['linking_matrix'][order[i]][order[j]]:
            raise ValueError('Pairwise linking number mismatch')
rpd = [[a+1 for a in row] for row in d['pd_code']]
rlink = regina.Link.fromPD(rpd)
if rlink.size() != d['crossing_count'] or rlink.countComponents() != d['component_count']:
    raise ValueError('Regina PD import has wrong size')
if any(rlink.crossing(i).sign() != c['sign'] for i,c in enumerate(d['crossings'])):
    raise ValueError('Regina PD crossing sign mismatch')
print(json.dumps({'pd_component_to_source': order,
                  'decorations': [('x' if d['components'][i]['dotted'] else d['components'][i]['framing']) for i in order],
                  'regina_pd_code': rpd,
                  'pd_validation': {'spherogram_version': importlib.metadata.version('spherogram'),
                                    'regina_distribution_version': importlib.metadata.version('regina'),
                                    'planarity_check': True,
                                    'all_crossing_signs_match': True,
                                    'all_marked_components_match': True,
                                    'linking_matrix_matches': True,
                                    'regina_import_and_crossing_signs_match': True}}))
'''
    p = subprocess.run([str(link_python), "-c", script], input=json.dumps(result),
                       text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError(f"PD library validation failed:\n{p.stderr}")
    result.update(json.loads(p.stdout))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--pdf", type=Path, default=root/"inputs/motoo_rollspun.pdf")
    parser.add_argument("--out", type=Path, default=root/"inputs")
    parser.add_argument("--scratch", type=Path, default=root/"scratch/diagram")
    parser.add_argument("--tolerance", type=float, default=1e-4)
    parser.add_argument("--link-python", type=Path, default=root/"venv/bin/python")
    args = parser.parse_args()
    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise SystemExit("Pinned source PDF checksum mismatch")
    args.out.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    audits = []
    for figure in (1,2):
        result = extract(args.pdf, figure, args.tolerance)
        finer = extract(args.pdf, figure, args.tolerance/2)
        if result["pd_code"] != finer["pd_code"]:
            raise ValueError("PD transcription changed at finer curve resolution")
        result["geometry"]["finer_tolerance_same_pd"] = args.tolerance/2
        validate_with_libraries(result, args.link_python)
        audits.append(audit_all_intersections(result))
        output = args.out/f"figure{figure}_link.json"
        output.write_text(json.dumps(result, indent=2)+"\n")
        annotate(args.pdf, result, args.scratch/f"figure{figure}_annotated.png")
        print(json.dumps({"figure": figure, "components": result["component_count"],
                          "crossings": result["crossing_count"], "empty_gaps": len(result["geometry"]["empty_gaps"]),
                          "file": str(output), "sha256": hashlib.sha256(output.read_bytes()).hexdigest()}))
    (args.scratch/"complete_intersection_audit.json").write_text(json.dumps(audits,indent=2)+"\n")


if __name__ == "__main__":
    main()

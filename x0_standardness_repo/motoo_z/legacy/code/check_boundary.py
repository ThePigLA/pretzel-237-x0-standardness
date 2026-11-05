#!/usr/bin/env python3
"""Recognise the 3D boundary of a dotted Kirby 0/1/2-handlebody.

This is a boundary calculation, not a construction of a closed 4-manifold.
Dots are replaced by zero-framed surgery; integer 2-handle framings are
filled in the meridian / preferred-longitude basis of Link.exterior().
The permanent filling is imported from its exact tetrahedron face data.
Regina's summands() supplies the native prime-decomposition trust boundary.
All component-to-cusp assignments are checked against the input PD arcs.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time

import regina
import snappy


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def versions():
    return {
        "python": sys.version,
        "regina_distribution": importlib.metadata.version("regina"),
        "regina_engine": regina.versionString(),
        "snappy": snappy.__version__,
        "spherogram": importlib.metadata.version("spherogram"),
    }


def import_link(document):
    """Return Link and checked cusp metadata, without assuming source order."""
    pd = document["pd_code"]
    require(isinstance(pd, list) and pd, "PD must have crossings")
    require(all(isinstance(row, list) and len(row) == 4 for row in pd),
            "every PD crossing must have four entries")
    require(all(type(a) is int and a >= 0 for row in pd for a in row),
            "PD arc labels must be nonnegative integers")
    counts = collections.Counter(a for row in pd for a in row)
    require(set(counts.values()) == {2}, "each PD arc must occur twice")
    components = document["components"]
    source_sets = [frozenset(c["arc_labels"]) for c in components]
    require(all(source_sets), "every component must have source arcs")
    require(len(set(source_sets)) == len(source_sets), "duplicate source component")
    require(sum(map(len, source_sets)) == len(set().union(*source_sets)),
            "source components overlap")
    require(set().union(*source_sets) == set(counts), "component arcs do not cover PD")
    link = snappy.Link(pd, check_planarity=True)
    require(len(link.crossings) == len(pd), "PD import changed crossing count")
    require(len(link.link_components) == len(components), "PD import changed component count")
    source_by_set = {s: i for i, s in enumerate(source_sets)}
    cusp_data = []
    for cusp, component in enumerate(link.link_components):
        arcs = set()
        visits = []
        for ep in component:
            crossing = int(ep.crossing.label)
            require(0 <= crossing < len(pd), "crossing identity was not preserved")
            # Spherogram may rotate a crossing by 180 degrees during orientation.
            # The opposite pair of input arcs is invariant under that rotation.
            parity = ep.strand_index % 2
            pair = [pd[crossing][parity], pd[crossing][parity + 2]]
            arcs.update(pair)
            visits.append([crossing, parity, pair])
        require(frozenset(arcs) in source_by_set,
                f"cusp {cusp} has no matching source component")
        source = source_by_set[frozenset(arcs)]
        c = components[source]
        handle = c["handle_index"]
        require(handle in (1, 2), "only 1- and 2-handle decorations are allowed")
        require(bool(c.get("dotted", handle == 1)) == (handle == 1),
                "inconsistent dotted/handle_index metadata")
        if handle == 1:
            require(c["framing"] is None, "dotted handle must have null framing")
            framing = 0
        else:
            require(type(c["framing"]) is int, "2-handle framing must be integer")
            framing = c["framing"]
        cusp_data.append({
            "cusp_index": cusp, "source_index": source,
            "source_id": c.get("id", str(source)),
            "source_name": c.get("name"), "handle_index": handle,
            "original_framing": c["framing"],
            "boundary_surgery_slope": [framing, 1],
            "matched_source_arc_labels": sorted(arcs),
            "crossing_pair_checks": visits,
        })
    require(sorted(x["source_index"] for x in cusp_data) == list(range(len(components))),
            "source-to-cusp assignment is not bijective")
    mapping = [x["source_index"] for x in cusp_data]
    if "pd_component_to_source" in document:
        supplied = document["pd_component_to_source"]
        if supplied and isinstance(supplied[0], str):
            require(supplied == [x["source_id"] for x in cusp_data],
                    "supplied component map disagrees with source arcs")
        else:
            require(supplied == mapping, "supplied component map disagrees with source arcs")
    if "decorations" in document:
        expected = ["x" if x["handle_index"] == 1 else x["original_framing"] for x in cusp_data]
        require(document["decorations"] == expected,
                "supplied decoration order disagrees with independently matched cusps")
    return link, cusp_data


def regina_from_exact_gluings(gluings):
    """Check the finite closed tetrahedral gluing table and import it directly."""
    count = len(gluings)
    require(count > 0, "empty closed filling")
    for i, (neighbours, permutations) in enumerate(gluings):
        require(len(neighbours) == len(permutations) == 4, "bad tetrahedron row")
        for f, (j, p) in enumerate(zip(neighbours, permutations)):
            require(type(j) is int and 0 <= j < count, "unfilled face in closed output")
            require(sorted(p) == list(range(4)), "face map is not a permutation")
            g = p[f]
            require((i, f) != (j, g), "face glued to itself")
            other_neighbours, other_permutations = gluings[j]
            require(other_neighbours[g] == i, "non-reciprocal neighbour")
            require(all(other_permutations[g][p[v]] == v for v in range(4)),
                    "face maps are not reciprocal inverses")
    tri = regina.Triangulation3()
    tetrahedra = tri.newTetrahedra(count)
    for i, (neighbours, permutations) in enumerate(gluings):
        for f, (j, p) in enumerate(zip(neighbours, permutations)):
            if (i, f) < (j, p[f]):
                tetrahedra[i].join(f, tetrahedra[j], regina.Perm4(*p))
    require(tri.isValid(), "permanent filling is not a valid triangulation")
    require(tri.isClosed(), "permanent filling is not closed")
    require(tri.isConnected(), "permanent filling is disconnected")
    require(tri.isOrientable(), "boundary filling is non-orientable")
    return tri


def save_regina(tri, path):
    packet = regina.PacketOfTriangulation3(tri)
    require(packet.save(str(path)), f"failed to save {path}")


def describe_tri(tri):
    return {
        "tetrahedra": tri.size(),
        "f_vector": [tri.countVertices(), tri.countEdges(), tri.countTriangles(), tri.size()],
        "isosig": tri.isoSig(), "homology_H1": str(tri.homology()),
        "valid": tri.isValid(), "closed": tri.isClosed(),
        "connected": tri.isConnected(), "orientable": tri.isOrientable(),
    }


def run_check(input_path, outdir, expected):
    started = time.monotonic()
    input_path, outdir = Path(input_path), Path(outdir)
    require(not outdir.exists(), f"output directory already exists: {outdir}")
    outdir.mkdir(parents=True)
    raw = input_path.read_bytes()
    document = json.loads(raw)
    (outdir / "input_link.json").write_bytes(raw)
    receipt = {
        "schema": "motoo.kirby.boundary-check.v1",
        "status": "RUNNING", "input_file_sha256": hashlib.sha256(raw).hexdigest(),
        "source_pdf_sha256": document.get("source_sha256"),
        "checker_sha256": sha256(__file__), "versions": versions(),
        "expected_s2xs1_summands": expected,
        "scope": "3-dimensional boundary of the 0/1/2-handlebody only",
        "explicit_4d_upper_handle_gluing_constructed": False,
        "trust_boundary": [
            "Spherogram planar link import and meridian/preferred-longitude exterior construction",
            "SnapPea exact combinatorial permanent Dehn filling",
            "Regina prime decomposition by normal-surface crushing and sphere recognition",
        ],
    }

    def stage(name):
        receipt["stage"] = name
        receipt["elapsed_seconds"] = time.monotonic() - started
        write_json(outdir / "progress.json", receipt)
        print(json.dumps({"input": input_path.name, "stage": name,
                          "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)

    try:
        link, cusps = import_link(document)
        receipt["component_mapping"] = cusps
        receipt["crossings"] = len(link.crossings)
        receipt["components"] = len(link.link_components)
        receipt["spherogram_PD_after_orientation"] = link.PD_code()
        stage("PD_COMPONENT_ASSIGNMENTS_CHECKED")
        exterior = link.exterior(with_hyperbolic_structure=False)
        require(exterior.num_cusps() == len(cusps), "exterior cusp count disagrees")
        (outdir / "link_exterior.tri").write_text(exterior._to_string())
        receipt["exterior_tetrahedra"] = exterior.num_tetrahedra()
        receipt["exterior_H1"] = str(exterior.homology())
        slopes = [x["boundary_surgery_slope"] for x in cusps]
        exterior.dehn_fill(slopes)
        receipt["filled_boundary_H1_snappy"] = str(exterior.homology())
        (outdir / "boundary_with_filling_coefficients.tri").write_text(exterior._to_string())
        stage("SURGERY_SLOPES_APPLIED")
        filled = exterior.filled_triangulation()
        require(filled.num_cusps() == 0, "permanent filling left cusps")
        (outdir / "boundary_closed.tri").write_text(filled._to_string())
        gluings = filled._get_tetrahedra_gluing_data()
        write_json(outdir / "boundary_closed_gluings.json", {
            "dimension": 3, "tetrahedra": len(gluings),
            "gluings": gluings, "source_input_sha256": receipt["input_file_sha256"],
            "surgery_slopes": slopes,
        })
        tri = regina_from_exact_gluings(gluings)
        other = regina.Triangulation3.fromSnapPea(filled._to_string())
        require(tri.isoSig() == other.isoSig(), "two exact import routes disagree")
        receipt["independent_face_table_import_agrees_with_snappea_import"] = True
        receipt["closed_boundary"] = describe_tri(tri)
        save_regina(tri, outdir / "boundary_closed.rga")
        stage("CLOSED_VALID_BOUNDARY_IMPORTED")
        before = tri.size()
        tri.simplify()
        receipt["native_simplification"] = {"tetrahedra_before": before,
                                             "tetrahedra_after": tri.size()}
        receipt["simplified_boundary"] = describe_tri(tri)
        save_regina(tri, outdir / "boundary_simplified.rga")
        stage("NATIVE_BOUNDARY_SIMPLIFICATION_FINISHED")
        summands = tri.summands()
        reference = regina.Example3.s2xs1()
        receipt["s2xs1_reference"] = describe_tri(reference)
        summaries = []
        all_s2xs1 = True
        for index, summand in enumerate(summands):
            info = describe_tri(summand)
            iso = summand.isIsomorphicTo(reference)
            info["isomorphic_to_standard_s2xs1"] = iso is not None
            if iso is not None:
                info["explicit_isomorphism"] = [
                    {"simplex": i, "image": iso.simpImage(i),
                     "vertex_permutation": [iso.facetPerm(i)[j] for j in range(4)]}
                    for i in range(summand.size())
                ]
            else:
                all_s2xs1 = False
            standard = regina.StandardTriangulation.recognise(summand)
            info["regina_standard_name"] = None if standard is None else str(standard.name())
            save_regina(summand, outdir / f"summand_{index:02}.rga")
            summaries.append(info)
        receipt["prime_summands"] = summaries
        receipt["recognized_s2xs1_summands"] = len(summands) if all_s2xs1 else None
        receipt["boundary_kind"] = (
            "S3" if all_s2xs1 and not summands else
            "connected_sum_of_S2xS1" if all_s2xs1 else
            "other_closed_orientable_3_manifold"
        )
        if expected is not None:
            receipt["expected_boundary_verified"] = all_s2xs1 and len(summands) == expected
            receipt["status"] = (
                "PASS_EXPECTED_BOUNDARY" if receipt["expected_boundary_verified"]
                else "FAIL_EXPECTED_BOUNDARY"
            )
            if receipt["expected_boundary_verified"]:
                h1 = sum(c["handle_index"] == 1 for c in document["components"])
                h2 = len(document["components"]) - h1
                receipt["conditional_closed_handle_counts"] = [1, h1, h2, expected, 1]
                receipt["conditional_closed_euler_characteristic"] = 2 - h1 + h2 - expected
                receipt["closure_statement"] = (
                    "Boundary recognition establishes existence and diffeomorphism uniqueness "
                    "of the usual 3/4-handle completion. This program does not construct its "
                    "4-dimensional attaching triangulations or prove its standardness."
                )
        else:
            receipt["status"] = "CLASSIFIED_BOUNDARY"
        stage("PRIME_DECOMPOSITION_CHECKED")
        receipt["artifacts"] = {
            p.name: {"sha256": sha256(p), "bytes": p.stat().st_size}
            for p in sorted(outdir.iterdir()) if p.is_file() and p.name != "progress.json"
        }
        write_json(outdir / "BOUNDARY_VERIFICATION.json", receipt)
        print(json.dumps({"status": receipt["status"], "out": str(outdir),
                          "boundary_kind": receipt["boundary_kind"],
                          "recognized_s2xs1_summands": receipt["recognized_s2xs1_summands"]}), flush=True)
        return receipt
    except Exception as exc:
        receipt["status"] = "ERROR"
        receipt["error"] = f"{type(exc).__name__}: {exc}"
        receipt["elapsed_seconds"] = time.monotonic() - started
        write_json(outdir / "BOUNDARY_ERROR.json", receipt)
        raise


def fixture(name, framings):
    pd = [[2*i, 2*i+1, 2*i+1, 2*i] for i in range(len(framings))]
    return {
        "schema": "motoo.boundary.fixture.v1", "name": name, "pd_code": pd,
        "components": [
            {"id": f"C{i}", "name": "unknot", "index": i,
             "handle_index": 2, "dotted": False, "framing": framing,
             "arc_labels": [2*i, 2*i+1]}
            for i, framing in enumerate(framings)
        ],
    }


def self_test(outdir):
    outdir = Path(outdir)
    require(not outdir.exists(), "self-test output already exists")
    outdir.mkdir(parents=True)
    results = []
    for name, framings, expected in [
        ("zero_unknot", [0], 1),
        ("plus_one_unknot", [1], 0),
        ("two_zero_unlinks", [0, 0], 2),
        ("plus_three_unknot_negative", [3], None),
    ]:
        source = outdir / f"{name}.json"
        write_json(source, fixture(name, framings))
        receipt = run_check(source, outdir / name, expected)
        if expected is not None:
            require(receipt["expected_boundary_verified"], f"fixture failed: {name}")
        else:
            require(receipt["recognized_s2xs1_summands"] is None,
                    "lens-space negative control was incorrectly accepted")
        results.append({"name": name, "status": "PASS", "result": receipt["boundary_kind"],
                        "verification_sha256": sha256(outdir / name / "BOUNDARY_VERIFICATION.json")})
    report = {"status": "PASS", "versions": versions(), "checks": results,
              "checker_sha256": sha256(__file__)}
    write_json(outdir / "SELF_TEST.json", report)
    print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", help="source figure*_link.json")
    parser.add_argument("--out", required=True, help="new output directory")
    parser.add_argument("--expected", type=int, default=2,
                        help="expected number of S2xS1 prime summands (default 2)")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(args.out)
        return 0
    require(args.input is not None, "input JSON is required")
    require(args.expected >= 0, "expected count must be nonnegative")
    report = run_check(args.input, args.out, args.expected)
    return 0 if report["status"] == "PASS_EXPECTED_BOUNDARY" else 2


if __name__ == "__main__":
    raise SystemExit(main())

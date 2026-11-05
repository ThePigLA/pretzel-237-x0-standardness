#!/usr/bin/env python3
"""Build or replay a genus-two upper-handle completion control.

The lower Kirby diagram is the 0-framed two-component unlink, displayed
with a connected Reidemeister-II projection. Katie, ideal truncation and
the initial finite-body simplification are native Regina trust boundaries.
Later 3D Pachner moves, actual 4D layerings, reverse shellings and boundary
gluings have explicit replay records. This script does not infer standardness
of an unknown manifold from homology or a trivial fundamental group.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

import regina

from boundary_pachner_path import replay as replay_boundary, search
from complete_upper_handles import (boundary_data, cap_audit, glue_pieces,
    layer_path, require, save3, save4, save_json)


PROJECT = Path(__file__).resolve().parents[1]
CAP_RELATIVE = "results/boundary/standard_cap_six_boundary_v1/cap_with_matched_boundary.rga"
CAP_SHA256 = "2b46b7e8d760ce52641f78cb78c998748a8e476410b9c6525ae1b28168ca31b7"
GOAL = "gLALQccbefefdpcghfc"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load4(path):
    packet = regina.open(str(path))
    require(isinstance(packet, regina.Triangulation4), "expected a 4D RGA file")
    return regina.Triangulation4(packet)


def verify_source(data):
    require(data["decorations"] == [0, 0], "control framings differ")
    proof = data["unlink_proof"]
    link = regina.Link.tightDecoding(proof["initial_tight_encoding"])
    require(link.pdData() == [[x + 1 for x in c] for c in data["pd_code"]],
            "unlink certificate PD export differs from Katie input")
    require(link.countComponents() == 2 and link.virtualGenus() == 0,
            "control is not a classical two-component diagram")
    require(link.tightEncoding() == proof["initial_tight_encoding"], "unlink source differs")
    require(link.r2(link.crossing(proof["r2_crossing"])), "control R2 guard failed")
    require(link.size() == 0 and link.countComponents() == 2,
            "control did not reduce to two crossing-free unknots")
    require(link.tightEncoding() == proof["final_tight_encoding"], "unlink endpoint differs")
    return {"status": "PASS", "guarded_R2_moves": 1, "crossing_free_components": 2}


def audit_ideal(tri):
    require(tri.isValid() and tri.isConnected() and tri.isOrientable() and tri.isIdeal(),
            "raw control must be one valid orientable ideal 4D body")
    require(tri.countBoundaryComponents() == 1, "raw control must have one ideal boundary")
    rows = []
    exceptional = 0
    for vertex in tri.vertices():
        link = vertex.buildLink()
        sphere = link.isSphere()
        row = {"vertex": vertex.index(), "tetrahedra": link.size(), "is_sphere": bool(sphere)}
        if not sphere:
            summands = link.summands()
            require(link.isValid() and link.isClosed() and len(summands) == 2 and
                    all(s.isIsomorphicTo(regina.Example3.s2xs1()) is not None for s in summands),
                    "nonspherical raw link is not #2(S1xS2)")
            row["summands"] = [s.isoSig() for s in summands]
            exceptional += 1
        rows.append(row)
    require(exceptional == 1, "expected exactly one nonspherical ideal vertex")
    return {"pentachora": tri.size(), "components": tri.countComponents(),
            "nonspherical_vertices": exceptional, "vertex_links": rows}


def replay_layerings(lower, path, events_file, target):
    """Replay the actual 4D face gluings, independently of the producer's layering routine."""
    current = regina.Triangulation4(lower)
    events = [json.loads(line) for line in gzip.open(events_file, "rt")]
    require(len(events) == len(path["moves"]), "layering event count differs")
    for i, (event, move) in enumerate(zip(events, path["moves"])):
        require(event["step"] == i and event["added_pentachoron"] == current.size(),
                "wrong layering simplex index")
        for key, source_key in [("before", "before"), ("after", "after"),
                                ("dimension", "dimension"), ("face", "face")]:
            require(event[key] == move[source_key], "3D/4D record binding differs")
        simplex = current.newSimplex()
        for gluing in event["gluings"]:
            old_index, old_facet = gluing["old_pentachoron"], gluing["old_facet"]
            require(0 <= old_index < simplex.index(), "layering target is not an old simplex")
            old = current.simplex(old_index)
            require(old.adjacentSimplex(old_facet) is None, "layering target is not boundary")
            permutation = regina.Perm5(*gluing["permutation"])
            require(permutation[gluing["new_facet"]] == old_facet, "wrong opposite-facet map")
            simplex.join(gluing["new_facet"], old, permutation)
        require(current.hasShellBoundary(simplex), "layering lacks guarded inverse shelling")
    require(current == target, "explicit layering replay differs from saved labelled piece")
    require(boundary_data(current)[0].isoSig() == path["end"], "wrong layered boundary")
    for i in reversed(range(lower.size(), current.size())):
        require(current.shellBoundary(current.simplex(i)), "reverse shelling failed")
    require(current == lower, "reverse shelling did not recover exact labelled lower body")
    return {"status": "PASS", "moves": len(events), "all_inverse_shelling_guards": True,
            "reverse_shellings_recover_exact_source": True}


def replay_gluing(first, second, maps_file, saved):
    maps = json.loads(Path(maps_file).read_text())
    first_boundary, first_locations = boundary_data(first)
    second_boundary, second_locations = boundary_data(second)
    require(len(maps) == first_boundary.size() == second_boundary.size(),
            "boundary gluing does not pair every tetrahedron")
    correspondence = {}
    destinations = set()
    for row in maps:
        i, j = row["source_boundary_tetrahedron"], row["target_boundary_tetrahedron"]
        require(0 <= i < first_boundary.size() and i not in correspondence and
                0 <= j < second_boundary.size() and j not in destinations,
                "boundary tetrahedron correspondence is not bijective")
        a, b = first_locations[i], second_locations[j]
        require((row["source_pentachoron"], row["source_facet"]) ==
                (a["pentachoron"], a["facet"]) and
                (row["target_pentachoron"], row["target_facet"]) ==
                (b["pentachoron"], b["facet"]), "wrong boundary embedding in gluing record")
        p = row["permutation"]
        require(sorted(p) == list(range(5)) and p[a["facet"]] == b["facet"],
                "wrong 4D facet gluing permutation")
        q = [b["vertices"].index(p[a["vertices"][v]]) for v in range(4)]
        correspondence[i] = (j, q)
        destinations.add(j)
    for i, (j, q) in correspondence.items():
        source_tet, target_tet = first_boundary.simplex(i), second_boundary.simplex(j)
        for facet in range(4):
            other = source_tet.adjacentSimplex(facet)
            require(other is not None, "source boundary is not closed")
            image_other, q_other = correspondence[other.index()]
            require(target_tet.adjacentSimplex(q[facet]).index() == image_other,
                    "boundary isomorphism fails adjacency")
            old_map = source_tet.adjacentGluing(facet)
            new_map = target_tet.adjacentGluing(q[facet])
            require(all(new_map[q[v]] == q_other[old_map[v]] for v in range(4)),
                    "boundary isomorphism fails vertex-map compatibility")
    result = regina.Triangulation4(first)
    offset = result.size()
    result.insertTriangulation(second)
    for row in maps:
        source = result.simplex(row["source_pentachoron"])
        target = result.simplex(offset + row["target_pentachoron"])
        require(source.adjacentSimplex(row["source_facet"]) is None,
                "gluing source is not a boundary facet")
        source.join(row["source_facet"], target, regina.Perm5(*row["permutation"]))
    require(result == saved, "boundary gluing replay differs from saved labelled triangulation")
    require(result.isValid() and result.isConnected() and result.isClosed() and
            result.isOrientable() and not result.isIdeal(), "closed control fails manifold flags")
    links = []
    for vertex in result.vertices():
        link = vertex.buildLink()
        print(f"Recognising closed-control vertex link {vertex.index()} ({link.size()} tetrahedra)", flush=True)
        require(link.isSphere(), "closed control has a nonspherical vertex link")
        links.append({"vertex": vertex.index(), "tetrahedra": link.size(), "is_sphere": True})
    require(result.homology().isTrivial() and result.homology(2).isTrivial() and
            result.eulerCharManifold() == 2, "closed control has wrong homology or Euler characteristic")
    print("Checking the closed-control fundamental group", flush=True)
    group = result.fundamentalGroup()
    group.simplify()
    require(group.countGenerators() == 0, "closed control group did not simplify to trivial")
    return {"status": "PASS", "pentachora": result.size(), "H1": str(result.homology()),
            "H2": str(result.homology(2)), "euler_characteristic": result.eulerCharManifold(),
            "boundary_isomorphism_independently_checked": True,
            "native_group_generators_after_simplify": 0, "vertex_links": links}


def audit(output):
    output = Path(output)
    print("Checking the source unlink and Katie input binding", flush=True)
    source_file = output / "connected_input.json"
    source = json.loads(source_file.read_text())
    source_check = verify_source(source)
    katie_receipt = json.loads((output / "katie_connected/RUN_RECEIPT.json").read_text())
    require(katie_receipt["input_sha256"] == sha(source_file), "Katie source binding differs")
    raw = load4(output / "katie_connected/raw_katie.rga")
    signature_file = output / "katie_connected/triangulation.isosig"
    require(sha(signature_file) == katie_receipt["isosig_sha256"] and
            raw.isoSig() == signature_file.read_text().strip(), "Katie output signature binding differs")
    raw_check = audit_ideal(raw)
    print("Checking ideal truncation and the finite lower-body reference", flush=True)
    truncated = regina.Triangulation4(raw)
    require(truncated.truncateIdeal(), "ideal truncation made no change")
    require(truncated == load4(output / "truncated_raw.rga"), "labelled truncation replay differs")
    require(not truncated.isIdeal() and truncated.countBoundaryComponents() == 1 and
            truncated.boundaryComponent(0).isReal(), "truncation did not produce one real boundary")
    lower_file = output / "lower_reference.rga"
    lower = load4(lower_file)
    preparation = json.loads((output / "LOWER_REFERENCE.json").read_text())
    require(preparation["source_truncated_sha256"] == sha(output / "truncated_raw.rga") and
            preparation["result_reference_sha256"] == sha(lower_file),
            "native lower-body simplification receipt binding differs")
    lower_check = cap_audit(lower)
    require(lower.homology().isTrivial() and lower.homology(2).rank() == 2 and
            lower.homology(2).isFree(2), "wrong lower-body homology")
    path_file = output / "selected_boundary_path.json"
    path = json.loads(path_file.read_text())
    require(boundary_data(lower)[0].isoSig() == path["start"], "boundary path starts elsewhere")
    require(path["end"] == GOAL, "boundary path has wrong endpoint")
    path_check = replay_boundary(path)
    print("Replaying actual 4D layerings and reverse boundary shellings", flush=True)
    layered = load4(output / "lower_layered/cap_with_matched_boundary.rga")
    layer_receipt = json.loads((output / "lower_layered/LAYERING_VERIFICATION.json").read_text())
    require(layer_receipt["source_cap_sha256"] == sha(lower_file) and
            layer_receipt["source_path_sha256"] == sha(path_file), "layering source binding differs")
    layer_check = replay_layerings(lower, path, output / "lower_layered/layering_events.jsonl.gz", layered)
    cap_file = PROJECT / CAP_RELATIVE
    require(sha(cap_file) == CAP_SHA256, "standard upper-cap checksum differs")
    cap = load4(cap_file)
    require(boundary_data(cap)[0].isoSig() == GOAL, "standard cap boundary differs")
    closed = load4(output / "closed/closed_triangulation.rga")
    print("Checking the boundary isomorphism, closed gluing and vertex links", flush=True)
    gluing_receipt = json.loads((output / "closed/CLOSED_GLUING_VERIFICATION.json").read_text())
    require(gluing_receipt["source_first_sha256"] == sha(output / "lower_layered/cap_with_matched_boundary.rga")
            and gluing_receipt["source_second_sha256"] == CAP_SHA256, "closed gluing source binding differs")
    closed_check = replay_gluing(layered, cap, output / "closed/boundary_gluing.json", closed)
    report = {"schema": "motoo.genus_two_completion_control.v1",
              "status": "PASS_EXPLICIT_GENUS_TWO_COMPLETION_CONTROL",
              "source_unlink": source_check, "raw_ideal_body": raw_check,
              "truncated_pentachora": truncated.size(), "labelled_truncation_reproduced": True,
              "finite_lower_body": lower_check, "boundary_path": path_check,
              "actual_4d_layering_replay": layer_check, "standard_cap_sha256": CAP_SHA256,
              "standard_cap_pentachora": cap.size(), "closed_control": closed_check,
              "native_trust_boundaries": ["Katie 7.4.1 diagram-to-triangulation algorithm",
                   "Regina ideal truncation", "Regina initial finite-body simplify",
                   "Regina guarded moves and 3-sphere recognition"],
              "initial_native_simplification_has_individual_move_log": False,
              "native_preparation_receipt": preparation,
              "standardness_not_inferred_from_invariants": True,
              "unknown_motoo_diagram_converted": False,
              "regina_python_distribution": importlib.metadata.version("regina"),
              "checker_sha256": sha(__file__)}
    negative = output / "katie/raw_katie.rga"
    if negative.exists():
        split = load4(negative)
        require(split.countComponents() == 2, "historical disconnected-projection control changed")
        report["disconnected_projection_negative_control"] = {
            "components": 2, "would_fail_single_body_admission": True,
            "warning": "A disconnected planar projection must not be used as if Katie automatically forms a boundary connected sum."}
    save_json(output / "COMPLETION_CONTROL.json", report)
    print(json.dumps({"status": report["status"], "closed_pentachora": closed.size(),
                      "layering_moves": len(path["moves"]), "output": str(output)}), flush=True)
    return report


def build(output, seconds, seed, katie=None):
    output = Path(output)
    require(not output.exists(), "choose a new output directory")
    output.mkdir(parents=True)
    link = regina.Link.fromData([1, -1], [[1, 2], [-1, -2]])
    end = regina.Link(link)
    require(end.r2(end.crossing(0)), "control unlink move failed")
    data = {"schema": "motoo.katie.control.v1", "name": "two_zero_unlinks_connected_projection",
            "pd_code": [[x - 1 for x in row] for row in link.pdData()], "decorations": [0, 0],
            "unlink_proof": {"initial_tight_encoding": link.tightEncoding(), "r2_crossing": 0,
                             "final_tight_encoding": end.tightEncoding()}}
    save_json(output / "connected_input.json", data)
    command = [sys.executable, str(PROJECT / "code/run_katie.py"), str(output / "connected_input.json"),
               "--out", str(output / "katie_connected")]
    if katie:
        command += ["--katie", str(katie)]
    subprocess.run(command, check=True)
    lower = load4(output / "katie_connected/raw_katie.rga")
    audit_ideal(lower)
    require(lower.truncateIdeal(), "control ideal truncation failed")
    save4(lower, output / "truncated_raw.rga")
    lower.simplify()
    save4(lower, output / "lower_reference.rga")
    save_json(output / "LOWER_REFERENCE.json", {
        "operation": "Regina Triangulation4.simplify()",
        "source_truncated_sha256": sha(output / "truncated_raw.rga"),
        "result_reference_sha256": sha(output / "lower_reference.rga"),
        "native_regina_simplification_used": True,
        "individual_move_log_available": False,
        "regina_python_distribution": importlib.metadata.version("regina"),
        "pentachora": lower.size(), "isosig": lower.isoSig(),
    })
    boundary = boundary_data(lower)[0]
    save3(boundary, output / "lower_boundary.rga")
    path = search(boundary.isoSig(), output / "boundary_search", seconds=seconds,
                  seed=seed, goals={GOAL})
    require(path["status"] == "FOUND_AND_REPLAYED", "boundary search unfinished; checkpoint retained")
    save_json(output / "selected_boundary_path.json", path)
    layer_path(output / "lower_reference.rga", output / "selected_boundary_path.json", output / "lower_layered")
    glue_pieces(output / "lower_layered/cap_with_matched_boundary.rga", PROJECT / CAP_RELATIVE, output / "closed")
    return audit(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--build", type=Path, help="new output directory for a fresh construction")
    choice.add_argument("--audit", type=Path, help="replay a saved completed control directory")
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--katie", type=Path)
    args = parser.parse_args()
    require(importlib.metadata.version("regina") == "7.4.1", "requires Regina Python 7.4.1")
    if args.build:
        build(args.build, args.seconds, args.seed, args.katie)
    else:
        audit(args.audit)

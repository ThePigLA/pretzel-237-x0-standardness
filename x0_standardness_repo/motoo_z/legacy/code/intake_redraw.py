#!/usr/bin/env python3
"""Check a marked six/eight-component redraw before Katie conversion.

This tool checks data consistency, classical planarity, the source-component
correspondence, and Katie position. It can separately rerun native group and
boundary checks. Neither matching invariants nor a declared cancellation is
accepted as proof that the new diagram represents the original manifold.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

import regina
from spherogram import Link

from orient_katie_diagram import components_of_graph, orient
from run_katie import inspect_input


SCHEMA = "motoo.incoming-marked-redraw.v1"
ROOT = Path(__file__).resolve().parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def zero_pd(pd):
    require(isinstance(pd, list) and pd, "PD must contain crossings")
    require(all(isinstance(c, list) and len(c) == 4 for c in pd), "PD tuple shape")
    require(all(type(x) is int and x >= 0 for c in pd for x in c), "PD arc label type")
    counts = Counter(x for c in pd for x in c)
    require(set(counts.values()) == {2}, "Each PD arc must occur exactly twice")
    origin = min(counts)
    require(origin in (0, 1) and set(counts) == set(range(origin, origin + 2 * len(pd))),
            "PD labels must be contiguous from zero or one")
    return [[x-origin for x in c] for c in pd], origin


def template(figure):
    original = ROOT / f"inputs/figure{figure}_link.json"
    data = json.loads(original.read_text())
    return {
        "schema": SCHEMA,
        "source": {"figure": figure, "original_input_sha256": sha(original)},
        "relationship": {"kind": "isotopy", "description": "Unchanged original PD; replace with an actual redraw and its component correspondence."},
        "pd_code": data["pd_code"], "decorations": data["decorations"],
        "components": [{"source_id": c["id"], "arc_labels": c["arc_labels"],
                        "handle_index": c["handle_index"], "framing": c["framing"]}
                       for c in data["components"]],
    }


def inspect(document):
    require(document.get("schema") == SCHEMA, "Wrong redraw input schema")
    source = document.get("source")
    require(isinstance(source, dict), "Source figure and original-input binding are required")
    figure = source.get("figure")
    require(type(figure) is int and figure in (1, 2), "Source figure must be 1 or 2")
    original_path = ROOT / f"inputs/figure{figure}_link.json"
    require(source.get("original_input_sha256") == sha(original_path), "Original input checksum/figure mismatch")
    original = json.loads(original_path.read_text())
    relationship = document.get("relationship")
    require(isinstance(relationship, dict), "Explicit relationship metadata is required")
    kind = relationship.get("kind")
    require(kind in ("isotopy", "lower_pair_cancelled"), "Supported relationship kinds: isotopy, lower_pair_cancelled")
    require(isinstance(relationship.get("description"), str) and relationship["description"].strip(),
            "Describe the claimed redraw or handle-move relationship")
    pd, origin = zero_pd(document["pd_code"])
    records = document.get("components")
    require(isinstance(records, list), "Explicit source-component records are required")
    expected_ids = {f"C{i}" for i in range(8)}
    if kind == "lower_pair_cancelled":
        expected_ids -= {"C4", "C7"}
        require(relationship.get("cancelled_source_components") == ["C4", "C7"],
                "The supported lower cancellation removes C4 and C7")
    require(len(records) == len(expected_ids), "Wrong component count for the declared relationship")
    ids = [c.get("source_id") for c in records]
    require(set(ids) == expected_ids and len(set(ids)) == len(ids), "Source-component correspondence is not the required bijection")
    old_by_id = {c["id"]: c for c in original["components"]}
    source_sets = {}
    for record in records:
        identity = record["source_id"]
        old = old_by_id[identity]
        handle = record.get("handle_index")
        framing = record.get("framing")
        require(type(handle) is int and handle == old["handle_index"], "A dotted/framed source role changed")
        require(framing is None if handle == 1 else type(framing) is int,
                "Dots require null framing; 2-handles require integer framing")
        if kind == "isotopy":
            require(framing == old["framing"], "An isotopy must preserve the integral framing")
        arcs = record.get("arc_labels")
        require(isinstance(arcs, list) and arcs and all(type(x) is int for x in arcs),
                "Every marked component needs its nonempty PD arc-label set")
        require(len(set(arcs)) == len(arcs), "Repeated arc in component metadata")
        arcs = frozenset(x-origin for x in arcs)
        require(arcs not in source_sets, "Duplicate component arc set")
        source_sets[arcs] = record
    require(sum(len(x) for x in source_sets) == len(set().union(*source_sets)), "Marked component arc sets overlap")
    require(set().union(*source_sets) == set(range(2 * len(pd))), "Marked component arcs do not cover the PD")

    # Two independent native PD imports; neither is asked to simplify the link.
    link = Link(pd, check_planarity=True)
    rlink = regina.Link.fromPD([[x+1 for x in c] for c in pd])
    require(len(link.link_components) == rlink.countComponents() == len(records), "Native PD component counts differ")
    require(len(link.crossings) == rlink.size() == len(pd), "Native PD crossing counts differ")
    order = []
    native_components = []
    owner = {}
    for i in range(rlink.countComponents()):
        start = rlink.component(i)
        at = start
        arcs = set()
        visits = []
        for _ in range(2*len(pd)+1):
            crossing, level = at.crossing().index(), at.strand()
            arcs.update((pd[crossing][level], pd[crossing][level+2]))
            visits.append({"crossing": crossing, "role": "over" if level else "under"})
            require((crossing, level) not in owner, "Repeated component traversal strand")
            owner[crossing, level] = i
            at = at.next()
            if at == start:
                break
        else:
            raise ValueError("Component traversal did not close")
        require(frozenset(arcs) in source_sets, "Native component has no matching source arc set")
        c = source_sets[frozenset(arcs)]
        order.append(c["source_id"])
        native_components.append({"index": i, "id": c["source_id"], "name": c["source_id"],
                                  "handle_index": c["handle_index"], "dotted": c["handle_index"] == 1,
                                  "framing": c["framing"], "arc_labels": sorted(arcs),
                                  "crossing_visits": visits})
    for i, component in enumerate(link.link_components):
        arcs = set()
        for ep in component:
            j, level = int(ep.crossing.label), ep.strand_index % 2
            arcs.update((pd[j][level], pd[j][level+2]))
        require(arcs == set(native_components[i]["arc_labels"]), "Spherogram/Regina component-order disagreement")
    crossings = []
    for i, c in enumerate(link.crossings):
        require(int(c.label) == i, "Spherogram relabelled crossings")
        require(c.sign == rlink.crossing(i).sign(), "Native crossing signs disagree")
        require(c.strand_components[0] == owner[i, 0] and c.strand_components[1] == owner[i, 1],
                "Native marked crossing ownership disagrees")
        crossings.append({"id": i, "sign": int(c.sign), "under_component": owner[i, 0],
                          "over_component": owner[i, 1]})
    decorations = ["x" if c["dotted"] else c["framing"] for c in native_components]
    require(document.get("decorations") == decorations, "Decoration vector disagrees with independently matched component records")
    n1 = sum(c["dotted"] for c in native_components)
    require((n1, len(records)-n1) == ((3, 5) if kind == "isotopy" else (2, 4)), "Wrong handle counts")

    # Projection connectivity is distinct from link component count and from
    # whether the represented smooth 4-manifold ought to be connected.
    occurrences = defaultdict(list)
    for i, c in enumerate(pd):
        for arc in c:
            occurrences[arc].append(i)
    adjacency = defaultdict(set)
    for a, b in occurrences.values():
        adjacency[a].add(b)
        adjacency[b].add(a)
    projection_components = components_of_graph(range(len(pd)), adjacency)
    native = {"schema": "motoo.native-redraw-data.v1", "figure": figure,
              "pd_code": pd, "decorations": decorations,
              "pd_component_to_source": list(range(len(records))),
              "original_source_ids_in_pd_order": order,
              "components": native_components, "crossings": crossings,
              "component_count": len(records), "crossing_count": len(pd)}
    signature = {"pd_code": pd, "decorations": decorations, "source_ids": order}
    original_signature = {"pd_code": zero_pd(original["pd_code"])[0],
                          "decorations": original["decorations"],
                          "source_ids": [original["components"][i]["id"] for i in original["pd_component_to_source"]]}
    exact = kind == "isotopy" and signature == original_signature
    shape = inspect_input(native)
    oriented = None
    orientation_error = None
    if len(projection_components) == 1:
        try:
            # Export only metadata that survives the PD relabelling.
            minimal = {"pd_code": pd, "decorations": decorations,
                       "original_source_ids_in_pd_order": order,
                       "source_equivalence_verified": exact}
            oriented = orient(minimal)
            after = inspect_input(oriented)
            require(not after["necessary_shape_failures"], "Oriented output still fails position")
        except ValueError as error:
            orientation_error = str(error)
            oriented = None
    else:
        orientation_error = "Disconnected projection needs normalization for this Katie route"
    report = {
        "schema": "motoo.redraw-intake-report.v1", "data_consistency": "PASS",
        "source_figure": figure, "original_input_sha256": sha(original_path),
        "relationship_declared": relationship,
        "source_ids_in_pd_order": order, "handle_counts_0_1_2": [1, n1, len(records)-n1],
        "crossings": len(pd), "crossing_signs_and_component_ownership_agree": True,
        "spherogram_classical_planarity_passed": True,
        "projection_connected": len(projection_components) == 1,
        "projection_component_count": len(projection_components),
        "disconnected_projection_is_not_a_topological_invalidity_claim": True,
        "source_equivalence_status": "EXACT_ORIGINAL_MARKED_PD" if exact else "PENDING_GEOMETRIC_RECORD",
        "source_equivalence_verified": exact,
        "source_equivalence_note": ("A declaration, component map, group or boundary match is not an equivalence certificate."
                                    " Isotopic redraws need a checked marked-isotopy chain; cancellations/slides need a recorded Kirby argument."),
        "initial_katie_position": shape,
        "katie_position_after_auxiliary_orientation": "PASS" if oriented else "NEEDS_NORMALIZATION",
        "orientation_error": orientation_error,
        "native_group_status": "NOT_RUN", "native_boundary_status": "NOT_RUN",
        "closed_4d_triangulation_constructed": False, "standardness_proved": False,
    }
    return native, oriented, report


def check(input_path, output, native_checks=False, timeout=180):
    input_path, output = Path(input_path), Path(output)
    require(not output.exists(), "Choose a new output directory")
    output.mkdir(parents=True)
    raw = input_path.read_bytes()
    (output / "incoming.json").write_bytes(raw)
    try:
        native, oriented, report = inspect(json.loads(raw))
    except Exception as error:
        save(output / "INTAKE.json", {"status": "REJECTED_INPUT_DATA", "input_sha256": hashlib.sha256(raw).hexdigest(),
                                      "error": str(error), "source_equivalence_verified": False,
                                      "closed_4d_triangulation_constructed": False})
        raise
    report["input_sha256"] = hashlib.sha256(raw).hexdigest()
    report["checker_sha256"] = sha(__file__)
    report["versions"] = {p: importlib.metadata.version(p) for p in ("regina", "snappy", "spherogram")}
    save(output / "native_input.json", native)
    report["native_input_sha256"] = sha(output / "native_input.json")
    if oriented:
        save(output / "oriented_katie_input.json", oriented)
        report["oriented_input_sha256"] = sha(output / "oriented_katie_input.json")
    if native_checks:
        for stage, command, result_path in [
            ("group", [sys.executable, str(ROOT / "code/check_handlebody_group.py"), "compute",
                       str((output / "native_input.json").resolve()), str((output / "group.json").resolve())], output / "group.json"),
            ("boundary", [sys.executable, str(ROOT / "code/check_boundary.py"), str((output / "native_input.json").resolve()),
                          "--out", str((output / "boundary").resolve()), "--expected", "2"], output / "boundary/BOUNDARY_VERIFICATION.json"),
        ]:
            logfile = output / f"{stage}.log"
            try:
                with logfile.open("w") as stream:
                    subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                   timeout=timeout, check=True)
                native_report = json.loads(result_path.read_text())
                report[f"native_{stage}_status"] = native_report["status"]
                report[f"native_{stage}_receipt_sha256"] = sha(result_path)
            except subprocess.TimeoutExpired:
                report[f"native_{stage}_status"] = "TIMEOUT_NOT_A_NEGATIVE_TOPOLOGICAL_RESULT"
            except subprocess.CalledProcessError:
                report[f"native_{stage}_status"] = "CHECK_FAILED_SEE_LOG"
    report["ready_for_project_conversion"] = bool(report["source_equivalence_verified"] and oriented
        and report["projection_connected"]
        and report["native_group_status"] == "TRIVIAL_HANDLEBODY_PI1"
        and report["native_boundary_status"] == "PASS_EXPECTED_BOUNDARY")
    report["status"] = "READY_FOR_PROJECT_CONVERSION" if report["ready_for_project_conversion"] else "CHECKED_DATA_WITH_OPEN_GATES"
    save(output / "INTAKE.json", report)
    print(json.dumps({"status": report["status"], "data_consistency": "PASS",
                      "source_equivalence": report["source_equivalence_status"],
                      "katie_position": report["katie_position_after_auxiliary_orientation"],
                      "group": report["native_group_status"], "boundary": report["native_boundary_status"],
                      "output": str(output)}), flush=True)
    return report


def selftest(output):
    output = Path(output)
    require(not output.exists(), "Choose a new output directory")
    output.mkdir(parents=True)
    base = template(1)
    save(output / "figure1_unchanged_template.json", base)
    _, _, positive = inspect(base)
    require(positive["source_equivalence_verified"], "Exact original was not recognized")
    require(positive["katie_position_after_auxiliary_orientation"] == "NEEDS_NORMALIZATION", "Original diagram unexpectedly admitted")
    cases = []
    changes = [
        ("missing_relationship", lambda d: d.pop("relationship")),
        ("missing_figure", lambda d: d["source"].pop("figure")),
        ("wrong_source_hash", lambda d: d["source"].update(original_input_sha256="0"*64)),
        ("wrong_component_count", lambda d: d["components"].pop()),
        ("swapped_dot_framing_vector", lambda d: d["decorations"].__setitem__(slice(0, 2), ["x", 0])),
        ("changed_source_handle_role", lambda d: d["components"][0].update(handle_index=1, framing=None)),
        ("changed_isotopy_framing", lambda d: d["components"][3].update(framing=1)),
        ("missing_component_arc", lambda d: d["components"][0]["arc_labels"].pop()),
    ]
    for name, mutate in changes:
        d = copy.deepcopy(base)
        mutate(d)
        try:
            inspect(d)
        except ValueError as error:
            cases.append({"case": name, "status": "REJECTED_AS_EXPECTED", "reason": str(error)})
        else:
            raise ValueError("Negative input was accepted: " + name)
    unproved = copy.deepcopy(base)
    unproved["pd_code"] = unproved["pd_code"][1:] + unproved["pd_code"][:1]
    unproved["relationship"]["description"] = "Relabelled crossing list supplied without an equivalence record"
    _, _, pending = inspect(unproved)
    require(not pending["source_equivalence_verified"], "An unverified equivalence declaration was promoted")
    save(output / "unproved_isotopy_claim.json", unproved)
    cases.append({"case": "unproved_isotopy_claim", "status": "PENDING_AS_EXPECTED"})

    # A six-component star link tests the reduced-count schema and the
    # successful position/orientation route. It is explicitly unrelated to
    # Motoo's diagram; its input declaration is not an equivalence proof.
    cycles = [sum(([-2*i-1, 2*i+2] for i in range(5)), [])]
    cycles += [[2*i+1, -2*i-2] for i in range(5)]
    fixture_link = regina.Link.fromData([1]*10, cycles)
    fixture_pd = [[x-1 for x in c] for c in fixture_link.pdData()]
    fixture_link = regina.Link.fromPD([[x+1 for x in c] for c in fixture_pd])
    fixture = copy.deepcopy(base)
    fixture["pd_code"] = fixture_pd
    fixture["decorations"] = [0, "x", "x", -1, 0, 0]
    fixture["components"] = []
    identities = ["C0", "C1", "C2", "C3", "C5", "C6"]
    for i in range(6):
        at = start = fixture_link.component(i)
        arcs = set()
        while True:
            j, level = at.crossing().index(), at.strand()
            arcs.update((fixture_pd[j][level], fixture_pd[j][level+2]))
            at = at.next()
            if at == start:
                break
        decoration = fixture["decorations"][i]
        fixture["components"].append({"source_id": identities[i], "arc_labels": sorted(arcs),
                                      "handle_index": 1 if decoration == "x" else 2,
                                      "framing": None if decoration == "x" else decoration})
    fixture["relationship"] = {"kind": "lower_pair_cancelled", "cancelled_source_components": ["C4", "C7"],
        "description": "UNRELATED star-link control; this is deliberately not a proved redraw of Motoo's diagram."}
    native_fixture, oriented_fixture, fixture_report = inspect(fixture)
    require(fixture_report["handle_counts_0_1_2"] == [1, 2, 4], "Six-component count fixture failed")
    require(oriented_fixture is not None, "Simple two-dot position fixture was refused")
    require(not fixture_report["source_equivalence_verified"], "Position success was mistaken for source equivalence")
    require("components" not in oriented_fixture and "crossings" not in oriented_fixture,
            "Stale arc/crossing metadata leaked into reoriented Katie input")
    save(output / "six_component_position_fixture.json", fixture)
    save(output / "six_component_position_report.json", fixture_report)
    save(output / "six_component_oriented_input.json", oriented_fixture)
    cases.append({"case": "six_component_position_pass_source_pending", "status": "PASS"})
    result = {"status": "PASS", "exact_original_case": positive,
              "negative_cases": cases, "cases": 1+len(cases),
              "scope": "Data and provenance gates, not a new Motoo redraw or 4D conversion"}
    save(output / "SELFTEST.json", result)
    print(json.dumps({"status": "PASS", "cases": result["cases"], "output": str(output)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    t = sub.add_parser("template")
    t.add_argument("--figure", type=int, choices=(1, 2), required=True)
    t.add_argument("--output", type=Path, required=True)
    c = sub.add_parser("check")
    c.add_argument("input", type=Path)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--native-checks", action="store_true")
    c.add_argument("--timeout", type=float, default=180, help="seconds per native group/boundary subprocess")
    s = sub.add_parser("selftest")
    s.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(importlib.metadata.version("regina") == "7.4.1", "Requires Regina Python distribution 7.4.1")
    if args.command == "template":
        require(not args.output.exists(), "Refusing to overwrite an existing template")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        save(args.output, template(args.figure))
    elif args.command == "check":
        check(args.input, args.out, args.native_checks, args.timeout)
    else:
        selftest(args.out)


if __name__ == "__main__":
    main()

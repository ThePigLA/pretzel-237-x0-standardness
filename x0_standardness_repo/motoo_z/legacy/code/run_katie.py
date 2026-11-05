#!/usr/bin/env python3
"""Run pinned Katie on an explicitly decorated PD diagram, retaining raw output.

Input JSON: {"pd_code": [[...], ...], "decorations": ["x", 0, ...]}.
The decoration order must be Regina's Link.fromPD component order. This runner
checks sufficient planar position conditions for dotted components, including
counterclockwise orientations. Source-diagram identity still requires a separate
isotopy/transcription admission record.
It does not fill ideal vertices, simplify the output, or claim standardness.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import time

import regina


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect_input(data: dict) -> dict:
    pd = data["pd_code"]
    decorations = data["decorations"]
    if not isinstance(pd, list) or not pd:
        raise ValueError("pd_code must be a nonempty list")
    if any(not isinstance(c, list) or len(c) != 4 or
           any(type(x) is not int or x < 0 for x in c) for c in pd):
        raise ValueError("Every crossing must be four nonnegative integer arc labels")
    occurrences = Counter(x for c in pd for x in c)
    if any(n != 2 for n in occurrences.values()):
        raise ValueError("Every arc label must occur exactly twice")
    if set(occurrences) not in (set(range(2 * len(pd))), set(range(1, 2 * len(pd) + 1))):
        raise ValueError("Arc labels must be contiguous, starting from zero or one")
    if any(not (type(d) is int or d in ("x", ".")) for d in decorations):
        raise ValueError("Decorations must be integer framings or x / .")
    # Katie's graph construction does not join disconnected diagram
    # projections into the one-0-handle Kirby interpretation used here.
    # In particular, two isolated zero-framed unknot projections produce
    # two disconnected 4D pieces. This is an unsupported input position,
    # not a claim that the underlying link or Kirby diagram is invalid.
    from orient_katie_diagram import components_of_graph
    arc_crossings = defaultdict(list)
    for crossing, row in enumerate(pd):
        for arc in row:
            arc_crossings[arc].append(crossing)
    projection_adjacency = {i: set() for i in range(len(pd))}
    for first, second in arc_crossings.values():
        projection_adjacency[first].add(second)
        projection_adjacency[second].add(first)
    projection_parts = components_of_graph(range(len(pd)), projection_adjacency)
    # Katie accepts SnapPy's zero-based codes and shifts them to Regina's
    # required one-based range. Mirror precisely that harmless label change.
    label_origin = min(occurrences)
    regina_pd = [[x + 1 for x in crossing] for crossing in pd] if label_origin == 0 else pd
    link = regina.Link.fromPD(regina_pd)
    if link.countComponents() != len(decorations):
        raise ValueError("Decoration count differs from PD component count")
    components = []
    owners = {}
    for i in range(link.countComponents()):
        start = link.component(i)
        at = start
        visits = []
        for _ in range(2 * link.size() + 1):
            crossing = at.crossing().index()
            strand = at.strand()
            visits.append([crossing, strand])
            if (crossing, strand) in owners:
                raise ValueError("Repeated strand in component traversal")
            owners[crossing, strand] = i
            at = at.next()
            if at == start:
                break
        else:
            raise ValueError("Component traversal did not close")
        over_under = [v[1] for v in visits]
        changes = sum(a != b for a, b in zip(over_under, over_under[1:] + over_under[:1]))
        components.append({"index": i, "decoration": decorations[i],
                           "visits_crossing_strand": visits,
                           "cyclic_over_under_transitions": changes})
    for c in components:
        i = c["index"]
        c["self_crossings"] = [j for j in range(link.size())
                                if owners[j, 0] == owners[j, 1] == i]
    dotted = {i for i, d in enumerate(decorations) if d in ("x", ".")}
    dotted_crossings = [j for j in range(link.size())
                       if owners[j, 0] in dotted and owners[j, 1] in dotted]
    failures = []
    if len(projection_parts) != 1:
        failures.append("Disconnected diagram projections are unsupported by this Katie conversion route; "
                        "supply a connected projection of the same marked link")
    for i in sorted(dotted):
        c = components[i]
        if c["self_crossings"]:
            failures.append(f"Dotted component {i} is not a simple projected loop")
        if c["cyclic_over_under_transitions"] != 2:
            failures.append(f"Dotted component {i} does not have one overpass block and one underpass block")
    # Disjoint projected dotted circles are a sufficient, deliberately strong
    # diagrammatic condition; more general diagrams need separate isotopies.
    if dotted_crossings:
        failures.append("Dotted-component projections have crossings with dotted components")
    planar = None
    if not failures and dotted:
        from orient_katie_diagram import admission as check_planar_position
        try:
            geometric = data.get("geometry_admission") or {}
            exterior = (geometric.get("after") or {}).get("exterior_dart")
            planar = check_planar_position(regina_pd, dotted, exterior)
            if not planar["all_dotted_counterclockwise"]:
                failures.append("Dotted components require orientation changes: " +
                                str(planar["reverse_components"]))
        except ValueError as error:
            failures.append("Planar dotted-disc check failed: " + str(error))
    return {"crossings": link.size(), "input_arc_label_origin": label_origin,
            "component_count": link.countComponents(),
            "projection_connected": len(projection_parts) == 1,
            "projection_components": [sorted(part) for part in projection_parts],
            "components": components, "dotted_indices": sorted(dotted),
            "dotted_projection_crossings": dotted_crossings,
            "necessary_shape_failures": failures,
            "counterclockwise_orientation_checked": bool(not dotted or
                planar and planar["all_dotted_counterclockwise"]),
            "planar_disc_admission": planar,
            "geometry_admission": data.get("geometry_admission"),
            "source_diagram_identity_checked_by_this_runner": False}


def triangulation_report(tri: regina.Triangulation4) -> dict:
    result = {"pentachora": tri.size(), "vertices": tri.countVertices(),
              "valid": tri.isValid(), "closed": tri.isClosed(),
              "ideal": tri.isIdeal(), "orientable": tri.isOrientable(),
              "connected": tri.isConnected(),
              "euler_characteristic_triangulation": tri.eulerCharTri(),
              "euler_characteristic_manifold": tri.eulerCharManifold()}
    if tri.isValid():
        result["homology"] = {str(k): str(tri.homology(k)) for k in (1, 2)}
    links = []
    for vertex in tri.vertices():
        link = vertex.buildLink()
        links.append({"vertex": vertex.index(), "tetrahedra": link.size(),
                      "valid": link.isValid(), "closed": link.isClosed(),
                      "sphere": link.isSphere(), "homology": str(link.homology())})
    result["vertex_links"] = links
    result["all_vertex_links_spheres"] = all(x["sphere"] for x in links)
    result["closed_manifold_admission"] = bool(tri.isValid() and tri.isClosed()
                                              and not tri.isIdeal()
                                              and result["all_vertex_links_spheres"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--katie", type=Path)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--exploratory", action="store_true",
                        help="Permit necessary dotted-position tests to fail; output remains unadmitted")
    args = parser.parse_args()
    if importlib.metadata.version("regina") != "7.4.1":
        raise RuntimeError("This runner requires the regina 7.4.1 Python distribution")
    project = Path(__file__).resolve().parents[1]
    local = project / "scratch/katie/prefix/usr/bin/katie"
    executable = (args.katie or local).resolve()
    if not executable.is_file():
        raise FileNotFoundError("Katie executable absent; set --katie or run setup_katie_linux.py")
    data = json.loads(args.input.read_text())
    admission = inspect_input(data)
    if admission["necessary_shape_failures"] and not args.exploratory:
        raise ValueError("Dotted-position checks failed: " + "; ".join(admission["necessary_shape_failures"]))
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "input.json").write_text(json.dumps(data, indent=2) + "\n")
    (args.out / "input_checks.json").write_text(json.dumps(admission, indent=2) + "\n")
    env = dict(os.environ)
    if executable == local.resolve():
        env["LD_LIBRARY_PATH"] = str(project / "scratch/katie/prefix/usr/lib/x86_64-linux-gnu")
    version = subprocess.run([str(executable), "--version"], env=env, capture_output=True,
                             text=True, timeout=30, check=True).stdout.strip()
    if "7.4.1" not in version:
        raise RuntimeError("Katie binary is not the pinned 7.4.1 release")
    command = [str(executable), json.dumps(data["pd_code"]),
               " ".join(str(x) for x in data["decorations"]), "--dim4", "--verbose"]
    start = time.monotonic()
    with (args.out / "katie.stdout").open("w") as stdout, (args.out / "katie.stderr").open("w") as stderr:
        completed = subprocess.run(command, env=env, stdout=stdout, stderr=stderr,
                                   timeout=args.timeout)
    receipt = {"input_sha256": sha256(args.input), "katie_sha256": sha256(executable),
               "katie_version": version, "regina_python_distribution": "7.4.1",
               "argv": command, "returncode": completed.returncode,
               "elapsed_seconds": time.monotonic() - start,
               "input_checks": admission, "simplification_performed": False,
               "three_handle_filling_performed": False}
    (args.out / "RUN_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if completed.returncode:
        raise RuntimeError(f"Katie failed with exit status {completed.returncode}; see katie.stderr")
    # In 7.4.1, --verbose emits one link-walk diagnostic to stdout as well
    # as the intended stderr diagnostics. The signature is the final line.
    stdout_lines = [line.strip() for line in (args.out / "katie.stdout").read_text().splitlines()
                    if line.strip()]
    if not stdout_lines:
        raise RuntimeError("Katie returned no output")
    signature = stdout_lines[-1]
    tri = regina.Triangulation4.fromIsoSig(signature)
    (args.out / "triangulation.isosig").write_text(signature + "\n")
    packet = regina.PacketOfTriangulation4(tri)
    packet.setLabel(data.get("name", args.input.stem) + " -- raw Katie output")
    if not packet.save(str(args.out / "raw_katie.rga")):
        raise RuntimeError("Could not save raw native triangulation")
    receipt["triangulation"] = triangulation_report(tri)
    receipt["status"] = ("RAW_FINITE_CLOSED_TRIANGULATION_GEOMETRY_REQUIRES_SEPARATE_ADMISSION"
                         if receipt["triangulation"]["closed_manifold_admission"] else
                         "RAW_IDEAL_OR_BOUNDED_TRIANGULATION_NOT_CLOSED_COMPLETION")
    receipt["isosig_sha256"] = sha256(args.out / "triangulation.isosig")
    (args.out / "RUN_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "output": str(args.out),
                      "pentachora": tri.size(), "closed": tri.isClosed(),
                      "ideal": tri.isIdeal()}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Reconstruct the published 24,170-pentachoron Z^circ triangulation.

Starts from the admitted marked PD code, runs Katie (or explicitly starts at
its saved raw output), and constructs the finite pieces and their gluing.
Recorded boundary paths are replayed; their discovery is not rerun.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from common import (ROOT, LEGACY, REFERENCE, check_integrity, fresh_output,
                    require, run_stage, sha256, write_json)

GUARDED_CAP = ROOT / "legacy/results/figure2_working_20261102/upper_cap_guarded_v1"
SMALL_CAP = ROOT / "legacy/results/boundary/upper_cap_selftest_v1/reference_cap.rga"
MATCHED_CAP = ROOT / "legacy/results/boundary/standard_cap_six_boundary_v1"
PD_SHA = "5aeaea4ffa70b3d4084fbae11cd144c256ae3af361cbebb10da07166cad3427f"
ORIENTED_SHA = "c6b4d9c7522a69b4df67d35029638480ff74b62e78ee5c3377a25c0b311f8acd"
PATH_SHA = "f1014e73417e65bfdf947183e9d339b4f7d6b7843d1c3f7a22f8216124b2f870"
HISTORICAL_KATIE_SHA = "f9aa41346517db1179463c92f53f22aff7e1b6f09d7677164bbe7b4d3fba4b60"

# The expected file hashes bind the archived witnesses. For generated files
# enforce exact labelled-state and uncompressed XML identity: zlib builds can
# produce different gzip bytes for identical serialized triangulations.
# Except in --from-raw mode, witnesses are never substituted for generated
# construction states.
CHECKPOINTS = {
    "raw": (1528, "f1590ec6de3e7a9d9fdeb2175f98a9cd70688027776d274de0f2da4718a119bb",
            "d4750392cf2661bc4af787afc2a253235fb62357b17e1ca0266c8de0e38f4478"),
    "finite": (6112, "3da4605e87047418b42f5464182aa928b55e5d04be6142b67b35395a4e501665",
               "9af38481e4f93cfe5a1505bb5a3e927e88b51f1aa6692e678d10db10dfa0302d"),
    "lower": (1392, "754122baa2a4a381f0c4859e2cd883cc8538d5513a83eed6134fa7b302f8287d",
              "a10761b045c3a22933f66200f192bacb8b1855a6748ba3eb134ebcc724276c54"),
    "upper24": (24, "5d435806c52d64b056bd8b02e33de0887613948a82acebe6025e04d92bf401f0",
                "02467c65e6bc171745b3304764fd63d40b40f7cc9846b56448fa214b66574d1f"),
    "upper9": (9, "6ff5b588fb2cf655e1a5a40d630312c5fe82da9a2327ea22f2ca5176ab1c0cb1",
               "dd1f0a488ffc68591e19f685d8b8fb729211125f2fa778bf5cd9d00f09f8fbe1"),
    "upper474": (474, "2b46b7e8d760ce52641f78cb78c998748a8e476410b9c6525ae1b28168ca31b7",
                 "9ddcd8590142bf08cb94ad7bbda02d38506cef96f47d23f1e6fde7c70ca58e78"),
    "layered": (23696, "e0d1b39192ddfc7fc8852d210658d7e0647a0ba0d0a36b8fe07c704483bc9a37",
                "e26233b6916921c97e95318f6ff73628617b303bc25a5af46a969a3436c9a95b"),
    "closed": (24170, "e6281897c358e315757a77f851433db1527f48d71b97ec0056174e0a7eb2c467",
               "1f698ed1061bffaf7f380ec5b4ed8075306e6cffd8f6c3cc281873f67003f92d"),
}
CHECKPOINT_FILES = {
    "raw": REFERENCE / "katie/raw_katie.rga",
    "finite": REFERENCE / "truncation/finite_staircase.rga",
    "lower": REFERENCE / "lower/lower_v1.rga",
    "upper24": GUARDED_CAP / "constructed_upper24.rga",
    "upper9": SMALL_CAP,
    "upper474": MATCHED_CAP / "cap_with_matched_boundary.rga",
    "layered": REFERENCE / "boundary/layered_piece.rga",
    "closed": REFERENCE / "closed/closed_candidate.rga",
}


def runtime():
    require(not sys.flags.optimize, "Run without -O/-OO: legacy assertions must remain enabled")
    require(importlib.metadata.version("regina") == "7.4.1", "Requires Regina 7.4.1")
    sys.path.insert(0, str(LEGACY))
    import regina
    return regina


def tight_sha(t):
    return hashlib.sha256(t.tightEncoding().encode()).hexdigest()


def read4(path):
    import regina
    packet = regina.open(str(path))
    require(isinstance(packet, regina.Triangulation4), "Not a 4D triangulation: " + str(path))
    return regina.Triangulation4(packet)


def save4(t, path, label=""):
    import regina
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    require(not path.exists(), "Refusing to overwrite " + str(path))
    packet = regina.PacketOfTriangulation4(t)
    packet.setLabel(label)
    require(packet.save(str(path)), "Could not save " + str(path))


def checkpoint(name, path, out):
    t = read4(path)
    n, expected_file, expected_tight = CHECKPOINTS[name]
    witness = CHECKPOINT_FILES[name]
    require(sha256(witness) == expected_file, name + ": archived witness hash differs")
    def payload(file):
        data = Path(file).read_bytes()
        return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data
    record = {"name": name, "file": str(Path(path).relative_to(out)),
              "pentachora": t.size(), "sha256": sha256(path),
              "tight_sha256": tight_sha(t), "archived_file_sha256": expected_file}
    require(record["pentachora"] == n, name + ": unexpected number of pentachora")
    require(record["tight_sha256"] == expected_tight,
            name + ": exact labelled triangulation differs; no saved replacement will be used")
    actual_payload = payload(path)
    require(actual_payload == payload(witness),
            name + ": uncompressed serialized triangulation differs from the published witness")
    record["uncompressed_payload_sha256"] = hashlib.sha256(actual_payload).hexdigest()
    record["compressed_file_bytes_match_archive"] = record["sha256"] == expected_file
    record["uncompressed_file_bytes_match_archive"] = True
    record["status"] = "PASS_EXACT_PUBLISHED_STATE_AND_SERIALIZED_PAYLOAD"
    write_json(out / "checks" / ("checkpoint_" + name + ".json"), record)
    return record


def step_input(context):
    from orient_katie_diagram import orient
    from run_katie import inspect_input
    out = Path(context["out"])
    source = REFERENCE / "diagram/fig2_published_katie.json"
    require(sha256(source) == PD_SHA, "Marked PD input hash differs")
    original = json.loads(source.read_text())
    # The historical runner used a smaller adapter JSON. Its original bytes
    # are not shipped; construct the same mathematical input with fresh,
    # accurate metadata instead of copying its obsolete source-file digest.
    adapter = {"pd_code": original["pd_code"], "decorations": original["decorations"],
               "component_names": original["order"], "source": "fig2_published_katie.json"}
    adapter_bytes = (json.dumps(adapter, indent=2) + "\n").encode()
    adapter_hash = hashlib.sha256(adapter_bytes).hexdigest()
    data = orient(adapter)
    data["geometry_admission"]["input_sha256"] = adapter_hash
    historical_file = REFERENCE / "katie/input_oriented.json"
    require(sha256(historical_file) == ORIENTED_SHA, "Historical oriented-input binding differs")
    historical = json.loads(historical_file.read_text())
    comparable = json.loads(json.dumps(data))
    comparable["geometry_admission"]["input_sha256"] = historical["geometry_admission"]["input_sha256"]
    require(comparable == historical, "Recomputed marked PD and orientation admission differ")
    require(data["decorations"] == [0, 0, 0, 0, -1, "x", "x", "x"],
            "Unexpected framed/dotted component sequence")
    require(data["geometry_admission"]["reversed_components"] == [5, 7],
            "Unexpected auxiliary orientation reversals")
    folder = out / "katie"
    folder.mkdir()
    (folder / "input_adapter.json").write_bytes(adapter_bytes)
    input_path = folder / "input_oriented.json"
    input_path.write_text(json.dumps(data, indent=2) + "\n")
    admission = inspect_input(data)
    require(not admission["necessary_shape_failures"], "Katie dotted-position admission failed")
    require(admission["counterclockwise_orientation_checked"], "Dotted orientations not admitted")
    require(admission["crossings"] == 188 and admission["component_count"] == 8,
            "Unexpected PD size")
    write_json(folder / "input_checks.json", admission)
    record = {"status": "PASS_MARKED_PD_AND_AUXILIARY_ORIENTATION",
              "source_sha256": PD_SHA, "adapter_sha256": adapter_hash,
              "oriented_sha256": sha256(input_path), "historical_oriented_sha256": ORIENTED_SHA,
              "marked_pd_decorations_and_orientation_checks_match_historical_input": True,
              "historical_adapter_bytes_unavailable_metadata_digest_recomputed": True,
              "crossings": 188, "components": 8, "reversed_components": [5, 7],
              "original_drawing_identification_reproved_here": False}
    write_json(out / "checks/input.json", record)


def step_katie(context):
    import regina
    out = Path(context["out"])
    folder = out / "katie"
    data = json.loads((folder / "input_oriented.json").read_text())
    if context["from_raw"]:
        shutil.copyfile(REFERENCE / "katie/raw_katie.rga", folder / "raw_katie.rga")
        receipt = {"status": "STARTED_FROM_PINNED_RAW_KATIE_OUTPUT",
                   "katie_executed": False, "source": str(REFERENCE / "katie/raw_katie.rga")}
    else:
        executable = Path(context["katie"])
        version_run = subprocess.run([str(executable), "--version"],
                                     capture_output=True, text=True, timeout=30, check=True)
        version = (version_run.stdout + "\n" + version_run.stderr).strip()
        require(re.search(r"(?<![\d.])7\.4\.1(?![\d.])", version) is not None,
                "Katie executable must identify itself as release 7.4.1")
        executable_hash = sha256(executable)
        command = [str(executable), json.dumps(data["pd_code"]),
                   " ".join(str(x) for x in data["decorations"]), "--dim4", "--verbose"]
        started = time.monotonic()
        with (folder / "katie.stdout").open("w") as stdout, \
                (folder / "katie.stderr").open("w") as stderr:
            process = subprocess.run(command, stdout=stdout, stderr=stderr,
                                     timeout=context["katie_timeout"], check=False)
        receipt = {"katie_executed": True, "argv": command, "returncode": process.returncode,
                   "seconds": time.monotonic() - started, "katie_version": version,
                   "katie_sha256": executable_hash,
                   "same_executable_bytes_as_historical_run": executable_hash == HISTORICAL_KATIE_SHA,
                   "included_katie_source_sha256": sha256(ROOT / "legacy/source/katie_7_4_1.cpp"),
                   "provided_binary_source_build_verified_by_this_wrapper": False,
                   "shared_library_path_inherited_from_environment": True}
        write_json(folder / "RUN_RECEIPT.json", receipt)
        require(process.returncode == 0, "Katie failed; inspect katie/katie.stderr")
        require(sha256(executable) == executable_hash, "Katie executable changed while running")
        lines = [line.strip() for line in (folder / "katie.stdout").read_text().splitlines()
                 if line.strip()]
        require(lines, "Katie returned no triangulation signature")
        signature = lines[-1]
        t = regina.Triangulation4.fromIsoSig(signature)
        (folder / "triangulation.isosig").write_text(signature + "\n")
        save4(t, folder / "raw_katie.rga", "input_oriented -- raw Katie output")
        receipt["status"] = "KATIE_REEXECUTED_EXACT_RAW_STATE_REPRODUCED"
    receipt["checkpoint"] = checkpoint("raw", folder / "raw_katie.rga", out)
    receipt["input_sha256"] = sha256(folder / "input_oriented.json")
    receipt["regina_distribution"] = "7.4.1"
    receipt["closed_completion_constructed"] = False
    write_json(folder / "RUN_RECEIPT.json", receipt)


def step_simplify(context):
    from boundary_layer_fast import boundary_from_gluings
    out = Path(context["out"])
    checkpoint("finite", out / "truncation/finite_staircase.rga", out)
    t = read4(out / "truncation/finite_staircase.rga")
    started = time.monotonic()
    changed = t.simplify()
    save4(t, out / "lower/lower_v1.rga")
    record = {"status": "PASS_RECOMPUTED_SINGLE_NATIVE_SIMPLIFY",
              "method": "Regina 7.4.1 Triangulation4.simplify() called exactly once",
              "changed": bool(changed), "seconds": time.monotonic() - started,
              "individual_native_simplification_moves_recorded": False,
              "saved_lower_substituted": False,
              "checkpoint": checkpoint("lower", out / "lower/lower_v1.rga", out)}
    boundary, _, _ = boundary_from_gluings(t)
    require(boundary.size() == 826, "Unexpected simplified lower boundary size")
    record["boundary_tetrahedra"] = boundary.size()
    write_json(out / "lower/SIMPLIFY_RECEIPT.json", record)


def step_cap(context):
    import regina
    from verify_upper_cap_guarded import construct_upper24
    from complete_upper_handles import layer_pachner, cap_audit, boundary_data
    out = Path(context["out"])
    folder = out / "cap"
    folder.mkdir()
    proof = json.loads((GUARDED_CAP / "SOURCE_CONSTRUCTION.json").read_text())
    t = construct_upper24(proof)
    save4(t, folder / "constructed_upper24.rga")
    checkpoint("upper24", folder / "constructed_upper24.rga", out)
    rows = [json.loads(line) for line in (GUARDED_CAP / "moves.jsonl").read_text().splitlines()]
    require(len(rows) == 15, "Expected fifteen guarded cap collapses")
    for i, row in enumerate(rows):
        require(row["index"] == i and row["before_sha256"] == tight_sha(t),
                "Discontinuous cap collapse chain")
        operation = row["operation"]
        require(len(operation) == 3 and operation[:2] == ["C", 1], "Unexpected cap move")
        index = operation[2]
        require(type(index) is int and 0 <= index < t.countEdges(), "Invalid cap edge index")
        require(t.hasCollapseEdge(t.edge(index)), "Cap edge-collapse guard failed")
        require(t.collapseEdge(t.edge(index)), "Cap edge collapse failed")
        require(t.size() == row["pentachora"] and tight_sha(t) == row["after_sha256"],
                "Cap collapse endpoint differs")
    save4(t, folder / "reference_cap9.rga")
    checkpoint("upper9", folder / "reference_cap9.rga", out)
    initial = regina.Triangulation4(t)
    supplied_events = MATCHED_CAP / "layering_events.jsonl.gz"
    with gzip.open(supplied_events, "rt") as stream:
        events = [json.loads(line) for line in stream]
    require(len(events) == 465, "Expected 465 upper boundary moves")
    before = boundary_data(t)[0].isoSig()
    derived_path = {"start": before, "moves": [],
                    "derived_from_recorded_layering_events_sha256": sha256(supplied_events)}
    with gzip.open(folder / "layering_events.jsonl.gz", "wt") as stream:
        for i, event in enumerate(events):
            require(event["step"] == i and event["before"] == before,
                    "Discontinuous upper boundary path")
            actual = layer_pachner(t, before, event["dimension"], event["face"],
                                   expected_after=event["after"], audit=False)
            actual["step"] = i
            for key, value in event.items():
                require(actual.get(key) == value,
                        "Recomputed upper affine attachment differs at move " + str(i) + ": " + key)
            stream.write(json.dumps(actual, sort_keys=True) + "\n")
            derived_path["moves"].append({k: event[k] for k in ("before", "dimension", "face", "after")})
            before = actual["after"]
            if (i + 1) % 100 == 0:
                print("Reconstructed upper boundary moves:", i + 1, flush=True)
    derived_path["end"] = before
    require(before == "gLALQccbefefdpcghfc", "Wrong upper boundary endpoint")
    shellback = regina.Triangulation4(t)
    while shellback.size() > initial.size():
        cell = shellback.simplex(shellback.size() - 1)
        require(shellback.hasShellBoundary(cell), "Upper reverse-shelling guard failed")
        require(shellback.shellBoundary(cell), "Upper reverse shelling failed")
    require(shellback == initial, "Upper shellback did not recover the exact generated nine-cell cap")
    audit = cap_audit(t)
    save4(t, folder / "cap_with_matched_boundary.rga")
    check = checkpoint("upper474", folder / "cap_with_matched_boundary.rga", out)
    write_json(folder / "derived_boundary_path.json", derived_path)
    write_json(folder / "CONSTRUCTION_RECEIPT.json", {
        "status": "PASS_COORDINATE_CAP_RECONSTRUCTION_AND_GUARDED_LAYERING",
        "coordinate_product_reconstructed": True, "guarded_edge_collapses": 15,
        "boundary_moves": 465, "every_recorded_affine_gluing_recomputed": True,
        "reverse_guarded_shellings_recover_generated_cap9": True,
        "saved_cap_substituted": False, "checkpoint": check, "manifold_audit": audit,
        "construction_sha256": sha256(GUARDED_CAP / "SOURCE_CONSTRUCTION.json"),
        "collapse_log_sha256": sha256(GUARDED_CAP / "moves.jsonl"),
        "layering_events_sha256": sha256(supplied_events)})


def step_path(context):
    out = Path(context["out"])
    with gzip.open(REFERENCE / "boundary/path_short.json.gz", "rb") as source:
        data = source.read()
    require(hashlib.sha256(data).hexdigest() == PATH_SHA, "Lower boundary path digest differs")
    path = out / "inputs/lower_boundary_path.json"
    path.parent.mkdir()
    path.write_bytes(data)
    require(len(json.loads(data)["moves"]) == 22304, "Wrong lower boundary path length")
    write_json(out / "checks/boundary_path.json", {
        "status": "PASS_PINNED_RECORDED_BOUNDARY_PATH", "moves": 22304,
        "uncompressed_sha256": PATH_SHA, "search_rerun": False})


def step_finish(context):
    out = Path(context["out"])
    checkpoint("layered", out / "layered_lower/layered_piece.rga", out)
    checkpoint("closed", out / "closed/closed_candidate.rga", out)
    admission = json.loads((out / "admission/CLOSED_MANIFOLD_ADMISSION.json").read_text())
    require(admission["status"] == "PASS_CLOSED_ORIENTABLE_PL_4_MANIFOLD", "Closed admission failed")
    require(admission["source_sha256"] == sha256(out / "closed/closed_candidate.rga")
            and admission["source_tight_sha256"] == CHECKPOINTS["closed"][2],
            "Closed admission binds another state")
    require(admission["f_vector"] == [19, 12136, 48374, 60425, 24170]
            and admission["euler_characteristic"] == 2, "Unexpected closed f-vector")
    write_json(out / "checks/closed_admission.json", {
        "status": "PASS_CLOSED_ADMISSION_BINDS_REBUILT_STATE",
        "f_vector": admission["f_vector"], "euler_characteristic": 2,
        "admission_sha256": sha256(out / "admission/CLOSED_MANIFOLD_ADMISSION.json")})


STEPS = {"input": step_input, "katie": step_katie, "simplify": step_simplify,
         "cap": step_cap, "path": step_path, "finish": step_finish}


def build(args):
    runtime()
    integrity = check_integrity()
    if args.katie:
        executable = args.katie.expanduser().resolve()
        require(executable.is_file() and os.access(executable, os.X_OK),
                "Katie path must name an executable file")
    else:
        executable = None
    require(args.katie_timeout > 0 and args.group_seconds > 0, "Timeouts must be positive")
    out = fresh_output(args.out)
    (out / "checks").mkdir()
    context = {"out": str(out), "from_raw": args.from_raw,
               "katie": str(executable) if executable else None,
               "katie_timeout": args.katie_timeout, "group_seconds": args.group_seconds}
    context_path = out / "BUILD_CONTEXT.json"
    write_json(context_path, context)
    report = {"schema": "motoo-z-construction-v1", "status": "RUNNING",
              "mode": "from_saved_raw_katie" if args.from_raw else "from_marked_pd_with_katie",
              "regina_distribution": "7.4.1", "python": sys.version,
              "released_file_integrity": integrity, "stages": [],
              "standardness_asserted_by_construction": False,
              "source_drawing_identification_reproved_here": False,
              "boundary_searches_rerun": False}
    started = time.monotonic()

    def internal(name):
        run_stage(name, [sys.executable, ROOT / "build.py", "_step", name, context_path], out, report)

    def legacy(name, script, *argv):
        run_stage(name, [sys.executable, LEGACY / script, *argv], out, report)

    try:
        internal("input")
        internal("katie")
        legacy("finite_truncation", "figure2_completion_fast.py", "--source",
               out / "katie/raw_katie.rga", "--out", out / "truncation")
        legacy("independent_truncation", "verify_staircase_truncation.py", "--source",
               out / "katie/raw_katie.rga", "--folder", out / "truncation", "--out",
               out / "checks/independent_truncation.json")
        internal("simplify")
        legacy("independent_cap24_to9", "verify_upper_cap_guarded.py", "--certificate",
               GUARDED_CAP, "--target", SMALL_CAP, "--out", out / "checks/cap24_to9.json")
        internal("cap")
        internal("path")
        legacy("lower_boundary_layering", "boundary_layer_fast.py", "--source",
               out / "lower/lower_v1.rga", "--path", out / "inputs/lower_boundary_path.json",
               "--out", out / "layered_lower")
        legacy("independent_collar", "audit_layering_collar.py", "--source",
               out / "lower/lower_v1.rga", "--path", out / "inputs/lower_boundary_path.json",
               "--layered", out / "layered_lower/layered_piece.rga", "--events",
               out / "layered_lower/layering_events.jsonl.gz", "--out", out / "checks/collar.json")
        legacy("glue_finite_pieces", "boundary_layer_fast.py", "--source",
               out / "layered_lower/layered_piece.rga", "--glue-second",
               out / "cap/cap_with_matched_boundary.rga", "--out", out / "closed")
        legacy("independent_gluing", "audit_fast_closed_gluing.py", "--folder",
               out / "closed", "--out", out / "checks/gluing.json")
        legacy("closed_manifold_admission", "audit_raw_algebra4.py",
               out / "closed/closed_candidate.rga", "--method", "primal", "--out",
               out / "admission", "--group-seconds", str(args.group_seconds))
        internal("finish")
        final_integrity = check_integrity()
        require(final_integrity == integrity, "Release manifest changed during construction")
        report["released_file_integrity_after_run"] = final_integrity
        report["status"] = ("PASS_RECONSTRUCTION_FROM_SAVED_RAW_KATIE" if args.from_raw
                            else "PASS_CONSTRUCTION_FROM_MARKED_PD_WITH_KATIE")
        report["katie_executed"] = not args.from_raw
        report["closed_triangulation"] = "closed/closed_candidate.rga"
        report["closed_file_sha256"] = sha256(out / "closed/closed_candidate.rga")
        report["archived_closed_file_sha256"] = CHECKPOINTS["closed"][1]
        report["closed_tight_sha256"] = CHECKPOINTS["closed"][2]
        report["checkpoints"] = [json.loads(p.read_text()) for p in sorted((out / "checks").glob("checkpoint_*.json"))]
    except BaseException as error:
        report["status"] = "FAILED_CONSTRUCTION"
        report["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        report["seconds"] = time.monotonic() - started
        write_json(out / "BUILD_REPORT.json", report)
    print(report["status"], flush=True)
    print("Closed triangulation:", out / "closed/closed_candidate.rga", flush=True)
    print("Construction is complete; the standardness move replay is a separate command.", flush=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_step":
        require(len(sys.argv) == 4 and sys.argv[2] in STEPS, "Invalid internal construction step")
        runtime()
        context = json.loads(Path(sys.argv[3]).read_text())
        STEPS[sys.argv[2]](context)
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", "--output", dest="out", required=True, type=Path,
                        help="new output directory; must not exist")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--katie", type=Path, help="Katie 7.4.1 executable; shared-library environment is inherited")
    mode.add_argument("--from-raw", action="store_true", help="explicitly skip Katie execution and begin with its pinned raw output")
    parser.add_argument("--katie-timeout", type=float, default=600, help="Katie subprocess budget in seconds (default: 600)")
    parser.add_argument("--group-seconds", type=float, default=120, help="per-vertex-link group search budget (default: 120)")
    build(parser.parse_args())


if __name__ == "__main__":
    main()

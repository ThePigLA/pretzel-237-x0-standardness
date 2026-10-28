#!/usr/bin/env python3
"""Fresh construction audits on isolated copies of the archived inputs.

The default suite checks finite source geometry and the complete X0 assembly.
--native additionally reruns native source transport and manifold admission.
--links replays the portable certificates for all X0 master vertex links.
No downstream 4D contraction or standardization search is run.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time

from common import verify_sources


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source"
REFERENCE = ROOT / "reference/x0"
SOURCE_MANIFEST_SHA256 = "b05c8a1c7af0c0dcf3092fc1d8492fd813aa7e831e487ef44c9116be4a869fe9"
MASTER_SHA256 = "251d7ad38235869fa268e9ee6ec9da2f397cf1ca71db267d0bf0ca9b434c5e4c"
FACETS_SHA256 = "f40b43c09373dff94cde542a3076a8c2f482b1632726b63eed2c5610ddd1e92b"
CANDIDATE_FILES = ("marked_candidate.json.gz", "cover_regina.json.gz",
                   "quotient_regina.json.gz")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def snapshot(folder):
    require(folder.is_dir(), "Missing input directory: " + str(folder))
    answer = {}
    for path in sorted(folder.rglob("*")):
        require(not path.is_symlink(), "Symlink in immutable input tree: " + str(path))
        if path.is_file():
            answer[path.relative_to(folder).as_posix()] = sha(path)
    return answer


def preflight(candidate, native, links):
    package_integrity = verify_sources(ROOT)
    manifest_path = SOURCE / "frozen_release/MANIFEST_SHA256.json"
    require(sha(manifest_path) == SOURCE_MANIFEST_SHA256,
            "Original source manifest identity differs")
    manifest = load(manifest_path)
    for name, row in manifest["files"].items():
        path = (manifest_path.parent / name).resolve()
        require(manifest_path.parent in path.parents, "Manifest path escapes source")
        require(path.stat().st_size == row["bytes"] and sha(path) == row["sha256"],
                "Original source integrity failure: " + name)
    require(sha(REFERENCE / "marked_candidate.json.gz") == MASTER_SHA256,
            "Reference master identity differs")
    files = {}
    for name in CANDIDATE_FILES:
        path = candidate / name
        reference = REFERENCE / name
        actual, expected = load(path), load(reference)
        require(actual == expected,
                "Selected candidate differs from the complete accepted JSON object: " + name)
        files[name] = {"file_sha256": sha(path), "canonical_sha256": digest(actual),
                       "reference_file_sha256": sha(reference),
                       "byte_identical_to_reference": sha(path) == sha(reference)}
    require(load(candidate / "marked_candidate.json.gz")["facets_sha256"] == FACETS_SHA256,
            "Wrong marked X0 facets")
    inherited = []
    for which in ("cover", "quotient"):
        reference_receipt = load(REFERENCE / (which + "_native_admission.json"))
        current_sha = files[which + "_regina.json.gz"]["file_sha256"]
        if not native:
            require(reference_receipt["source_sha256"] == current_sha,
                    "Stored native receipt does not bind selected export bytes; use --native: " + which)
            inherited.append({"which": which,
                              "receipt_sha256": sha(REFERENCE / (which + "_native_admission.json")),
                              "export_sha256": current_sha,
                              "mode": "INHERITED_HASH_BOUND_NATIVE_ADMISSION"})
    if links:
        require(files["marked_candidate.json.gz"]["file_sha256"] == MASTER_SHA256,
                "Portable link certificate requires the exact archived master bytes")
    return {"package_integrity": package_integrity,
            "original_source_payload_files_checked": len(manifest["files"]),
            "original_source_manifest_sha256": SOURCE_MANIFEST_SHA256,
            "selected_candidate_files": files,
            "inherited_native_admissions": inherited}


@dataclass(frozen=True)
class Stage:
    name: str
    purpose: str
    arguments: tuple[str, ...]
    receipts: tuple[str, ...]
    success: tuple[str, object]
    tier: str = "geometry"


def stages(native=False, links=False):
    helper = str(ROOT / "audit_support.py")
    fresh = "_fresh"
    base = [
        Stage("source_coarsening", "Replay all 588 paired source ball replacements and reconstruct the complete descendant",
              (helper, "source-coarsening", "--work", ".", "--receipt", fresh + "/source_coarsening.json"),
              (fresh + "/source_coarsening.json",), ("status", "PASS")),
        Stage("diagram", "Check the signed diagram, spatial knot, puncture and peripheral markings",
              ("frozen_release/diagram/verify_diagram.py", "--receipt", fresh + "/diagram.json"),
              (fresh + "/diagram.json",), ("result", "PASS")),
        Stage("marked_endpoint", "Recompute exact marked exterior isomorphisms from the saved gluing and peripheral arrays",
              ("frozen_release/diagram/verify_marked_isomorphism.py", "--write-certificate", "--certificate", fresh + "/marked_endpoint.json"),
              (fresh + "/marked_endpoint.json",), ("verified", True)),
        Stage("source_geometry", "Check source orientation transport, exact clipping partitions and rational torus charts",
              ("frozen_release/verify_source_geometry.py", "--mutations"),
              ("frozen_release/coarsened/source_orientation.json", "frozen_release/peripheral/cut_orientation.json",
               "frozen_release/peripheral/chart_independent_verification.json"), ("passed", True)),
        Stage("torus_adaptor", "Replay exact torus overlay/collar geometry and all adaptor coarsening steps",
              ("frozen_release/verify_torus_adaptor.py",),
              ("frozen_release/adaptor/independent_verification.json",), ("passed", True)),
        Stage("central_roll", "Check finite roll charts, straightening, inverse maps, full winding and source witnesses",
              ("frozen_release/verify_central_roll.py", "--mutations"),
              ("frozen_release/central_roll/independent_verification.json",), ("passed", True)),
        Stage("mixed_profile", "Check the original simultaneous source-profile transfer and enlarged-neighborhood data",
              ("frozen_release/verify_mixed_profile_transfer.py",),
              ("frozen_release/central_roll/mixed_profile_transfer.json",), ("passed", True)),
        Stage("local_branch_inputs", "Check fixed-circle normal links, branched source vertex links and puncture hypotheses",
              ("frozen_release/verify_twist_localization_inputs.py",),
              ("frozen_release/central_roll/localization_inputs.json",), ("passed", True)),
        Stage("caps", "Check both finite caps, attachment maps and full unwrapped winding",
              ("frozen_release/cap/verify_caps.py", "--output", fresh + "/caps.json"),
              (fresh + "/caps.json",), ("result", "PASS")),
        Stage("bridge", "Check the exact source-to-cap common refinements and relative collars",
              ("frozen_release/cap/verify_bridge.py", "--output", fresh + "/bridge.json"),
              (fresh + "/bridge.json",), ("result", "PASS")),
        Stage("control_transfer", "Check the inherited m=2 control transfer; supplementary to the direct X0 zero-cap argument",
              ("frozen_release/cap/verify_transfer.py", "--output", fresh + "/control_transfer.json"),
              (fresh + "/control_transfer.json",), ("result", "PASS")),
        Stage("localization_review", "Check the finite inputs to the stated relative-PL localization and local branch arguments",
              ("frozen_release/cap/verify_review_inputs.py",),
              ("frozen_release/cap/geometric_review.json",), ("finite_inputs_result", "PASS")),
        Stage("source_interfaces", "Recheck signed source normal orientations, exact peripheral endpoint and interface bindings",
              ("frozen_release/diagram/audit_source_chain.py",),
              ("frozen_release/diagram/source_audit_certificate.json",), ("passed", True)),
        Stage("zero_cap", "Check all m=0 spherical identity cells and compare the actual marked master",
              ("review/check_zero_cap_geometry.py", "--candidate", "candidate/marked_candidate.json.gz", "--output", fresh + "/zero_cap.json"),
              (fresh + "/zero_cap.json",), ("result", "PASS")),
    ]
    if native:
        base.extend([
            Stage("native_diagram", "Reconstruct the explicit rational-tangle diagram with pinned Spherogram",
                  ("frozen_release/diagram/verify_diagram.py", "--with-library", "--receipt", fresh + "/native_diagram.json"),
                  (fresh + "/native_diagram.json",), ("result", "PASS"), "native"),
            Stage("native_source_transport", "Rerun native source peripheral transport and independently check its fresh marked endpoint",
                  (helper, "native-transport", "--work", ".", "--receipt", fresh + "/native_source_transport.json"),
                  (fresh + "/native_source_transport.json", fresh + "/native_transport_endpoint.json"),
                  ("status", "PASS"), "native"),
            Stage("native_source_manifolds", "Rebuild and admit the original 3D source cover and quotient; recognize the quotient as S3",
                  ("frozen_release/seifert/native_admit.py",),
                  ("frozen_release/seifert/data/native_admission.json",), ("cases", None), "native"),
            Stage("native_cover", "Fresh native manifold admission of the exact X0 cover export",
                  ("scripts/admit_native_candidate.py", "cover", "--output", fresh + "/cover_native_admission.json"),
                  (fresh + "/cover_native_admission.json",), ("passed", True), "native"),
            Stage("native_quotient", "Fresh native manifold admission of the exact X0 quotient export",
                  ("scripts/admit_native_candidate.py", "quotient", "--output", fresh + "/quotient_native_admission.json"),
                  (fresh + "/quotient_native_admission.json",), ("passed", True), "native"),
        ])
    base.append(Stage("x0_assembly", "Independently regenerate all marked X0 assembly cells and check every cover/quotient face square",
                      ("scripts/verify_candidate_assembly.py", "--output", fresh + "/x0_assembly.json"),
                      (fresh + "/x0_assembly.json",), ("passed", True)))
    if links:
        base.append(Stage("vertex_links", "Replay all 4310 portable master vertex-link certificates",
                          ("portable_links/replay_x0_links.py", "--receipt", fresh + "/vertex_links.json"),
                          (fresh + "/vertex_links.json",), ("status", "PASS"), "links"))
    return base


def environment(native):
    result = {"python": platform.python_version(), "implementation": platform.python_implementation(),
              "platform": platform.platform(), "machine": platform.machine(),
              "subprocess_python_flags": ["-E", "-B"], "assertions_enabled": True}
    packages = {}
    for name in ("regina", "snappy", "spherogram", "sympy", "psutil"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    result["installed_distributions"] = packages
    if native:
        for name, version in (("regina", "7.4.1"), ("snappy", "3.3.2"),
                              ("spherogram", "2.4.1"), ("sympy", "1.14.0")):
            require(packages[name] == version,
                    "Native audit requires " + name + "==" + version)
    return result


def run_process(argv, cwd, log, timeout):
    started = time.monotonic()
    process = None
    env = dict(os.environ)
    env.pop("PYTHONOPTIMIZE", None)
    try:
        with log.open("w") as handle:
            process = subprocess.Popen(argv, cwd=cwd, env=env,
                                       stdout=handle, stderr=subprocess.STDOUT,
                                       start_new_session=(os.name == "posix"))
            while True:
                try:
                    code = process.wait(timeout=20)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - started
                    if elapsed > timeout:
                        raise TimeoutError("Construction stage exceeded its time limit")
                    print("RUNNING", log.parent.name, round(elapsed), "seconds", flush=True)
    except BaseException:
        if process is not None and process.poll() is None:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGTERM)
            else:
                process.terminate()
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
                process.wait()
        raise
    require(code == 0, "Stage failed; see " + str(log))
    return time.monotonic() - started


def validate_receipts(stage, work):
    first = load(work / stage.receipts[0])
    if stage.name == "native_source_manifolds":
        for which in ("cover", "quotient"):
            case = first["cases"][which]
            require(all(case[k] is True for k in ("valid", "closed", "orientable", "connected")),
                    "Native source admission failed: " + which)
            require(case["H1_on_unmarked_simplified_copy"] == "0",
                    "Native source first homology differs: " + which)
        require(first["cases"]["quotient"]["S3_recognition"] is True,
                "Native source quotient was not recognized as S3")
    else:
        key, value = stage.success
        require(first.get(key) == value, "Stage did not produce a successful receipt: " + stage.name)
    if stage.name == "source_geometry":
        require(all(load(work / p).get("passed") is True for p in stage.receipts),
                "A source-geometry subcheck failed")
    if stage.name == "marked_endpoint":
        require(first == load(work / "frozen_release/diagram/marked_isomorphism_certificate.json"),
                "Recomputed marked endpoint certificate differs from the archived exact certificate")
    if stage.name == "native_source_transport":
        require(load(work / stage.receipts[1])["verified"] is True,
                "Fresh native transport endpoint certificate failed")
    if stage.name == "native_diagram":
        require(first["optional_rational_tangle_regeneration"] == "PASS",
                "Native diagram reconstruction was not performed")
    return first


def run_stage(stage, output, candidate, fresh_native, timeout, keep_work):
    folder = output / "stages" / stage.name
    folder.mkdir(parents=True)
    work = folder / "work"
    shutil.copytree(SOURCE, work)
    shutil.copytree(REFERENCE, work / "candidate")
    for name in CANDIDATE_FILES:
        shutil.copy2(candidate / name, work / "candidate" / name)
    for which, receipt in fresh_native.items():
        shutil.copy2(receipt, work / "candidate" / (which + "_native_admission.json"))
    (work / "_fresh").mkdir()
    # Some legacy scripts write to fixed receipt locations. Remove each target
    # in the disposable copy so an old PASS cannot masquerade as a fresh one.
    for name in stage.receipts:
        target = work / name
        if target.exists():
            target.unlink()
    argv = [sys.executable, "-E", "-B", *stage.arguments]
    log = folder / "stage.log"
    row = {"name": stage.name, "purpose": stage.purpose, "tier": stage.tier,
           "status": "RUNNING", "executed": True,
           "command": ["python3", "-E", "-B", *[
               "audit_support.py" if x == str(ROOT / "audit_support.py") else x
               for x in stage.arguments]],
           "worktree_policy": "Pristine isolated copy; original source/reference files are never edited"}
    atomic_json(folder / "stage.json", row)
    try:
        row["wall_seconds"] = run_process(argv, work, log, timeout)
        validate_receipts(stage, work)
        receipt_folder = folder / "receipts"
        receipt_folder.mkdir()
        row["receipts"] = []
        for name in stage.receipts:
            path = work / name
            target = receipt_folder / path.name
            require(not target.exists(), "Two stage receipts have the same basename")
            shutil.copy2(path, target)
            row["receipts"].append({"path": target.relative_to(output).as_posix(),
                                     "sha256": sha(target), "recomputed_this_run": True})
        if stage.name in ("native_cover", "native_quotient"):
            fresh_native[stage.name.removeprefix("native_")] = output / row["receipts"][0]["path"]
        row["status"] = "PASS"
        if stage.name == "x0_assembly":
            row["manifold_admission_evidence"] = {
                which: "FRESH_NATIVE_REPLAY_THIS_RUN" if which in fresh_native
                else "INHERITED_HASH_BOUND_NATIVE_RECEIPT"
                for which in ("cover", "quotient")}
        if not keep_work:
            shutil.rmtree(work)
        row["worktree_retained"] = keep_work
    except BaseException as exc:
        row["status"] = "FAILED"
        row["error"] = str(exc)
        row["worktree_retained"] = True
        raise
    finally:
        row["log"] = log.relative_to(output).as_posix()
        row["log_sha256"] = sha(log) if log.exists() else None
        atomic_json(folder / "stage.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="New or empty directory for reports, logs and disposable working copies")
    parser.add_argument("--candidate", type=Path, default=REFERENCE,
                        help="Directory containing the three exact generated X0 master/export JSON.gz files; default reference/x0")
    parser.add_argument("--native", action="store_true", help="Also rerun pinned native source transport and source/X0 manifold admission")
    parser.add_argument("--links", action="store_true", help="Also replay every portable X0 master vertex-link certificate")
    parser.add_argument("--stage", action="append", default=[],
                        help="Run named stage(s) only; the report explicitly marks the audit as partial")
    parser.add_argument("--timeout", type=float, default=900,
                        help="Per-stage wall-time limit in seconds (default900)")
    parser.add_argument("--keep-work", action="store_true", help="Keep successful disposable worktrees as well as reports")
    parser.add_argument("--plan", action="store_true", help="List selected stages without executing or writing")
    args = parser.parse_args()
    require(not sys.flags.optimize, "Run audit.py without -O: assertions are proof checks")
    require(args.timeout > 0, "Per-stage timeout must be positive")
    available = stages(args.native, args.links)
    names = {s.name for s in available}
    require(set(args.stage) <= names,
            "Unknown or disabled stage(s): " + ", ".join(sorted(set(args.stage) - names)))
    selected = [s for s in available if not args.stage or s.name in args.stage]
    require(selected, "No audit stages selected")
    if args.plan:
        print(json.dumps({"scope": "SELECTED_STAGES" if args.stage else "COMPLETE_REQUESTED_SUITE",
                          "stages": [{"name": s.name, "tier": s.tier, "purpose": s.purpose}
                                     for s in selected]}, indent=2))
        return
    candidate = args.candidate.expanduser().resolve()
    output = args.output.expanduser().resolve()
    for protected in (ROOT, SOURCE, REFERENCE, candidate):
        require(output != protected and output not in protected.parents,
                "Output would contain an input tree")
    for protected in (SOURCE, ROOT / "reference", ROOT / "provenance", candidate):
        require(protected not in output.parents, "Output is inside an immutable input tree")
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())),
            "Output directory must be new or empty")
    before = {"source": snapshot(SOURCE), "reference": snapshot(ROOT / "reference")}
    entry_points = {name: sha(ROOT / name)
                    for name in ("audit.py", "audit_support.py", "common.py")}
    bound = preflight(candidate, args.native, args.links)
    env = environment(args.native)
    if args.native and args.stage and "x0_assembly" in args.stage:
        require({"native_cover", "native_quotient"} <= set(args.stage),
                "A selected --native assembly audit also requires --stage native_cover and --stage native_quotient")
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    report = {"schema": "x0-construction-audit-v1", "status": "RUNNING",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "SELECTED_STAGES" if args.stage else "COMPLETE_REQUESTED_SUITE",
              "requested_native": args.native, "requested_links": args.links,
              "candidate_source": "reference/x0" if candidate == REFERENCE.resolve() else "user-supplied generated candidate",
              "environment": env, "input_bindings": bound,
              "entry_point_file_sha256": entry_points,
              "input_inventory_sha256": digest(before),
              "requested_stages": [s.name for s in selected], "stages": [],
              "inherited_evidence_and_mathematical_premises": [
                  "Recorded compact-source and peripheral data are exact inputs; the default endpoint stage rechecks their saved marked comparison without rerunning native transport.",
                  "Native diagram constructors, finite-vertex removal and peripheral move semantics remain trusted when --native reruns them; no elementary native3D move transcript is claimed.",
                  "The stated relative PL neighborhood, product, collar, cone and gluing arguments retain their mathematical meaning; these programs do not export a global pointwise conjugacy.",
                  "Default assembly admission binds stored native cover/quotient receipts. --native replaces these by fresh admissions; --links additionally proves the cover's vertex-link premise by portable moves."
              ],
              "excluded_work": ["4D contraction traces", "fundamental-group elimination", "terminal standardization", "randomized simplification searches"]}
    summary_path = output / "audit_report.json"
    atomic_json(summary_path, report)
    fresh_native = {}
    try:
        for stage in selected:
            print("CHECK", stage.name, flush=True)
            try:
                row = run_stage(stage, output, candidate, fresh_native,
                                args.timeout, args.keep_work)
            except BaseException:
                stage_report = output / "stages" / stage.name / "stage.json"
                if stage_report.exists():
                    report["stages"].append(load(stage_report))
                raise
            report["stages"].append(row)
            atomic_json(summary_path, report)
            print("PASS", stage.name, round(row["wall_seconds"], 2), "seconds", flush=True)
        after = {"source": snapshot(SOURCE), "reference": snapshot(ROOT / "reference")}
        require(after == before, "An immutable source/reference file changed during the audit")
        require(all(sha(ROOT / name) == expected for name, expected in entry_points.items()),
                "An audit entry point changed during the audit")
        for name, row in bound["selected_candidate_files"].items():
            require(sha(candidate / name) == row["file_sha256"],
                    "Selected candidate changed during the audit: " + name)
        report["immutable_inputs_unchanged"] = True
        report["status"] = "PASS"
    except BaseException as exc:
        report["status"] = "INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "FAILED"
        report["error"] = str(exc)
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - started
        atomic_json(summary_path, report)
    print("PASS:", report["scope"], "with", len(report["stages"]), "fresh stages", flush=True)
    print("Report:", summary_path, flush=True)


if __name__ == "__main__":
    main()

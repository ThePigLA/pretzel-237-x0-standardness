#!/usr/bin/env python3
"""Package the selected Motoo-readiness checkpoint; no topology is recomputed.

Run --check-selection before building.  Pass completed additional result
directories with repeated --extra-results PATH arguments.  Paths must remain
inside this project, and result directories with a best_certificate.json must
also contain their finished REPLAY.json receipt.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


# Keep the earlier conversion checkpoint intact.  Selection is deliberate:
# native binaries, environments, abandoned searches and duplicate cap builds
# do not belong to this portable source-and-certificate archive.
BASELINE = (
    "README.md", "NORMALIZATION.md", "requirements.txt", "code", "inputs", "source",
    "results/figure1_handlebody_group.json", "results/figure2_handlebody_group.json",
    "results/fresh_input_verification_v1",
    "results/boundary/figure1_v1", "results/boundary/figure2_v1",
    "results/boundary/selftest_v1",
    "results/boundary/INDEPENDENT_GROUP_EXTRACTION_REVIEW.json",
    "results/boundary/upper_cap_selftest_v3",
    "results/boundary/reference_walk_seed0/SUCCESS.json",
    "results/boundary/standard_cap_six_boundary_v1",
    "results/boundary/upper_cap_double_fixture_v1",
    "scratch/diagram/complete_intersection_audit.json",
    "scratch/diagram/figure1_annotated.png", "scratch/diagram/figure2_annotated.png",
    "scratch/diagram/figure1_weighted_pickup.json",
    "scratch/diagram/figure1_weighted_then_r3.json",
    "scratch/diagram/figure1_pickup_walk.json",
    "scratch/diagram/figure1_final_bounded_pickup.json",
    "scratch/diagram/figure1_remaining_dotted_R2_targets.json",
    "scratch/diagram/figure2_weighted_pickup.json",
    "scratch/katie/controls_admitted",
    "scratch/katie/DEPENDENCY_RECEIPT.json", "scratch/katie/EXECUTABLE_CHECK.json",
    "scratch/katie/PACKAGE_RECEIPT.json",
    "scratch/katie/figure1_final_CHAIN_CHECK.json",
    "scratch/katie/figure2_combined_CHAIN_CHECK.json",
)
READINESS = (
    "READINESS.md",
    "results/readiness/intake_v1",
    "results/readiness/runner_admission_v1",
    "results/readiness/INDEPENDENT_READINESS_REVIEW.json",
)
CONTROL = "results/readiness/completion_control_v1"
CONTROL_FILES = (
    "README.md", "connected_input.json", "katie_connected",
    "truncated_raw.rga", "lower_reference.rga", "LOWER_REFERENCE.json",
    "lower_boundary.rga", "selected_boundary_path.json", "lower_layered",
    "closed", "COMPLETION_CONTROL.json",
    # Retain the historical diagnostic that prompted projection admission.
    "input.json", "katie",
)
EXCLUDED_TREES = (
    "results/readiness/intake_v1/selftest",
    "results/readiness/completion_control_reduction_v1",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text())


def select(root, extra_results, include_exploration_streams=False):
    selected = set()
    included_extra = []

    def include(relative):
        path = Path(relative)
        if not path.is_absolute():
            path = root / path
        path = path.resolve()
        require(path.is_relative_to(root), "Selection escapes project: " + str(path))
        if not path.exists():
            raise FileNotFoundError(path)
        candidates = path.rglob("*") if path.is_dir() else [path]
        for candidate in candidates:
            if not candidate.is_file():
                continue
            relative = candidate.relative_to(root)
            name = relative.as_posix()
            if "__pycache__" in relative.parts or candidate.suffix == ".pyc":
                continue
            # The chosen certificate embeds every before/after state and move
            # needed for its proof. Exploration streams additionally contain
            # discarded walks and restarts; omit that diagnostic bulk by default.
            if candidate.name == "accepted_moves.jsonl.gz" and not include_exploration_streams:
                continue
            if any(name == p or name.startswith(p + "/") for p in EXCLUDED_TREES):
                continue
            # Only selected_boundary_path.json is needed for the control.
            if name.startswith(CONTROL + "/boundary_path/") or name.startswith(CONTROL + "/boundary_seed"):
                continue
            require(not candidate.is_symlink(), "Symlink is not a portable package input: " + name)
            selected.add(candidate)

    for relative in BASELINE + READINESS:
        include(relative)
    for relative in CONTROL_FILES:
        include(CONTROL + "/" + relative)

    # Preserve the original static review and any later review/regression.
    for path in sorted((root / "results/readiness").glob("INDEPENDENT*REVIEW*.json")):
        include(path)
    for followup in sorted((root / "results/readiness").glob("tight_replay_followup_v*")):
        include(followup)
    for followup in sorted((root / "results/readiness").glob("INTAKE_METADATA_FOLLOWUP*.json")):
        include(followup)

    # Reduction attempts are opt-in, including partial or superseded attempts.
    # Never silently package an in-progress attempt because its folder exists.
    for supplied in extra_results:
        path = Path(supplied)
        if not path.is_absolute():
            path = root / path
        path = path.resolve()
        require(path.is_relative_to(root / "results"), "--extra-results must be inside project results/")
        require(path.is_dir(), "--extra-results is not a directory: " + str(path))
        relative = path.relative_to(root).as_posix()
        require(not any(relative == p or relative.startswith(p + "/") for p in EXCLUDED_TREES),
                "Explicitly excluded result directory: " + relative)
        if (path / "best_certificate.json").is_file():
            require((path / "REPLAY.json").is_file(),
                    "Reduction has no completed replay receipt: " + relative)
        include(path)
        included_extra.append(relative)
    omitted = sorted(p.relative_to(root).as_posix()
                     for p in (root / "results/readiness").glob("completion_control_reduction_*")
                     if p.is_dir() and p.relative_to(root).as_posix() not in included_extra)
    return selected, sorted(set(included_extra)), omitted


def validate_receipt_bindings(root, selected, extra_results):
    """Only file bindings and receipt consistency; this is not a new replay."""
    for figure, relative in ((1, "scratch/katie/figure1_final_CHAIN_CHECK.json"),
                             (2, "scratch/katie/figure2_combined_CHAIN_CHECK.json")):
        receipt = read_json(root / relative)
        require(receipt["status"] == "PASS", f"Figure {figure} normalization has not passed")
        require(sha256(root / f"inputs/figure{figure}_link.json") == receipt["original_source_sha256"],
                "Original normalization source changed")
        for segment in receipt["segments"]:
            path = Path(segment["file"])
            if not path.is_absolute():
                path = root / path
            path = path.resolve()
            require(path in selected, "Checked normalization segment omitted: " + str(path))
            require(sha256(path) == segment["sha256"], "Normalization segment changed after replay")

    control = root / CONTROL
    audit = read_json(control / "COMPLETION_CONTROL.json")
    require(audit["status"] == "PASS_EXPLICIT_GENUS_TWO_COMPLETION_CONTROL", "Control audit did not pass")
    require(audit["checker_sha256"] == sha256(root / "code/verify_completion_control.py"),
            "Completion checker changed after its recorded audit")
    cap = root / "results/boundary/standard_cap_six_boundary_v1/cap_with_matched_boundary.rga"
    require(audit["standard_cap_sha256"] == sha256(cap), "Standard cap changed after control audit")
    preparation = read_json(control / "LOWER_REFERENCE.json")
    require(preparation["source_truncated_sha256"] == sha256(control / "truncated_raw.rga"),
            "Truncated source changed")
    require(preparation["result_reference_sha256"] == sha256(control / "lower_reference.rga"),
            "Saved lower reference changed")
    require(preparation["individual_move_log_available"] is False,
            "Initial native simplification trust boundary unexpectedly changed")

    reductions = []
    for relative in extra_results:
        path = root / relative
        if not (path / "best_certificate.json").is_file():
            continue
        certificate = read_json(path / "best_certificate.json")
        replay = read_json(path / "REPLAY.json")
        require(replay["status"] in ("PASS_PARTIAL_PATH", "PASS_STANDARDIZATION"),
                "Reduction replay did not pass: " + relative)
        require(replay["moves"] == len(certificate["moves"]), "Replay/certificate move counts differ")
        require(replay["standard_endpoint"] == certificate["standard_endpoint"],
                "Replay/certificate standardness flags differ")
        require(replay["source_binding_checked"] is True, "Reduction replay did not bind a source file")
        sources = [p for p in selected if p.name == Path(certificate["source_file"]).name]
        require(any(sha256(p) == certificate["source_sha256"] for p in sources),
                "Reduction source file is missing or changed: " + relative)
        reductions.append({"directory": relative, "status": replay["status"],
                           "moves": replay["moves"], "source_binding_checked": True})
    return {"normalization_receipt_file_bindings_checked": True,
            "completion_receipt_file_bindings_checked": True,
            "completed_reduction_receipts": reductions,
            "topological_audits_rerun_by_builder": False}


def build(root, target, selected, details):
    target.parent.mkdir(parents=True, exist_ok=True)
    require(not target.exists(), "Output already exists; choose a new archive filename")
    with tempfile.TemporaryDirectory(prefix="motoo_readiness_") as temporary:
        stage_parent = Path(temporary)
        stage = stage_parent / "motoo_kirby_v1"
        manifest = []
        for original in sorted(selected):
            relative = original.relative_to(root)
            destination = stage / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            before = sha256(original)
            shutil.copy2(original, destination)
            require(sha256(destination) == before == sha256(original),
                    "Input changed during archive staging: " + relative.as_posix())
            manifest.append(f"{before}  {relative.as_posix()}")
        (stage / "SHA256SUMS").write_text("\n".join(manifest) + "\n")
        pending_zip = stage_parent / "readiness.zip"
        with zipfile.ZipFile(pending_zip, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=9, allowZip64=True) as archive:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(stage_parent).as_posix())
        extracted = stage_parent / "extracted"
        with zipfile.ZipFile(pending_zip) as archive:
            bad = archive.testzip()
            require(bad is None, "ZIP CRC failure: " + str(bad))
            require(len(archive.namelist()) == len(selected) + 1, "Unexpected ZIP entry count")
            archive.extractall(extracted)
        project = extracted / "motoo_kirby_v1"
        subprocess.run([sys.executable, "-O", str(project / "code/check_archive_hashes.py")], check=True)
        extracted_files = {p.relative_to(project).as_posix() for p in project.rglob("*") if p.is_file()}
        expected_files = {p.relative_to(root).as_posix() for p in selected} | {"SHA256SUMS"}
        require(extracted_files == expected_files, "Extracted inventory differs from selection")
        # Only publish the artifact after both archive and extracted-file checks.
        shutil.copy2(pending_zip, target)
    receipt = {**details, "file": str(target), "bytes": target.stat().st_size,
               "sha256": sha256(target), "zip_crc_passed": True,
               "extracted_hashes_passed": True, "extracted_inventory_passed": True,
               "status": "CHECKED_READINESS_ARCHIVE_NOT_A_MOTOO_STANDARDNESS_CERTIFICATE"}
    target.with_suffix(".receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="New output ZIP; required unless --check-selection")
    parser.add_argument("--extra-results", type=Path, action="append", default=[],
                        help="Completed project result directory; may be repeated")
    parser.add_argument("--check-selection", action="store_true", help="Check selected files without creating an archive")
    parser.add_argument("--include-exploration-streams", action="store_true",
                        help="Also package raw search streams, beyond complete selected path certificates")
    args = parser.parse_args()
    if args.out is None and not args.check_selection:
        parser.error("--out is required unless --check-selection is used")
    root = Path(__file__).resolve().parents[1]
    selected, extras, omitted = select(root, args.extra_results, args.include_exploration_streams)
    details = {"distributed_files_in_manifest": len(selected),
               "uncompressed_bytes": sum(p.stat().st_size for p in selected),
               "extra_results": extras, "unselected_reduction_directories_omitted": omitted,
               "exploration_streams_included": args.include_exploration_streams,
               **validate_receipt_bindings(root, selected, extras)}
    if args.check_selection:
        print(json.dumps({"status": "PASS_SELECTION_ONLY_NO_ARCHIVE_CREATED", **details}, indent=2))
    else:
        print(json.dumps(build(root, args.out.resolve(), selected, details), indent=2))


if __name__ == "__main__":
    main()

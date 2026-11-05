#!/usr/bin/env python3
"""Create the checked research-checkpoint archive from the working directory."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = args.out.resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    selected = set()

    def include(relative):
        path = root / relative
        if not path.exists():
            raise FileNotFoundError(path)
        candidates = path.rglob("*") if path.is_dir() else [path]
        selected.update(p for p in candidates if p.is_file()
                        and "__pycache__" not in p.parts and p.suffix != ".pyc")

    for relative in [
        "README.md", "NORMALIZATION.md", "requirements.txt", "code", "inputs", "source",
        "results/figure1_handlebody_group.json", "results/figure2_handlebody_group.json",
        "results/fresh_input_verification_v1",
        "results/boundary/figure1_v1", "results/boundary/figure2_v1",
        "results/boundary/selftest_v1", "results/boundary/INDEPENDENT_GROUP_EXTRACTION_REVIEW.json",
        "results/boundary/upper_cap_selftest_v3", "results/boundary/reference_walk_seed0/SUCCESS.json",
        "results/boundary/standard_cap_six_boundary_v1", "results/boundary/upper_cap_double_fixture_v1",
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
    ]:
        include(relative)

    # Do not silently package a producer-only normalization endpoint.
    for figure, relative in [(1, "scratch/katie/figure1_final_CHAIN_CHECK.json"),
                             (2, "scratch/katie/figure2_combined_CHAIN_CHECK.json")]:
        receipt = json.loads((root / relative).read_text())
        if receipt["status"] != "PASS":
            raise ValueError(f"Figure {figure} normalization chain has not passed")
        source = root / f"inputs/figure{figure}_link.json"
        if hashlib.sha256(source.read_bytes()).hexdigest() != receipt["original_source_sha256"]:
            raise ValueError("Original normalization source mismatch")
        for segment in receipt["segments"]:
            path = Path(segment["file"])
            if not path.is_absolute():
                path = root / path
            if hashlib.sha256(path.read_bytes()).hexdigest() != segment["sha256"]:
                raise ValueError("Normalization segment changed after its replay")
            if path not in selected:
                raise ValueError("Checked normalization segment omitted from archive")

    with tempfile.TemporaryDirectory(prefix="motoo_distribution_") as temporary:
        stage_parent = Path(temporary)
        stage = stage_parent / "motoo_kirby_v1"
        for original in sorted(selected):
            destination = stage / original.relative_to(root)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, destination)
        manifest = []
        for path in sorted(stage.rglob("*")):
            if path.is_file():
                manifest.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(stage).as_posix()}")
        (stage / "SHA256SUMS").write_text("\n".join(manifest) + "\n")
        with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=9, allowZip64=True) as archive:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(stage_parent).as_posix())
        extracted = stage_parent / "extracted"
        with zipfile.ZipFile(target) as archive:
            bad = archive.testzip()
            if bad is not None:
                raise ValueError("ZIP CRC check failed: " + bad)
            archive.extractall(extracted)
        project = extracted / "motoo_kirby_v1"
        subprocess.run([sys.executable, str(project / "code/check_archive_hashes.py")], check=True)
        # These proof checks only require the standard library.
        for figure in (1, 2):
            subprocess.run([sys.executable, "-O", str(project / "code/check_handlebody_group.py"),
                            "replay", str(project / f"results/figure{figure}_handlebody_group.json")], check=True)
    result = {"file": str(target), "bytes": target.stat().st_size,
              "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
              "distributed_files_in_manifest": len(selected),
              "zip_crc_passed": True, "extracted_hashes_passed": True,
              "extracted_standard_library_group_replays_passed": True,
              "status": "VERIFIED_INPUT_AND_CONVERSION_CHECKPOINT_NOT_COMPLETED_4D_CONVERSION"}
    target.with_suffix(".receipt.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

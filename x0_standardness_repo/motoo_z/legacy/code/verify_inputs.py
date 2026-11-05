#!/usr/bin/env python3
"""Rebuild the PD transcription, pi1 certificates and 3D boundary checks.

This verifies the input/admission phase of the Kirby conversion project.
It deliberately does not declare a closed 4D triangulation or standardness.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("verification_fresh"))
    parser.add_argument("--skip-pdf-rebuild", action="store_true",
                        help="check shipped PD inputs rather than re-extracting the PDF")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.out.resolve()
    require(not output.exists(), "choose a new output directory")
    output.mkdir(parents=True)
    source = root / "inputs/motoo_rollspun.pdf"
    receipt = json.loads((root / "inputs/source_receipt.json").read_text())
    require(sha(source) == receipt["sha256"], "source PDF checksum mismatch")
    stages = []

    def run(label, script, *extra):
        command = [sys.executable, "-E", "-B", str(root / "code" / script), *map(str, extra)]
        logfile = output / f"{label}.log"
        print(label, flush=True)
        with logfile.open("w") as stream:
            subprocess.run(command, cwd=root, stdout=stream, stderr=subprocess.STDOUT, check=True)
        stages.append({"label": label, "log": logfile.name, "sha256": sha(logfile)})

    if not args.skip_pdf_rebuild:
        run("diagram_rebuild", "extract_diagrams.py", "--pdf", source,
            "--out", output / "inputs", "--scratch", output / "diagram",
            "--link-python", sys.executable)
        for figure in (1, 2):
            rebuilt = json.loads((output / f"inputs/figure{figure}_link.json").read_text())
            original = json.loads((root / f"inputs/figure{figure}_link.json").read_text())
            for key in ("pd_code", "regina_pd_code", "pd_component_to_source", "decorations",
                        "linking_matrix", "crossing_count", "component_count"):
                require(rebuilt[key] == original[key], f"Figure {figure}: rebuilt {key} differs")

    for figure in (1, 2):
        source_link = root / f"inputs/figure{figure}_link.json"
        shipped_group = root / f"results/figure{figure}_handlebody_group.json"
        shipped = json.loads(shipped_group.read_text())
        require(shipped["input_sha256"] == sha(source_link), "stored group/input binding mismatch")
        run(f"figure{figure}_saved_group_replay", "check_handlebody_group.py", "replay", shipped_group)
        fresh_group = output / f"figure{figure}_group.json"
        run(f"figure{figure}_fresh_group", "check_handlebody_group.py", "compute", source_link, fresh_group)
        group = json.loads(fresh_group.read_text())
        require(group["status"] == "TRIVIAL_HANDLEBODY_PI1", "fresh group proof did not finish")
        run(f"figure{figure}_fresh_boundary", "check_boundary.py", source_link,
            "--out", output / f"figure{figure}_boundary", "--expected", 2)
        boundary = json.loads((output / f"figure{figure}_boundary/BOUNDARY_VERIFICATION.json").read_text())
        require(boundary["status"] == "PASS_EXPECTED_BOUNDARY", "wrong boundary result")
        require(boundary["conditional_closed_handle_counts"] == [1, 3, 5, 2, 1], "handle counts")
        require(boundary["conditional_closed_euler_characteristic"] == 2, "closed Euler characteristic")

    summary = {
        "status": "PASS_DIAGRAM_INPUTS_AND_HOMOTOPY_SPHERE_PREMISES",
        "source_pdf_sha256": sha(source),
        "pdf_transcription_rebuilt": not args.skip_pdf_rebuild,
        "checks": stages,
        "closed_4d_triangulation_constructed_by_this_verifier": False,
        "standardness_proved_by_this_verifier": False,
        "qualification": "The diagram is interpreted with ordinary integral Kirby framings and the usual two-3-handle/one-4-handle completion. The PD transcription is explicitly checked against the supplied vector PDF."
    }
    (output / "VERIFICATION.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": summary["status"], "output": str(output)}), flush=True)


if __name__ == "__main__":
    main()

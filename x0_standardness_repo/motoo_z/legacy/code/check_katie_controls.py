#!/usr/bin/env python3
"""Check positive Kirby controls and an ideal-boundary negative control."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import regina


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path,
                        help="New output directory")
    parser.add_argument("--katie", type=Path,
                        help="Path to a pinned Katie 7.4.1 executable (for example in Regina.app)")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    runner = Path(__file__).with_name("run_katie.py")
    cases = [
        # Exercise both zero-based SnapPy and one-based Regina PD labels.
        ("cp2", [[0, 0, 1, 1]], [1]),
        ("zero_unknot", [[1, 1, 2, 2]], [0]),
        ("cancelling_pair", [[3, 2, 4, 1], [1, 4, 2, 3]], ["x", 0]),
    ]
    reports = []
    for name, pd, decorations in cases:
        input_path = args.out / (name + ".json")
        data = {"name": name, "pd_code": pd, "decorations": decorations}
        if "x" in decorations:
            from orient_katie_diagram import orient
            data = orient(data)
        input_path.write_text(json.dumps(data, indent=2) + "\n")
        command = [sys.executable, str(runner), str(input_path),
                   "--out", str(args.out / name)]
        if args.katie is not None:
            command.extend(["--katie", str(args.katie)])
        subprocess.run(command, check=True)
        report = json.loads((args.out / name / "RUN_RECEIPT.json").read_text())
        tri = report["triangulation"]
        if name == "cp2":
            require(tri["closed_manifold_admission"] and
                    tri["euler_characteristic_manifold"] == 3 and
                    tri["homology"] == {"1": "0", "2": "Z"},
                    "The +1-framed-unknot control did not produce CP2 invariants")
        elif name == "zero_unknot":
            require(tri["ideal"] and not tri["closed"] and
                    tri["euler_characteristic_manifold"] == 2 and
                    tri["homology"] == {"1": "0", "2": "Z"},
                    "The zero-framed-unknot control incorrectly became a closed sphere")
            require(sum(not x["sphere"] for x in tri["vertex_links"]) == 1,
                    "Expected one non-spherical ideal vertex")
        else:
            require(tri["closed_manifold_admission"] and
                    tri["euler_characteristic_manifold"] == 2 and
                    tri["homology"] == {"1": "0", "2": "0"},
                    "The cancelling 1/2-pair control did not produce sphere invariants")
            t = regina.Triangulation4.fromIsoSig(
                (args.out / name / "triangulation.isosig").read_text().strip())
            t.simplify()
            same = t.isIsomorphicTo(regina.Example4.sphere()) is not None
            report["test_only_simplification"] = {
                "pentachora": t.size(), "isosig": t.isoSig(),
                "matches_example_sphere": same,
                "certificate": False}
            group = t.fundamentalGroup()
            group.simplify()
            report["test_only_simplification"]["group_generators_after_simplify"] = group.countGenerators()
            require(group.countGenerators() == 0,
                    "The cancellation control did not yield a trivial group presentation")
            # This is a construction control with an explicit geometric
            # cancellation, not a requirement on a heuristic simplifier. A
            # plateau at a larger triangulation is not a failed Kirby test.
        reports.append({"name": name, "pass": True, "report": report})
    receipt = {"status": "PASS", "cases": reports,
               "note": "Construction smoke tests only; these do not admit the Motoo diagrams."}
    (args.out / "CONTROLS_PASS.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("KATIE CONTROLS PASS")


if __name__ == "__main__":
    main()

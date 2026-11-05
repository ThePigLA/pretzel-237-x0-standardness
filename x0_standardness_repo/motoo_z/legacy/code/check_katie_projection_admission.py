#!/usr/bin/env python3
"""Regression: two unlink projections, only the connected one is supported."""
import argparse
import hashlib
import json
from pathlib import Path

import regina

from run_katie import inspect_input


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), "choose a fresh output directory")
    args.out.mkdir(parents=True)
    fixtures = [
        {"name": "disconnected_two_zero_unlinks", "pd_code": [[0, 1, 1, 0], [2, 3, 3, 2]],
         "decorations": [0, 0], "expected_supported": False},
        {"name": "connected_two_zero_unlinks", "pd_code": [[2, 0, 3, 1], [3, 0, 2, 1]],
         "decorations": [0, 0], "expected_supported": True},
    ]
    cases = []
    for fixture in fixtures:
        inspection = inspect_input(fixture)
        supported = not inspection["necessary_shape_failures"]
        require(supported == fixture["expected_supported"], "wrong projection admission")
        # Both diagrams are the SAME underlying two-component unlink:
        # guarded R1/R2 removals reduce each to two crossing-free circles.
        link = regina.Link.fromPD([[a + 1 for a in row] for row in fixture["pd_code"]])
        operations = []
        while link.size():
            changed = False
            for operation in ("r1", "r2"):
                for i in range(link.size()):
                    before = link.tightEncoding()
                    if getattr(link, operation)(link.crossing(i)):
                        operations.append({"operation": operation, "crossing": i,
                                           "before": before, "after": link.tightEncoding()})
                        changed = True
                        break
                if changed:
                    break
            require(changed, "could not verify unlink fixture")
        require(link.countComponents() == 2, "unlink fixture lost a component")
        cases.append({"fixture": fixture, "inspection": inspection,
                      "guarded_unlink_moves": operations, "final_components": 2})
    result = {"status": "PASS_SUPPORTED_PROJECTION_ADMISSION",
              "runner_sha256": hashlib.sha256(Path(__file__).with_name("run_katie.py").read_bytes()).hexdigest(),
              "cases": cases,
              "interpretation": "Disconnected projection is a software-route restriction, not an invalid-link claim."}
    (args.out / "PROJECTION_ADMISSION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "cases": len(cases), "output": str(args.out)}))


if __name__ == "__main__":
    main()

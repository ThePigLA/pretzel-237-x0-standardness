#!/usr/bin/env python3
"""Compare the archived original Figure-2 PD with the marked Katie input.

This is a positive SnapPy marked-exterior comparison. It does not certify
the transcription of a PDF, the paper's Kirby-slide argument, or a 4D move path.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib.metadata
import json
from pathlib import Path
import sys
import time

from common import ROOT, REFERENCE, require, sha256, write_json

ORIGINAL_SHA = "45808e31280cfa002b8a7849683647f291a647909df5967df01ff0e93ef21767"
KATIE_SHA = "5aeaea4ffa70b3d4084fbae11cd144c256ae3af361cbebb10da07166cad3427f"
EXPECTED_IMAGES = [1, 5, 6, 4, 7, 3, 0, 2]


def load_marked(path, expected_hash, snappy):
    require(sha256(path) == expected_hash, "Unexpected marked PD input: " + str(path))
    data = json.loads(path.read_text())
    pd, groups = data["pd_code"], data["component_labels"]
    require(len(groups) == len(data["order"]) == len(data["decorations"]) == 8,
            "Wrong component metadata lengths")
    owner = {}
    for component, group in enumerate(groups):
        require(group and all(type(label) is int for label in group), "Malformed component labels")
        for label in group:
            require(label not in owner, "Repeated component label")
            owner[label] = component
    require([min(group) for group in groups] == sorted(min(group) for group in groups),
            "Declared component order differs from the PD constructor's minimum-label order")
    occurrences = Counter()
    for row in pd:
        require(isinstance(row, list) and len(row) == 4 and
                all(type(label) is int and label in owner for label in row), "Malformed PD row")
        require(owner[row[0]] == owner[row[2]] and owner[row[1]] == owner[row[3]],
                "Opposite half-edges have inconsistent component owners")
        occurrences.update(row)
    require(set(occurrences) == set(owner) and set(occurrences.values()) == {2},
            "PD labels do not each occur twice")
    link = snappy.Link(pd)
    require(len(link.link_components) == 8 and len(link.crossings) == len(pd),
            "Native PD input has wrong component or crossing counts")
    # Spherogram preserves input crossing labels. Its orientation operation
    # rotates by 180 degrees, which preserves the under/over component pair.
    require({crossing.label for crossing in link.crossings} == set(range(len(pd))),
            "Native crossing correspondence differs")
    for crossing in link.crossings:
        row = pd[crossing.label]
        require(list(crossing.strand_components) == [owner[label] for label in row],
                "Native component indices do not match the declared original PD owners")
    exterior = link.exterior()
    require(exterior.num_cusps() == 8, "Wrong exterior cusp count")
    return data, exterior


def audit():
    require(importlib.metadata.version("snappy") == "3.3.2", "Requires snappy==3.3.2")
    import snappy
    first = ROOT / "diagram/original_figure2_pd.json"
    second = REFERENCE / "diagram/fig2_published_katie.json"
    a, ma = load_marked(first, ORIGINAL_SHA, snappy)
    b, mb = load_marked(second, KATIE_SHA, snappy)
    require(a["order"] == ["C" + str(i) for i in range(8)], "Wrong original component names")
    require(b["order"] == ["yellow", "red", "cyan", "purple", "blue", "green", "salmon", "navy"],
            "Wrong Katie component names")
    require(all(a["decorations"][i] == b["decorations"][j]
                for i, j in enumerate(EXPECTED_IMAGES)), "The prescribed map changes a framing/dotted role")
    sa = ma.isometry_signature(of_link=True)
    sb = mb.isometry_signature(of_link=True)
    require(sa == sb, "Marked exterior signatures do not agree")
    accepted = []
    for isometry in ma.is_isometric_to(mb, return_isometries=True):
        images = list(isometry.cusp_images())
        matrices = [[[int(matrix[i, j]) for j in range(2)] for i in range(2)]
                    for matrix in isometry.cusp_maps()]
        if images != EXPECTED_IMAGES:
            continue
        if all(matrix in ([[1, 0], [0, 1]], [[-1, 0], [0, -1]]) for matrix in matrices):
            accepted.append({"cusp_images": images, "peripheral_matrices": matrices})
    require(accepted, "No positive isometry with the required component map and plus/minus identity peripheral maps")
    require(sha256(first) == ORIGINAL_SHA and sha256(second) == KATIE_SHA, "Input changed during comparison")
    return {"schema": "motoo-z-marked-diagram-comparison-v1",
            "status": "PASS_POSITIVE_MARKED_EXTERIOR_COMPARISON",
            "original_pd_sha256": ORIGINAL_SHA, "katie_pd_sha256": KATIE_SHA,
            "original_crossings": len(a["pd_code"]), "katie_crossings": len(b["pd_code"]),
            "original_order": a["order"], "katie_order": b["order"],
            "actual_native_component_ownership_checked": True,
            "accepted_isometries": accepted, "decorated_signature": sa,
            "framings_and_dotted_roles_preserved": True,
            "snappy_distribution": importlib.metadata.version("snappy"),
            "spherogram_distribution": importlib.metadata.version("spherogram"),
            "pdf_transcription_reverified": False, "kirby_slide_argument_reproved": False,
            "four_dimensional_standardness_proved": False,
            "trust": "SnapPy/Spherogram PD exteriors, decorated signatures and positive marked isometry computation; not an independent interval hyperbolicity proof"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.expanduser().resolve()
    require(not out.exists(), "Output file already exists")
    require(not any(out.is_relative_to(ROOT / folder) for folder in
                    ("reference", "legacy", "certificates", "diagram", "provenance", "verification", "environment")),
            "Use an output location outside the released data")
    started = time.monotonic()
    code_before = sha256(Path(__file__))
    result = audit()
    require(sha256(Path(__file__)) == code_before, "Audit code changed during comparison")
    result.update(seconds=time.monotonic() - started, checker_sha256=code_before)
    write_json(out, result)
    print(result["status"])
    print("Receipt: " + str(out))


if __name__ == "__main__":
    main()

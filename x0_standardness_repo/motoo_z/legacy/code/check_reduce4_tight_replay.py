#!/usr/bin/env python3
"""Small exact-state regression for reduce4_checked; no large search required."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import regina

import reduce4_checked as reducer


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rejected(label, function):
    try:
        function()
    except (ValueError, regina.InvalidArgument) as error:
        return {"case": label, "rejected": True, "error": str(error)}
    raise ValueError(f"Negative fixture was accepted: {label}")


def check(output):
    reducer.pinned()
    output = Path(output)
    reducer.require(not output.exists(), "Choose a new output directory")
    output.mkdir(parents=True)
    sphere = reducer.standard_sphere()
    reducer.perform(sphere, ["P", 4, 0])
    source = output / "subdivided_double_simplex.isosig"
    source.write_text(sphere.isoSig() + "\n")
    # Use exactly the labelling that the public source reader imports.
    start_tri = reducer.read4(source)
    state_format = "regina-tight4-v1"
    start = reducer.encode(start_tri, state_format)
    reducer.require(reducer.decode(start, state_format) == start_tri,
                    "Tight encoding failed to preserve exact labelling")
    row = None
    for operation in reducer.legal(start_tri, "P", 0):
        after_tri = reducer.decode(start, state_format)
        reducer.perform(after_tri, operation)
        if reducer.is_standard(after_tri):
            row = {"before": start, "operation": operation,
                   "after": reducer.encode(after_tri, state_format)}
            break
    reducer.require(row is not None, "Fixture has no inverse Pachner move")
    certificate = {"schema": "motoo.checked-4d-path.v1", "state_format": state_format,
                   "source_file": source.name, "source_sha256": sha(source),
                   "start": start, "end": row["after"], "moves": [row],
                   "standard_endpoint": True}
    positive = reducer.replay(certificate, source)
    reducer.require(positive["source_binding_checked"] and positive["standard_endpoint"],
                    "Bound tight-format inverse path failed")
    reducer.save_json(output / "tight_certificate.json", certificate)

    checks = []
    altered_bytes = output / "same_triangulation_altered_bytes.isosig"
    altered_bytes.write_text(source.read_text() + "\n")
    reducer.require(reducer.read4(altered_bytes) == start_tri,
                    "Byte-mutation fixture changed the triangulation")
    checks.append(rejected("source_bytes_changed_only",
                           lambda: reducer.replay(certificate, altered_bytes)))

    relabelling = regina.Isomorphism4.identity(start_tri.size())
    relabelling.setFacetPerm(0, regina.Perm5(0, 1))
    relabelled = relabelling(start_tri)
    relabelled_state = reducer.encode(relabelled, state_format)
    reducer.require(relabelled_state != start and
                    relabelled.isIsomorphicTo(start_tri) is not None,
                    "Relabelling fixture must change labels only")
    bad = copy.deepcopy(certificate)
    bad["start"] = relabelled_state
    checks.append(rejected("isomorphic_but_different_labelled_source",
                           lambda: reducer.replay(bad, source)))

    bad = copy.deepcopy(certificate)
    bad["moves"][0]["after"] = start
    checks.append(rejected("wrong_exact_move_endpoint",
                           lambda: reducer.replay(bad, source)))
    bad = copy.deepcopy(certificate)
    bad["moves"][0]["operation"][2] = 10**9
    checks.append(rejected("invalid_face_index",
                           lambda: reducer.replay(bad, source)))
    bad = copy.deepcopy(certificate)
    bad["standard_endpoint"] = False
    checks.append(rejected("incorrect_standardness_flag",
                           lambda: reducer.replay(bad, source)))
    checks.append(rejected("trailing_encoding_characters",
                           lambda: reducer.decode(start + " ", state_format)))

    partial = copy.deepcopy(certificate)
    partial.update({"end": start, "moves": [], "standard_endpoint": False})
    partial_result = reducer.replay(partial, source)
    reducer.require(partial_result["status"] == "PASS_PARTIAL_PATH" and
                    not partial_result["standard_endpoint"], "Partial path was misclassified")

    # Face indices and guard outcomes must agree whether the preceding move's
    # triangulation is retained or decoded afresh from its exact tight state.
    # Exercise successive growth and reduction, including warmed link caches.
    retained = reducer.decode(start, state_format)
    reducer.admit(retained)
    retained_rows = []
    compared_states = 0
    move_families = [("P", d) for d in range(5)] + [
        ("Z", 0), ("Z", 1), ("Z", 2), ("C", 1), ("F", 1)]
    for step in range(10):
        before = reducer.encode(retained, state_format)
        fresh = reducer.decode(before, state_format)
        choices = {}
        for kind, dimension in move_families:
            current = reducer.legal(retained, kind, dimension)
            rebuilt = reducer.legal(fresh, kind, dimension)
            reducer.require(current == rebuilt,
                            f"Retained/fresh face guard indexing differs: {kind}{dimension}")
            choices[kind, dimension] = current
        compared_states += 1
        if step >= 2 and reducer.is_standard(retained):
            break
        if step < 2:
            operation = choices["P", 4][0]
        else:
            reducing = next((choices[family] for family in
                             (("C", 1), ("P", 0), ("P", 1),
                              ("Z", 2), ("Z", 1), ("Z", 0))
                             if choices[family]), [])
            reducer.require(reducing, "Small growth fixture failed to offer a reduction")
            operation = reducing[0]
        reducer.perform(retained, operation)
        reducer.perform(fresh, operation)
        after = reducer.encode(retained, state_format)
        reducer.require(after == reducer.encode(fresh, state_format),
                        "Retained/fresh move outcomes differ")
        retained_rows.append({"before": before, "operation": operation, "after": after})
    reducer.require(reducer.is_standard(retained), "Small growth/reduction path failed to finish")
    retained_certificate = copy.deepcopy(certificate)
    retained_certificate.update({"moves": retained_rows,
                                 "end": reducer.encode(retained, state_format)})
    retained_result = reducer.replay(retained_certificate, source)
    reducer.require(retained_result["start_pentachora"] == 6 and
                    retained_result["end_pentachora"] == 2,
                    "Retained replay mutated the recorded source-size report")
    reducer.save_json(output / "retained_state_certificate.json", retained_certificate)
    search_result = reducer.search(source, output / "small_search", seconds=2, seed=0)
    reducer.require(search_result["standard_endpoint"], "Small exact-state search failed")
    report = {"schema": "motoo.tight4-replay-regression.v1", "status": "PASS",
              "checker_sha256": sha(__file__), "reducer_sha256": sha(reducer.__file__),
              "state_format": state_format, "positive_replay": positive,
              "partial_path_replay": partial_result, "negative_cases": checks,
              "retained_vs_fresh_states_compared": compared_states,
              "retained_state_replay": retained_result,
              "small_search": search_result,
              "large_control_rerun": False}
    reducer.save_json(output / "TIGHT_REPLAY_REGRESSION.json", report)
    print(json.dumps({"status": "PASS", "negative_cases": len(checks),
                      "source_pentachora": start_tri.size(), "moves": 1}))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    check(args.out)

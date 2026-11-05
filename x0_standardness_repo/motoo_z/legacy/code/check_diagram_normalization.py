#!/usr/bin/env python3
"""Replay a hash-bound chain of marked Reidemeister and strict pickup records."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import regina

from check_pickup_macros import check, verify as verify_pickups
from search_katie_position_growth import apply_checked, snapshot
from run_katie import inspect_input


def canonical_state(link):
    s, _ = snapshot(link)
    components = []
    for c in s["cycles"]:
        if c:
            k = min(range(len(c)), key=lambda i: (abs(c[i]), c[i]))
            c = c[k:]+c[:k]
        components.append(c)
    return {"signs": s["signs"], "components": components}


def digest(state):
    return hashlib.sha256(json.dumps(state, separators=(",", ":")).encode()).hexdigest()


def component_order(data):
    return data.get("pd_component_to_source", data.get("source_component_identity",
                        list(range(len(data["decorations"])))))


def verify_reidemeister(source, certificate):
    pd = source["pd_code"]
    if min(x for c in pd for x in c) == 0:
        pd = [[x+1 for x in c] for c in pd]
    link = regina.Link.fromPD(pd)
    check(canonical_state(link) == certificate["initial_state"], "Initial marked diagram differs from source PD")
    check(certificate["decorations"] == source["decorations"], "Framings/dots changed between records")
    check(component_order(certificate) == component_order(source), "Source component order changed")
    counts = {}
    for i, record in enumerate(certificate["moves"]):
        op = record["op"]
        check(digest(canonical_state(link)) == record["before_sha256"], f"Wrong state before move {i}")
        if op[0] == "R1-":
            removed = {op[1]}
        elif op[0] == "R2-":
            removed = {op[1], link.crossing(op[1]).upper().next().crossing().index()}
        elif op[0] == "R3":
            removed = set()
        else:
            raise ValueError(f"Unsupported move type {op[0]}")
        check(sorted(removed) == record["removed_crossings"], f"Wrong support at move {i}")
        link = apply_checked(link, op)
        check(digest(canonical_state(link)) == record["after_sha256"], f"Wrong state after move {i}")
        counts[op[0]] = counts.get(op[0], 0)+1
    check(canonical_state(link) == certificate["endpoint_state"], "Endpoint marked diagram differs")
    expected_pd = [[x-1 for x in c] for c in link.pdData()]
    check(expected_pd == certificate["pd_code"], "PD export does not encode replayed endpoint")
    if "regina_pd_code" in certificate:
        check(link.pdData() == certificate["regina_pd_code"], "Regina PD export differs")
    return {"status": "PASS", "moves": len(certificate["moves"]), "move_counts": counts,
            "endpoint_crossings": link.size(), "component_order": component_order(certificate),
            "trust": "Guarded Regina 7.4.1 R1/R2/R3; independently checked persistent component ownership"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("certificates", type=Path, nargs="+")
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    check(importlib.metadata.version("regina") == "7.4.1", "Requires Regina Python distribution 7.4.1")
    raw = a.source.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    current = json.loads(raw)
    reports = []
    for path in a.certificates:
        cert_raw = path.read_bytes()
        cert = json.loads(cert_raw)
        check(cert["source_json_sha256"] == hashlib.sha256(raw).hexdigest(), f"Broken source hash at {path}")
        if cert["schema"] == "motoo.marked_reidemeister.v1":
            report = verify_reidemeister(current, cert)
        elif cert["schema"] == "motoo.strict_pickup_macros.v1":
            report = verify_pickups(current, cert)
        else:
            raise ValueError(f"Unsupported certificate schema {cert['schema']}")
        reports.append({"file": str(path), "sha256": hashlib.sha256(cert_raw).hexdigest(), **report})
        raw, current = cert_raw, cert
    shape = inspect_input(current)
    result = {"schema": "motoo.marked_diagram_chain_check.v1", "status": "PASS",
              "original_source_sha256": source_hash, "segments": reports,
              "endpoint_shape": shape,
              "katie_position_admitted": False,
              "orientation_and_planar_discs_require_separate_admission": True,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "PASS", "segments": len(reports),
                      "endpoint_crossings": shape["crossings"],
                      "shape_failures": shape["necessary_shape_failures"],
                      "output": str(a.output)}))


if __name__ == "__main__":
    main()

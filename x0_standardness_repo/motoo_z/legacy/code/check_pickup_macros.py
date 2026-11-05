#!/usr/bin/env python3
"""Independently check strict over/under pickup macros from PD port graphs.

This checker does not invoke Spherogram's pickup, deletion, insertion, or
Reidemeister implementation.  Removing the old and replacement arcs must leave
the same marked planar tangle, with the same two endpoints.  Both removed arcs
have distinct crossings and lie uniformly above (or uniformly below) that
tangle.  Lifting either arc into the corresponding half-space gives an ambient
isotopy between them, relative to the complementary tangle.  The checker also
requires a sphere rotation system on both complete diagrams.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
import json
from pathlib import Path

from orient_katie_diagram import components_of_graph, plane_data


def check(condition, message):
    if not condition:
        raise ValueError(message)


def source_at(data, pd, source_order, port):
    return source_order[data["owner"][pd[port[0]][port[1]]]]


def cut_arc(data, pd, arc, start_out, end_in):
    """Smooth the transverse strand at each removed crossing; leave two ends."""
    removed = {c for c, _ in arc}
    parity = {c: p % 2 for c, p in arc}
    endpoints = {start_out, end_in}
    check(len(endpoints) == 2, "Pickup endpoints coincide")
    complement = {}
    visited_removed_ports = set()
    for c in range(len(pd)):
        if c in removed:
            continue
        for p in range(4):
            start = (c, p)
            if start in endpoints:
                complement[start] = None
                continue
            at = data["twin"][start]
            seen = set()
            while at[0] in removed:
                check(at not in seen, "Complement strand failed to reach an endpoint")
                seen.add(at)
                check(at[1] % 2 != parity[at[0]], "Complement entered the deleted arc")
                opposite = (at[0], (at[1]+2) % 4)
                visited_removed_ports.update((at, opposite))
                at = data["twin"][opposite]
            check(at not in endpoints, "Complement hit a deleted-arc endpoint")
            complement[start] = at
    for p, q in complement.items():
        if q is not None:
            check(complement[q] == p, "Smoothed complement is not paired")
    # A component may become a crossing-free circle when the picked arc is
    # deleted. Retain these circles instead of silently discarding them.
    circles = []
    for c, p in arc:
        first = (c, (p+1) % 4)
        if first in visited_removed_ports:
            continue
        at = first
        seen = set()
        while at not in seen:
            check(at[0] in removed and at[1] % 2 != parity[at[0]],
                  "Unexpected component in cut tangle")
            seen.add(at)
            opposite = (at[0], (at[1]+2) % 4)
            visited_removed_ports.update((at, opposite))
            at = data["twin"][opposite]
        check(at == first, "Closed complement traversal did not close")
        circles.append(first)
    return complement, circles


def same_marked_tangle(old_adj, new_adj, crossing_map, old_data, new_data,
                       old_pd, new_pd, old_sources, new_sources, fixed=None):
    """Solve the allowed rotations (0 or 180 degrees) at surviving crossings."""
    fixed = fixed or {}
    shifts = dict(fixed)
    check(len(set(crossing_map.values())) == len(crossing_map), "Crossing map not injective")
    old_vertices = set(crossing_map)
    check(set(p[0] for p in old_adj) == old_vertices, "Wrong old survivor support")
    check(set(p[0] for p in new_adj) == set(crossing_map.values()),
          "Wrong new survivor support")

    def propagate(seed, initial, existing):
        assigned = dict(existing)
        if seed in assigned and assigned[seed] != initial:
            raise ValueError("Inconsistent port rotation")
        assigned[seed] = initial
        q = deque([seed])
        seen = set()
        while q:
            c = q.popleft()
            if c in seen:
                continue
            seen.add(c)
            s = assigned[c]
            check(s in (0, 2), "A crossing rotation switches over and under")
            for p in range(4):
                oldp = (c, p)
                newp = (crossing_map[c], (p+s) % 4)
                check(source_at(old_data, old_pd, old_sources, oldp) ==
                      source_at(new_data, new_pd, new_sources, newp),
                      "A component identity changed")
                a, b = old_adj[oldp], new_adj[newp]
                check((a is None) == (b is None), "A tangle endpoint changed")
                if a is None:
                    continue
                check(b[0] == crossing_map[a[0]], "Complementary strand connection changed")
                required = (b[1]-a[1]) % 4
                check(required in (0, 2), "Crossing handedness/height changed")
                if a[0] in assigned:
                    check(assigned[a[0]] == required, "Inconsistent complementary rotation")
                else:
                    assigned[a[0]] = required
                q.append(a[0])
        return assigned

    for c in list(shifts):
        shifts = propagate(c, shifts[c], shifts)
    for c in sorted(old_vertices-set(shifts)):
        if c in shifts:
            continue
        success = None
        for s in (0, 2):
            try:
                success = propagate(c, s, shifts)
                break
            except ValueError:
                pass
        check(success is not None, "No marked planar tangle isomorphism exists")
        shifts = success
    return shifts


def check_macro(record):
    before, after = record["before_pd_code"], record["after_pd_code"]
    a, b = plane_data(before), plane_data(after)
    old_sources, new_sources = record["before_component_to_source"], record["after_component_to_source"]
    check(len(old_sources) == len(a["components"]), "Old component identity count mismatch")
    check(len(new_sources) == len(b["components"]), "New component identity count mismatch")
    check(sorted(old_sources) == sorted(new_sources), "Source component bijection changed")
    level = {"under": 0, "over": 1}[record["kind"]]
    arc = [tuple(x) for x in record["old_strand"]]
    check(bool(arc), "Empty original pickup arc")
    check(len(set(c for c, p in arc)) == len(arc), "Original picked arc self-intersects")
    check(all(p % 2 == level for _, p in arc), "Original pickup arc has mixed heights")
    for u, v in zip(arc, arc[1:]):
        check(a["twin"][(u[0], (u[1]+2) % 4)] == v,
              "Original pickup arc is not a consecutive path")
    start = tuple(record["start_endpoint"])
    end = tuple(record["end_endpoint"])
    start_out = (start[0], (start[1]+2) % 4)
    check(a["twin"][arc[0]] == start_out, "Original start endpoint mismatch")
    check(a["twin"][(arc[-1][0], (arc[-1][1]+2) % 4)] == end,
          "Original end endpoint mismatch")
    removed = {c for c, p in arc}
    check(start[0] not in removed and end[0] not in removed, "Removed endpoint crossing")
    check(removed == set(record["removed_crossings"]), "Deletion support mismatch")
    survivor_order = [c for c in range(len(before)) if c not in removed]
    # The producer retains survivors in order and appends the new crossings.
    mapping = {old: new for new, old in enumerate(survivor_order)}
    if "surviving_crossing_map" in record:
        supplied_raw = record["surviving_crossing_map"]
        supplied_items = supplied_raw.items() if isinstance(supplied_raw, dict) else supplied_raw
        supplied = {int(k): v for k, v in supplied_items}
        check(mapping == supplied, "Surviving-crossing map differs from record")
    new_crossings = set(range(len(survivor_order), len(after)))
    check(len(new_crossings) == record["new_crossings"], "New-crossing support mismatch")
    if "new_crossing_indices" in record:
        check(new_crossings == set(record["new_crossing_indices"]),
              "New-crossing indices differ from record")
    picked_source = source_at(a, before, old_sources, arc[0])
    check(picked_source == record["picked_component_source"], "Picked component mismatch")
    old_cut, old_circles = cut_arc(a, before, arc, start_out, end)
    # A rotation system specifies the sphere embedding without extra nesting
    # data only when this complementary projection is connected.  Reject the
    # more general case rather than silently forgetting relative placement.
    projection_adj = defaultdict(set)
    for u, v in old_cut.items():
        if v is not None:
            projection_adj[u[0]].add(v[0])
    check(not old_circles and
          len(components_of_graph({p[0] for p in old_cut}, projection_adj)) == 1,
          "Disconnected complementary projections need an additional embedding certificate")
    success = None
    errors = []
    for start_shift in (0, 2):
        try:
            ns = (mapping[start[0]], (start_out[1]+start_shift) % 4)
            out = ns
            new_arc = []
            while True:
                at = b["twin"][out]
                if at[0] not in new_crossings:
                    ne = at
                    break
                check(at[0] not in {c for c, p in new_arc}, "Replacement arc self-intersects")
                check(at[1] % 2 == level, "Replacement arc has mixed crossing heights")
                check(source_at(b, after, new_sources, at) == picked_source,
                      "Replacement arc lies on a different component")
                new_arc.append(at)
                out = (at[0], (at[1]+2) % 4)
            check({c for c, p in new_arc} == new_crossings, "Replacement omits new crossings")
            check(ne[0] == mapping[end[0]], "Replacement ends at another crossing")
            end_shift = (ne[1]-end[1]) % 4
            check(end_shift in (0, 2), "Replacement changes the endpoint strand")
            check(start[0] != end[0] or start_shift == end_shift,
                  "Inconsistent rotations at a shared endpoint crossing")
            new_cut, new_circles = cut_arc(b, after, new_arc, ns, ne)
            old_circle_labels = Counter(source_at(a, before, old_sources, p) for p in old_circles)
            new_circle_labels = Counter(source_at(b, after, new_sources, p) for p in new_circles)
            check(old_circle_labels == new_circle_labels, "A crossing-free tangle component changed")
            rotations = same_marked_tangle(old_cut, new_cut, mapping, a, b, before, after,
                                          old_sources, new_sources,
                                          {start[0]: start_shift, end[0]: end_shift})
            if "surviving_crossing_port_shifts" in record:
                raw_shifts = record["surviving_crossing_port_shifts"]
                supplied_items = raw_shifts.items() if isinstance(raw_shifts, dict) else raw_shifts
                supplied_shifts = {int(k): v for k, v in supplied_items}
                check(rotations == supplied_shifts, "Surviving port rotations differ from record")
            success = {"new_arc": [list(x) for x in new_arc],
                       "surviving_crossing_rotations": rotations,
                       "crossing_free_complement_components": dict(old_circle_labels)}
            break
        except ValueError as error:
            errors.append(str(error))
    check(success is not None, "Pickup tangle identity failed: "+"; ".join(errors))
    return {"status": "PASS", "kind": record["kind"],
            "crossings_before": len(before), "crossings_after": len(after),
            "removed_crossings": len(removed), "new_crossings": len(new_crossings),
            "picked_component_source": picked_source, **success}


def same_full_diagram(old_pd, new_pd, old_sources, new_sources):
    check(len(old_pd) == len(new_pd), "Initial/final diagram crossing count changed")
    a, b = plane_data(old_pd), plane_data(new_pd)
    mapping = {i: i for i in range(len(old_pd))}
    return same_marked_tangle(a["twin"], b["twin"], mapping, a, b,
                             old_pd, new_pd, old_sources, new_sources)


def verify(source, bundle):
    expected_pd = source["pd_code"]
    expected_order = source.get("pd_component_to_source", source.get("source_component_identity",
                                list(range(len(source["decorations"])))))
    records = []
    for i, macro in enumerate(bundle["macros"]):
        same_full_diagram(expected_pd, macro["before_pd_code"], expected_order,
                          macro["before_component_to_source"])
        checked = check_macro(macro)
        records.append({"index": i, **checked})
        expected_pd, expected_order = macro["after_pd_code"], macro["after_component_to_source"]
    same_full_diagram(expected_pd, bundle["pd_code"], expected_order, bundle["pd_component_to_source"])
    source_decorations = dict(zip(source.get("pd_component_to_source", source.get("source_component_identity",
                                           range(len(source["decorations"])))),
                                  source["decorations"]))
    check([source_decorations[i] for i in bundle["pd_component_to_source"]] == bundle["decorations"],
          "An integer framing or dotted decoration changed")
    return {"schema": "motoo.strict_pickup_independent_check.v1", "status": "PASS",
            "macro_count": len(records), "macro_checks": records,
            "trust": "Pure Python PD rotation and marked-tangle checks; no Spherogram move operations",
            "source_diagram_transcription_audit_inherited": True,
            "katie_standard_position_not_implied": True}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("bundle", type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    source_raw, bundle_raw = a.source.read_bytes(), a.bundle.read_bytes()
    bundle = json.loads(bundle_raw)
    check(hashlib.sha256(source_raw).hexdigest() == bundle["source_json_sha256"], "Source hash mismatch")
    report = verify(json.loads(source_raw), bundle)
    report.update({"source_sha256": hashlib.sha256(source_raw).hexdigest(),
                   "bundle_sha256": hashlib.sha256(bundle_raw).hexdigest(),
                   "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "macros": report["macro_count"],
                      "output": str(a.output)}))


if __name__ == "__main__":
    main()

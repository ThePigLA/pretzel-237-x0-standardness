#!/usr/bin/env python3
"""Check the planar discs of dotted components and orient them for Katie.

This is an admission check, not a diagram simplifier.  It requires a connected
classical projection, simple disjoint dotted loops, and one overpass block and
one underpass block per dot.  It constructs faces directly from the PD rotation
system, selects a common exterior outside disjoint dotted discs, and changes
only the auxiliary orientations of dotted components.  Integer framings and
unoriented embedded components are unchanged.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
import hashlib
import json
from pathlib import Path

def components_of_graph(vertices, adjacency):
    remaining = set(vertices)
    result = []
    while remaining:
        initial = min(remaining)
        found = {initial}
        q = deque([initial])
        remaining.remove(initial)
        while q:
            for v in adjacency[q.popleft()]:
                if v in remaining:
                    remaining.remove(v)
                    found.add(v)
                    q.append(v)
        result.append(found)
    return result


def plane_data(pd):
    """Pure-Python oriented rotation-system reconstruction from CCW PD data."""
    if not pd or any(len(c) != 4 for c in pd):
        raise ValueError("Nonempty four-entry PD tuples are required")
    ends = defaultdict(list)
    for c, ports in enumerate(pd):
        for p, label in enumerate(ports):
            ends[label].append((c, p))
    if any(len(v) != 2 for v in ends.values()):
        raise ValueError("Every PD arc must have exactly two ends")
    twin = {}
    for a, b in ends.values():
        twin[a] = b
        twin[b] = a
    # A directed half-edge d travels to twin(d), then turns clockwise in
    # the local CCW rotation.  The resulting face remains to its left.
    face_next = {d: (e[0], (e[1]-1) % 4) for d, e in twin.items()}
    faces = []
    face_of = {}
    for start in sorted(twin):
        if start in face_of:
            continue
        face = []
        at = start
        while at not in face_of:
            face_of[at] = len(faces)
            face.append(at)
            at = face_next[at]
        if at != start:
            raise ValueError("Invalid face permutation")
        faces.append(face)
    vertex_adj = defaultdict(set)
    arc_adj = defaultdict(set)
    for a, b in ends.values():
        vertex_adj[a[0]].add(b[0])
        vertex_adj[b[0]].add(a[0])
    if len(components_of_graph(range(len(pd)), vertex_adj)) != 1:
        raise ValueError("Disconnected projections require additional placement data")
    if len(pd)-len(ends)+len(faces) != 2:
        raise ValueError("PD rotation system is not a sphere embedding")
    for c in pd:
        for a, b in ((c[0], c[2]), (c[1], c[3])):
            arc_adj[a].add(b)
            arc_adj[b].add(a)
    components = sorted(components_of_graph(ends, arc_adj), key=min)
    owner = {a: i for i, comp in enumerate(components) for a in comp}
    walks = []
    for comp in components:
        # Slot 0 is incoming under, so slot 2 is outgoing under.  This
        # recovers the encoded orientation whenever there is an underpass.
        first = next(((c, 2) for c in range(len(pd)) if pd[c][0] in comp), None)
        if first is None:
            first = min(d for a in comp for d in ends[a])
        walk = []
        at = first
        seen = set()
        while at not in seen:
            seen.add(at)
            walk.append(at)
            target = twin[at]
            at = (target[0], (target[1]+2) % 4)
        if at != first or len(walk) != len(comp):
            raise ValueError("Component traversal failed")
        if any(p == 0 for _, p in walk):
            raise ValueError("PD undercrossing orientations are inconsistent")
        walks.append(walk)
    return {"ends": ends, "twin": twin, "faces": faces,
            "face_of": face_of, "owner": owner, "walks": walks,
            "components": components}


def admission(pd, dotted, exterior_dart=None):
    data = plane_data(pd)
    owner = data["owner"]
    dotted = set(dotted)
    if not dotted:
        return {"status": "NO_DOTTED_COMPONENTS", "reverse_components": []}
    if not dotted <= set(range(len(data["components"]))):
        raise ValueError("Dotted component index is outside the PD component list")
    for c in pd:
        a, b = owner[c[0]], owner[c[1]]
        if a in dotted and b in dotted:
            raise ValueError("Dotted projections have a self-crossing or cross each other")
    nfaces = len(data["faces"])
    records = []
    possible_exteriors = set(range(nfaces))
    for dot in sorted(dotted):
        walk = data["walks"][dot]
        levels = [p % 2 for _, p in walk]
        transitions = sum(a != b for a, b in zip(levels, levels[1:]+levels[:1]))
        if transitions != 2:
            raise ValueError(f"Dotted component {dot} has {transitions} O/U transitions")
        adj = defaultdict(set)
        for arc, ends in data["ends"].items():
            if owner[arc] == dot:
                continue
            a, b = [data["face_of"][d] for d in ends]
            adj[a].add(b)
            adj[b].add(a)
        regions = components_of_graph(range(nfaces), adj)
        if len(regions) != 2:
            raise ValueError("A simple dotted loop did not separate the projection sphere")
        region_of = {f: i for i, region in enumerate(regions) for f in region}
        left = {region_of[data["face_of"][d]] for d in walk}
        right = {region_of[data["face_of"][data["twin"][d]]] for d in walk}
        if len(left) != 1 or len(right) != 1 or left == right:
            raise ValueError("Dotted-loop left and right sides are inconsistent")
        occupied = set()
        for arc, ends in data["ends"].items():
            if owner[arc] in dotted-{dot}:
                occupied.update(region_of[data["face_of"][d]] for d in ends)
        # The disc selected for this dot must contain no other dotted loop.
        allowed_interiors = set(range(2))-occupied
        if not allowed_interiors:
            raise ValueError("Nested dotted circles do not bound disjoint planar discs")
        allowed_exteriors = set().union(*(regions[1-i] for i in allowed_interiors))
        possible_exteriors &= allowed_exteriors
        records.append({"component": dot, "transitions": transitions,
                        "regions": regions, "region_of": region_of,
                        "left_region": next(iter(left))})
    if not possible_exteriors:
        raise ValueError("No common exterior outside disjoint dotted discs")
    if exterior_dart is None:
        # Choosing a different exterior on S^2 does not change the link in S^3.
        exterior = min(possible_exteriors,
                       key=lambda f: (-len(data["faces"][f]), data["faces"][f]))
        exterior_dart = data["faces"][exterior][0]
    else:
        exterior_dart = tuple(exterior_dart)
        exterior = data["face_of"][exterior_dart]
        if exterior not in possible_exteriors:
            raise ValueError("The supplied exterior lies inside a selected dotted disc")
    reverse = []
    serial_records = []
    interiors = []
    for rec in records:
        outside = rec["region_of"][exterior]
        interior = rec["regions"][1-outside]
        if any(interior & other for other in interiors):
            raise ValueError("Selected dotted discs overlap")
        interiors.append(interior)
        ccw = rec["left_region"] != outside
        if not ccw:
            reverse.append(rec["component"])
        serial_records.append({"component": rec["component"],
                               "cyclic_over_under_transitions": rec["transitions"],
                               "interior_face_ids": sorted(interior),
                               "counterclockwise": ccw})
    return {"status": "SHAPE_AND_PLANAR_DISC_CHECKS_PASS",
            "faces": [[list(d) for d in f] for f in data["faces"]],
            "exterior_face": exterior, "exterior_dart": list(exterior_dart),
            "dotted_components": serial_records,
            "reverse_components": reverse,
            "all_dotted_counterclockwise": not reverse,
            "source_isotopy_checked_by_this_tool": False}


def orient(data):
    import regina

    pd = data["pd_code"]
    if min(x for c in pd for x in c) == 0:
        pd = [[x+1 for x in c] for c in pd]
    dotted = [i for i, d in enumerate(data["decorations"]) if d in ("x", ".")]
    before = admission(pd, dotted)
    link = regina.Link.fromPD(pd)
    if link.countComponents() != len(data["decorations"]):
        raise ValueError("Decoration/component count mismatch")
    original_owners = plane_data(pd)["owner"]
    reverse = before["reverse_components"]
    for i in reverse:
        link.reverse(link.component(i))
    # pdData orders crossings by lower-strand encounters in component order.
    crossing_order = []
    for i in range(link.countComponents()):
        start = link.component(i)
        at = start
        while True:
            if at.strand() == 0:
                crossing_order.append(at.crossing().index())
            at = at.next()
            if at == start:
                break
    crossing_map = {old: new for new, old in enumerate(crossing_order)}
    newpd = [list(c) for c in link.pdData()]
    if dotted:
        oldc, oldp = before["exterior_dart"]
        rotation = 2 if original_owners[pd[oldc][0]] in reverse else 0
        new_exterior = (crossing_map[oldc], (oldp+rotation) % 4)
        after = admission(newpd, dotted, new_exterior)
        if not after["all_dotted_counterclockwise"]:
            raise ValueError("Reoriented PD did not retain the expected CCW discs")
    else:
        after = admission(newpd, dotted)
    result = dict(data)
    result["pd_code"] = newpd
    result["geometry_admission"] = {
        "schema": "motoo.katie.planar_position.v1",
        "operation": "auxiliary dotted-component orientation reversals only",
        "reversed_components": reverse, "before": before, "after": after,
        "dotted_position_admitted": True,
        "source_isotopy_must_be_verified_separately": True,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    result = orient(json.loads(raw))
    result["geometry_admission"]["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"output": str(args.output), "dotted_position_admitted": True,
                      "reversed_components": result["geometry_admission"]["reversed_components"]}))


if __name__ == "__main__":
    main()

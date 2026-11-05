#!/usr/bin/env python3
"""Search for Katie's diagram position using checked full-link Reidemeister moves.

This exploratory producer uses ordinary (classical) R2 births, never virtual
R2 moves. Every accepted move is logged and component identities are checked
using surviving crossings. A zero-defect result still needs the dotted-circle
counterclockwise and original drawing admission checks before using Katie.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import time

import regina


def snapshot(link):
    cycles = []
    owners = {}
    for i in range(link.countComponents()):
        first = link.component(i)
        if not first:
            cycles.append([])
            continue
        at = first
        cycle = []
        while True:
            owners[at.id()] = i
            cycle.append((1 if at.strand() else -1) * (at.crossing().index() + 1))
            at = at.next()
            if at == first:
                break
        k = min(range(len(cycle)), key=cycle.__getitem__)
        cycles.append(cycle[k:] + cycle[:k])
    return {"signs": [c.sign() for c in link.crossings()], "cycles": cycles}, owners


def digest(state):
    return hashlib.sha256(json.dumps(state, separators=(",", ":")).encode()).hexdigest()


def defects(link, dotted):
    state, owners = snapshot(link)
    self_cross = 0
    inter_cross = 0
    for i in range(link.size()):
        a, b = owners[2*i], owners[2*i+1]
        if a in dotted and b in dotted:
            if a == b:
                self_cross += 1
            else:
                inter_cross += 1
    excess = 0
    transitions = {}
    for i in sorted(dotted):
        seq = [x > 0 for x in state["cycles"][i]]
        t = sum(a != b for a, b in zip(seq, seq[1:] + seq[:1])) if seq else 0
        transitions[i] = t
        excess += abs(t - 2)
    return (self_cross, inter_cross, excess, link.size()), state, owners, transitions


def value(d):
    return 10*d[0] + 7*d[1] + 2*d[2] + 0.025*d[3]


def faces(link):
    """Return face boundaries as (arc StrandRef id, side), side0=left."""
    next_face = {}
    for c in link.crossings():
        low = 2*c.index()
        up = low+1
        incoming_low = (link.strand(low).prev().id(), 1)
        incoming_up = (link.strand(up).prev().id(), 1)
        outgoing_low = (low, 0)
        outgoing_up = (up, 0)
        rotation = ([incoming_low, outgoing_up, outgoing_low, incoming_up]
                    if c.sign() > 0 else
                    [incoming_low, incoming_up, outgoing_low, outgoing_up])
        for j, dart in enumerate(rotation):
            next_face[(dart[0], 1-dart[1])] = rotation[(j-1) % 4]
    seen = set()
    result = []
    for dart in next_face:
        if dart in seen:
            continue
        current = dart
        face = []
        while current not in seen:
            seen.add(current)
            face.append(current)
            current = next_face[current]
        if current != dart:
            raise RuntimeError("Face traversal is not a permutation")
        result.append(face)
    return result


def candidates(link, rng, grow, dotted, owners, max_births=24):
    result = []
    for c in link.crossings():
        i = c.index()
        if link.hasR1(c):
            result.append(["R1-", i])
        if link.hasR2(c):
            result.append(["R2-", i])
        for side in (0, 1):
            if link.hasR3(c, side):
                result.append(["R3", i, side])
    if grow:
        fs = faces(link)
        possible = []
        for f in fs:
            for j, (a, sa) in enumerate(f):
                for b, sb in f[j+1:]:
                    if a == b:
                        continue
                    if owners[a] not in dotted and owners[b] not in dotted and rng.random() < .85:
                        continue
                    possible.append((a, sa, b, sb))
        rng.shuffle(possible)
        for a, sa, b, sb in possible:
            # The combinatorial face routine is an optimization only. Regina
            # always checks the classical planarity precondition independently.
            if rng.randrange(2):
                a, sa, b, sb = b, sb, a, sa
            if link.hasR2(link.strand(a), sa, link.strand(b), sb):
                result.append(["R2+", a, sa, b, sb])
                max_births -= 1
                if not max_births:
                    break
    return result


def apply_checked(link, move, expected_owners=None):
    """Mutate a copy; compare component memberships of every surviving strand."""
    t = regina.Link(link)
    _, before_owners = snapshot(t)
    removed = set()
    if move[0] == "R1-":
        removed.add(move[1])
    elif move[0] == "R2-":
        removed.update((move[1], t.strand(2*move[1]+1).next().crossing().index()))
    survivors = [(c, c.index()) for c in t.crossings() if c.index() not in removed]
    if move[0] == "R1-":
        ok = t.r1(t.crossing(move[1]))
    elif move[0] == "R2-":
        ok = t.r2(t.crossing(move[1]))
    elif move[0] == "R3":
        ok = t.r3(t.crossing(move[1]), move[2])
    elif move[0] == "R2+":
        ok = t.r2(t.strand(move[1]), move[2], t.strand(move[3]), move[4])
    else:
        raise ValueError(move)
    if not ok:
        raise RuntimeError("Logged Reidemeister move is not legal")
    if t.countComponents() != link.countComponents():
        raise RuntimeError("Reidemeister move changed component count")
    state, after_owners = snapshot(t)
    for crossing, old in survivors:
        new = crossing.index()
        for strand in (0, 1):
            if after_owners[2*new+strand] != before_owners[2*old+strand]:
                raise RuntimeError("Reidemeister move changed a component identity")
    return t


def write_json(path, data):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    temp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--seconds", type=float, default=180)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--ceiling", type=int, default=120)
    p.add_argument("--temperature", type=float, default=1.6)
    args = p.parse_args()
    data = json.loads(args.input.read_text())
    pd = data["pd_code"]
    if min(x for c in pd for x in c) == 0:
        pd = [[x+1 for x in c] for c in pd]
    link = regina.Link.fromPD(pd)
    decorations = data["decorations"]
    if link.countComponents() != len(decorations):
        raise RuntimeError("Decoration/component mismatch")
    dotted = {i for i, d in enumerate(decorations) if d in ("x", ".")}
    args.out.mkdir(parents=True, exist_ok=False)
    initial, _ = snapshot(link)
    write_json(args.out/"initial.json", {"input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
                                        "state": initial, "decorations": decorations,
                                        "pd_component_to_source": data.get("pd_component_to_source"),
                                        "native_moves": "Regina 7.4.1 classical R1/R2/R3"})
    rng = random.Random(args.seed)
    current, state, owners, transitions = defects(link, dotted)
    best = current
    step = 0
    epoch = time.monotonic()
    seen = {digest(state): 0}
    last_birth = -1000
    last_improvement = 0
    print("START", current, transitions, flush=True)
    with (args.out/"moves.jsonl").open("w", buffering=1) as log:
        while time.monotonic()-epoch < args.seconds:
            if current[:3] == (0, 0, 0):
                break
            grow = link.size()+2 <= args.ceiling and (step-last_improvement > 40 or rng.random() < .08)
            moves = candidates(link, rng, grow, dotted, owners)
            rng.shuffle(moves)
            options = []
            for move in moves:
                t = apply_checked(link, move)
                d, s, own, trans = defects(t, dotted)
                h = digest(s)
                age = step-seen.get(h, -10000)
                if age < 12:
                    continue
                # An R2 birth needs several R3 opportunities before its inverse
                # can be selected; the visited-state taboo prevents instant undo.
                energy = value(d)
                if move[0] == "R2+":
                    energy -= 1.3 if step-last_improvement > 80 else .2
                options.append((energy, rng.random(), move, t, d, s, own, trans, h))
            if not options:
                # Forget local taboo only; no topological operation occurs.
                seen = {digest(state): step}
                last_improvement = min(last_improvement, step-100)
                continue
            minimum = min(x[0] for x in options)
            temperature = args.temperature * (1.5 if step-last_improvement > 300 else 1)
            weights = [math.exp(-(x[0]-minimum)/temperature) for x in options]
            chosen = rng.choices(options, weights=weights, k=1)[0]
            _, _, move, t, d, s, own, trans, h = chosen
            previous = digest(state)
            step += 1
            log.write(json.dumps({"step": step, "move": move, "before": previous,
                                  "after": h, "defects": d}, separators=(",", ":"))+"\n")
            link, current, state, owners, transitions = t, d, s, own, trans
            seen[h] = step
            if len(seen) > 12000:
                seen = {k: v for k, v in seen.items() if step-v < 5000}
            if value(current) < value(best):
                best = current
                last_improvement = step
                checkpoint = {"step": step, "defects": current, "transitions": transitions,
                              "state": state, "state_sha256": h,
                              "decorations": decorations,
                              "counterclockwise_admission": False,
                              "status": "SHAPE_GOAL" if current[:3] == (0,0,0) else "SEARCH_CHECKPOINT"}
                write_json(args.out/"best.json", checkpoint)
                print("BEST", step, current, transitions, flush=True)
            if step % 100 == 0:
                write_json(args.out/"progress.json", {"step": step, "current": current,
                                                      "best": best, "seconds": time.monotonic()-epoch})
        write_json(args.out/"final.json", {"step": step, "current": current, "best": best,
                                           "state": state, "decorations": decorations,
                                           "seconds": time.monotonic()-epoch,
                                           "status": "SHAPE_GOAL" if current[:3] == (0,0,0) else "BUDGET_EXHAUSTED"})
    print("DONE", step, current, "BEST", best, flush=True)


if __name__ == "__main__":
    main()

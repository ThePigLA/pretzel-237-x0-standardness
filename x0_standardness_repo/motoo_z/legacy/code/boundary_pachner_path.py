#!/usr/bin/env python3
"""Find and independently replay an explicit 3D Pachner path.

Search is heuristic.  Each retained step starts in its canonical isosig
state, performs exactly one guarded move, and checks the next isosig.
Only 4-1, 3-2, 2-3 and 1-4 are used, so every retained step can be
realised by one 4-dimensional layering simplex.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import time

import regina

from complete_upper_handles import all_faces, get_face, require, save_json


def replay(path):
    before = path["start"]
    moves = path["moves"]
    for i, move in enumerate(moves):
        require(move["before"] == before, f"discontinuous path at {i}")
        tri = regina.Triangulation3.fromIsoSig(before)
        require(tri.pachner(get_face(tri, move["dimension"], move["face"])),
                f"illegal move at {i}")
        require(tri.isoSig() == move["after"], f"wrong move endpoint at {i}")
        before = move["after"]
    require(before == path["end"], "wrong complete endpoint")
    return {"status": "PASS", "moves": len(moves), "end": before}


def search(start, output, seconds=60, seed=0, goal_tetrahedra=6, goals=None):
    output = Path(output)
    require(not output.exists(), "output already exists")
    output.mkdir(parents=True)
    rng = random.Random(seed)
    now = time.monotonic()
    tri = regina.Triangulation3.fromIsoSig(start)
    require(tri.isValid() and tri.isClosed() and tri.isConnected(), "invalid source")
    states, moves, positions = [start], [], {start: 0}
    best_score = (tri.countVertices(), tri.size())
    best_path = {"start": start, "end": start, "moves": []}
    mode = "reduce"
    growth_target = tri.size()
    steps = 0
    max_tetrahedra = tri.size()
    counters = [0, 0, 0, 0]
    last_save = now
    restarts = 0

    def checkpoint(status):
        result = {"schema": "motoo.boundary.pachner-path.v1", "status": status,
                  "seed": seed, "elapsed_seconds": time.monotonic() - now,
                  "attempted_moves": steps, "move_counts_by_face_dimension": counters,
                  "max_tetrahedra": max_tetrahedra, "best_score_vertices_tetrahedra": best_score,
                  "restarts": restarts, **best_path}
        result["replay"] = replay(result)
        save_json(output / "best_path.json", result)
        return result

    while time.monotonic() - now < seconds:
        sig = states[-1]
        target_hit = sig in goals if goals else tri.size() <= goal_tetrahedra
        if target_hit:
            best_path = {"start": start, "end": sig, "moves": list(moves)}
            best_score = (tri.countVertices(), tri.size())
            result = checkpoint("FOUND_AND_REPLAYED")
            save_json(output / "SUCCESS.json", result)
            print(json.dumps({"status": result["status"], "steps": steps,
                              "retained_moves": len(moves), "end": sig}), flush=True)
            return result
        reductions = []
        for dimension in (0, 1):
            reductions.extend((dimension, f.index()) for f in all_faces(tri, dimension)
                              if tri.hasPachner(f))
        growers = [(2, f.index()) for f in tri.triangles() if tri.hasPachner(f)]
        cap = best_score[1] + rng.choice((8, 12, 16))
        if mode == "reduce" and reductions:
            # A vertex deletion is preferred whenever possible.
            vertex_deletions = [p for p in reductions if p[0] == 0]
            choices = vertex_deletions or reductions
        else:
            if mode == "reduce":
                if tri.size() >= cap or tri.countVertices() > best_score[0] + 1:
                    states = [best_path["start"]] + [m["after"] for m in best_path["moves"]]
                    moves = list(best_path["moves"])
                    positions = {s: i for i, s in enumerate(states)}
                    tri = regina.Triangulation3.fromIsoSig(states[-1])
                    mode = "reduce"
                    restarts += 1
                    continue
                mode = "grow"
                growth_target = min(cap, tri.size() + rng.choice((4, 6, 8, 12)))
            if tri.size() >= growth_target:
                mode = "reduce"
                choices = reductions or growers
            else:
                choices = growers
                if tri.size() + 3 <= growth_target and rng.random() < 0.015:
                    choices = [(3, rng.randrange(tri.size()))]
        if not choices:
            result = checkpoint("NO_LEGAL_SELECTED_MOVE")
            return result
        dimension, index = rng.choice(choices)
        require(tri.pachner(get_face(tri, dimension, index)), "selected move failed guard")
        after = tri.isoSig()
        steps += 1
        counters[dimension] += 1
        max_tetrahedra = max(max_tetrahedra, tri.size())
        if after in positions:
            where = positions[after]
            for old in states[where + 1:]:
                del positions[old]
            del states[where + 1:]
            del moves[where:]
        else:
            moves.append({"before": sig, "dimension": dimension, "face": index, "after": after})
            positions[after] = len(states)
            states.append(after)
        tri = regina.Triangulation3.fromIsoSig(after)
        score = (tri.countVertices(), tri.size())
        if score < best_score:
            best_score = score
            best_path = {"start": start, "end": after, "moves": list(moves)}
            print(json.dumps({"best": score, "steps": steps, "path": len(moves),
                              "seconds": time.monotonic() - now}), flush=True)
            mode = "reduce"
        if time.monotonic() - last_save > 20:
            checkpoint("RUNNING")
            last_save = time.monotonic()
    result = checkpoint("TIME_BUDGET_EXHAUSTED")
    print(json.dumps({"status": result["status"], "best": best_score,
                      "steps": steps, "retained_moves": len(best_path["moves"])}), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--goal-size", type=int, default=6)
    parser.add_argument("--goal", action="append")
    args = parser.parse_args()
    search(args.start, args.out, args.seconds, args.seed, args.goal_size, set(args.goal or []))

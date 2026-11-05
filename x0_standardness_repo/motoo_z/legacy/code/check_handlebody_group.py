#!/usr/bin/env python3
"""Compute a dotted-link handlebody group and save elementary Tietze steps.

The input is an audited PD code with a component decoration vector.  The
topological interpretation requires the dotted sublink to be the unlink.
For every framed component, its meridian and longitude are killed in the
full link-exterior group.  The longitude framing is immaterial after its
meridian has been killed.  This script also verifies the dotted unlink by
guarded Reidemeister moves on the extracted dotted subdiagram.

`replay REPORT` requires only the Python standard library, and checks the
saved elementary group eliminations.  It does not re-extract the group.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def reduced(word):
    out = []
    for x in word:
        require(isinstance(x, int) and x != 0, "invalid generator letter")
        if out and out[-1] == -x:
            out.pop()
        else:
            out.append(x)
    return out


def inverse(word):
    return [-x for x in reversed(word)]


def cyclic_reduced(word):
    w = reduced(word)
    i, j = 0, len(w)
    while j - i >= 2 and w[i] == -w[j - 1]:
        i += 1
        j -= 1
    return w[i:j]


def cleaned(relators):
    out = []
    seen = set()
    for word in relators:
        w = cyclic_reduced(word)
        if w and tuple(w) not in seen:
            out.append(w)
            seen.add(tuple(w))
    return out


def digest_state(active, relators):
    b = json.dumps([sorted(active), relators], separators=(",", ":")).encode()
    return hashlib.sha256(b).hexdigest()


def apply_elimination(active, relators, generator, index):
    require(generator in active, "generator is not active")
    require(0 <= index < len(relators), "relator index out of range")
    pivot = relators[index]
    positions = [i for i, x in enumerate(pivot) if abs(x) == generator]
    require(len(positions) == 1, "pivot must contain the generator exactly once")
    j = positions[0]
    tail = pivot[j + 1:] + pivot[:j]
    replacement = reduced(inverse(tail) if pivot[j] > 0 else tail)
    require(all(abs(x) != generator for x in replacement), "self-dependent replacement")
    rep_inverse = inverse(replacement)
    changed = []
    for i, word in enumerate(relators):
        if i == index:
            continue
        expanded = []
        for x in word:
            expanded.extend(replacement if x == generator else rep_inverse if x == -generator else [x])
        changed.append(expanded)
    return active - {generator}, cleaned(changed), replacement


def find_shortening(relators):
    """Find a Dehn replacement using a cyclic conjugate of another relator."""
    best = None
    for si, source in enumerate(relators):
        m = len(source)
        for invert in (False, True):
            oriented = inverse(source) if invert else source
            for rotation in range(m):
                s = oriented[rotation:] + oriented[:rotation]
                for ti, target in enumerate(relators):
                    if ti == si:
                        continue
                    n = len(target)
                    doubled = target + target
                    for amount in range(min(m, n), m // 2, -1):
                        prefix = s[:amount]
                        for start in range(n):
                            if doubled[start:start + amount] != prefix:
                                continue
                            rotated_target = target[start:] + target[:start]
                            result = cyclic_reduced(inverse(s[amount:]) + rotated_target[amount:])
                            gain = n - len(result)
                            if gain <= 0:
                                continue
                            key = (-gain, len(result), si, ti, invert, rotation, amount, start)
                            if best is None or key < best[0]:
                                best = (key, {"source_index": si, "target_index": ti,
                                              "source_inverted": invert, "source_rotation": rotation,
                                              "target_rotation": start, "matched_length": amount,
                                              "source_word": source, "target_word": target,
                                              "new_target": result})
    return None if best is None else best[1]


def apply_shortening(relators, step):
    si, ti = step["source_index"], step["target_index"]
    require(0 <= si < len(relators) and 0 <= ti < len(relators) and si != ti, "shortening indices")
    source, target = relators[si], relators[ti]
    require(source == step["source_word"] and target == step["target_word"], "shortening source mismatch")
    s = inverse(source) if step["source_inverted"] else source[:]
    a, b, k = step["source_rotation"], step["target_rotation"], step["matched_length"]
    require(0 <= a < len(s) and 0 <= b < len(target), "shortening rotation")
    s = s[a:] + s[:a]
    t = target[b:] + target[:b]
    require(1 <= k <= min(len(s), len(t)), "shortening matched length")
    require(s[:k] == t[:k], "shortening prefix mismatch")
    result = cyclic_reduced(inverse(s[k:]) + t[k:])
    require(result == step["new_target"] and len(result) < len(target), "shortening result mismatch")
    changed = [result if i == ti else w for i, w in enumerate(relators)]
    return cleaned(changed)


def eliminate(generators, original_relators):
    active = set(range(1, generators + 1))
    relators = cleaned(original_relators)
    steps = []
    while active:
        choices = []
        for i, w in enumerate(relators):
            counts = {}
            for x in w:
                counts[abs(x)] = counts.get(abs(x), 0) + 1
            for g, count in counts.items():
                if count == 1:
                    choices.append((len(w), g, i))
        if not choices:
            shortening = find_shortening(relators)
            if shortening is None:
                break
            before = digest_state(active, relators)
            relators = apply_shortening(relators, shortening)
            steps.append({"operation": "shorten_relator", **shortening,
                          "before_sha256": before, "after_sha256": digest_state(active, relators)})
            continue
        _, g, i = min(choices)
        before = digest_state(active, relators)
        pivot = list(relators[i])
        active, relators, replacement = apply_elimination(active, relators, g, i)
        steps.append({"operation": "eliminate_generator", "generator": g, "relator_index": i, "pivot": pivot,
                      "replacement": replacement, "before_sha256": before,
                      "after_sha256": digest_state(active, relators)})
    return {"initial_generator_count": generators, "initial_relators": original_relators,
            "steps": steps, "final_active_generators": sorted(active),
            "final_relators": relators, "trivial": not active and not relators}


def replay(proof):
    n = proof["initial_generator_count"]
    require(n >= 0, "invalid generator count")
    original = proof["initial_relators"]
    require(all(1 <= abs(x) <= n for w in original for x in w), "generator range")
    active = set(range(1, n + 1))
    relators = cleaned(original)
    for k, step in enumerate(proof["steps"]):
        require(digest_state(active, relators) == step["before_sha256"], f"pre-state mismatch at step {k}")
        if step.get("operation") == "shorten_relator":
            relators = apply_shortening(relators, step)
            require(digest_state(active, relators) == step["after_sha256"], f"post-state mismatch at step {k}")
            continue
        i = step["relator_index"]
        require(0 <= i < len(relators), "invalid pivot index")
        require(relators[i] == step["pivot"], f"pivot mismatch at step {k}")
        active, relators, replacement = apply_elimination(active, relators, step["generator"], i)
        require(replacement == step["replacement"], f"replacement mismatch at step {k}")
        require(digest_state(active, relators) == step["after_sha256"], f"post-state mismatch at step {k}")
    require(sorted(active) == proof["final_active_generators"], "final generator mismatch")
    require(relators == proof["final_relators"], "final relator mismatch")
    result = not active and not relators
    require(result == proof["trivial"], "incorrect triviality flag")
    return result


def check_dotted_unlink(data):
    import regina
    dotted = [c["index"] for c in data["components"] if c["dotted"]]
    keep = [c for c in data["crossings"] if c["under_component"] in dotted and c["over_component"] in dotted]
    renumber = {c["id"]: i + 1 for i, c in enumerate(keep)}
    signs = [c["sign"] for c in keep]
    cycles = []
    for j in dotted:
        visits = data["components"][j]["crossing_visits"]
        cycles.append([(1 if v["role"] == "over" else -1) * renumber[v["crossing"]]
                       for v in visits if v["crossing"] in renumber])
    link = regina.Link.fromData(signs, cycles)
    require(link.countComponents() == len(dotted), "lost dotted component")
    require(link.virtualGenus() == 0, "dotted sublink is not a planar diagram")
    initial = link.tightEncoding()
    steps = []
    while link.size():
        found = False
        for kind in ("r1", "r2"):
            for i in range(link.size()):
                before = link.tightEncoding()
                if getattr(link, kind)(link.crossing(i)):
                    steps.append({"move": kind, "crossing": i, "before": before, "after": link.tightEncoding()})
                    found = True
                    break
            if found:
                break
        if not found:
            break
    # Freshly replay the recorded, individually guarded Reidemeister moves.
    replayed = regina.Link.fromData(signs, cycles)
    for step in steps:
        require(replayed.tightEncoding() == step["before"], "unlink replay pre-state")
        require(getattr(replayed, step["move"])(replayed.crossing(step["crossing"])), "unlink move rejected")
        require(replayed.tightEncoding() == step["after"], "unlink replay post-state")
    require(replayed.tightEncoding() == link.tightEncoding(), "unlink endpoint mismatch")
    result = link.size() == 0 and link.countComponents() == len(dotted)
    return {"source_components": dotted, "crossing_signs": signs, "component_cycles": cycles,
            "initial_crossings": len(keep), "initial_encoding": initial, "moves": steps,
            "final_crossings": link.size(), "final_components": link.countComponents(),
            "final_encoding": link.tightEncoding(), "proved_unlink": result,
            "trust": "Regina guarded classical Reidemeister moves"}


def compute(path, output):
    import snappy
    import regina
    data = json.loads(path.read_text())
    snappy.set_rand_seed(0)
    decorations = data["decorations"]
    pd = data["pd_code"]
    require(len(decorations) == data["component_count"], "decoration count")
    require(data["pd_component_to_source"] == list(range(data["component_count"])),
            "this initial checker requires explicit identity source mapping")
    dotted = check_dotted_unlink(data)
    link = snappy.Link(pd)
    require(len(link.link_components) == len(decorations), "PD component count mismatch")
    exterior = link.exterior()
    require(exterior.num_cusps() == len(decorations), "cusp count mismatch")
    group = exterior.fundamental_group(simplify_presentation=False)
    n = group.num_generators()
    raw_relators = group.relators(as_int_list=True)
    periph = group.peripheral_curves(as_int_list=True)
    relators = list(raw_relators)
    added = []
    for j, decoration in enumerate(decorations):
        if decoration in ("x", "."):
            continue
        for kind, word in zip(("meridian", "longitude"), periph[j]):
            relators.append(list(word))
            added.append({"component": j, "source_component": data["pd_component_to_source"][j],
                          "curve": kind, "word": list(word)})
    proof = eliminate(n, relators)
    replay(proof)
    result = {
        "schema": "motoo.handlebody-group.v1", "input_filename": path.name,
        "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "versions": {"regina_distribution": importlib.metadata.version("regina"),
                     "regina_engine": regina.versionString(), "snappy": snappy.__version__},
        "pd_code": pd, "decorations": decorations,
        "dotted_sublink": dotted,
        "link_exterior": {"tetrahedra": exterior.num_tetrahedra(), "cusps": exterior.num_cusps(),
                          "isomorphism_signature": exterior.triangulation_isosig(),
                          "snappy_triangulation": exterior._to_string(),
                          "random_seed": 0,
                          "generator_count": n, "relators": raw_relators, "peripheral_curves": periph},
        "added_relations": added, "tietze_proof": proof,
        "interpretation": "4D handlebody pi1 is the link-exterior group modulo each framed component's meridian and longitude, once dotted unlink is established",
        "status": ("TRIVIAL_HANDLEBODY_PI1" if dotted["proved_unlink"] and proof["trivial"]
                   else "NOT_RESOLVED_BY_THIS_CHECK"),
        "qualification": "Conditional on the supplied diagram transcription; does not establish the upper-handle completion or smooth standardness."
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(output), "status": result["status"],
                      "initial_generators": n, "group_steps": len(proof["steps"]),
                      "dotted_unlink": dotted["proved_unlink"],
                      "dotted_reidemeister_moves": len(dotted["moves"])}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("compute")
    c.add_argument("input", type=Path)
    c.add_argument("output", type=Path)
    r = sub.add_parser("replay")
    r.add_argument("report", type=Path)
    sub.add_parser("selftest")
    args = parser.parse_args()
    if args.command == "compute":
        compute(args.input, args.output)
    elif args.command == "replay":
        report = json.loads(args.report.read_text())
        result = replay(report["tietze_proof"])
        print(json.dumps({"tietze_replay_passed": True, "trivial_quotient": result}))
    else:
        fixtures = [(3, [[1, -2], [2], [1, 3, -1]], True),
                    (1, [], False), (1, [[1, 1]], False),
                    (2, [[1, 1], [2, 2, 2], [1, 2, 1, 2]], False)]
        for n, rels, expected in fixtures:
            proof = eliminate(n, rels)
            require(replay(proof) == expected, "fixture failed")
        bad = eliminate(1, [[1]])
        bad["steps"][0]["replacement"] = [1]
        try:
            replay(bad)
        except ValueError:
            pass
        else:
            raise ValueError("corrupt proof accepted")
        print("PASS: trivial, free, torsion and S3 examples; altered certificate rejected")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Bounded 4D search with explicit guarded moves and separate certificate replay.

This is not an exhaustive search. A failed run has no classification meaning.
All move guards and manifold recognition use pinned Regina 7.4.1. The search
never calls Regina's randomized simplify(). Its standard endpoint is an
explicit double of a 4-simplex, built below from five identity face gluings.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import random
import time

import regina


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pinned():
    require(importlib.metadata.version("regina") == "7.4.1", "requires Regina distribution 7.4.1")


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def standard_sphere():
    tri = regina.Triangulation4()
    first, second = tri.newSimplices(2)
    for facet in range(5):
        first.join(facet, second, regina.Perm5())
    return tri


def read4(path):
    path = Path(path)
    if path.suffix == ".rga":
        packet = regina.open(str(path))
        require(isinstance(packet, regina.Triangulation4), "file root must be a 4D triangulation")
        return regina.Triangulation4(packet)
    return regina.Triangulation4.fromIsoSig(path.read_text().strip())


def decode(state, state_format):
    if state_format == "regina-tight4-v1":
        return regina.Triangulation4.tightDecoding(state)
    require(state_format == "regina-isosig4", "unknown state encoding")
    return regina.Triangulation4.fromIsoSig(state)


def encode(tri, state_format):
    return tri.tightEncoding() if state_format == "regina-tight4-v1" else tri.isoSig()


def is_standard(tri):
    return tri.size() == 2 and tri.isIsomorphicTo(standard_sphere()) is not None


def admit(tri):
    require(tri.size() > 0 and tri.isValid() and tri.isConnected() and
            not tri.hasBoundaryFacets() and tri.isOrientable(),
            "source must be a nonempty valid connected finite closed orientable triangulation")
    # Spherical links below establish that no ideal vertices are present.
    # Calling isClosed()/isIdeal() first can repeat the same recognitions.
    links = []
    for vertex in tri.vertices():
        link = vertex.buildLink()
        recognised = link.isSphere()
        require(recognised, f"source vertex {vertex.index()} has a nonspherical link")
        links.append({"vertex": vertex.index(), "tetrahedra": link.size(), "sphere": True})
    require(tri.homology().isTrivial() and tri.homology(2).isTrivial()
            and tri.eulerCharManifold() == 2, "source fails homology-sphere filter")
    return {"status": "PASS_CLOSED_PL_HOMOLOGY_SPHERE",
            "pentachora": tri.size(), "vertices": tri.countVertices(),
            "H1": str(tri.homology()), "H2": str(tri.homology(2)),
            "euler_characteristic": tri.eulerCharManifold(), "vertex_links": links,
            "simply_connectedness_not_inferred_from_homology": True}


def faces(tri, dimension):
    return (tri.vertices, tri.edges, tri.triangles, tri.tetrahedra, tri.pentachora)[dimension]()


def face(tri, dimension, index):
    require(type(dimension) is int and dimension in range(5), "bad face dimension")
    require(type(index) is int and 0 <= index < len(faces(tri, dimension)), "bad face index")
    return (tri.vertex, tri.edge, tri.triangle, tri.tetrahedron, tri.pentachoron)[dimension](index)


def perform(tri, operation):
    require(isinstance(operation, list) and len(operation) == 3, "bad operation record")
    kind, dimension, index = operation
    target = face(tri, dimension, index)
    before_size = tri.size()
    if kind == "P":
        require(tri.pachner(target), "Pachner move rejected")
        require(tri.size() - before_size == 2 * dimension - 4, "wrong Pachner size change")
    elif kind == "Z":
        require(dimension in (0, 1, 2), "bad 2-0 face dimension")
        require(tri.move20(target), "2-0 move rejected")
        require(tri.size() == before_size - 2, "wrong 2-0 size change")
    elif kind == "C":
        require(dimension == 1, "collapse is not on an edge")
        old_vertices = tri.countVertices()
        require(tri.collapseEdge(target), "edge collapse rejected")
        require(tri.countVertices() == old_vertices - 1 and tri.size() < before_size,
                "wrong collapse size change")
    elif kind == "F":
        require(dimension == 1, "4-4 move is not on an edge")
        require(tri.move44(target), "4-4 move rejected")
        require(tri.size() == before_size, "wrong 4-4 size change")
    else:
        raise ValueError("unknown move code")
    # Starting manifold admission plus guarded topology-preserving moves
    # propagates finite/closed status. Do not re-run S3 recognition at every
    # step: 4D isClosed()/isIdeal() can trigger those expensive recognitions.
    require(tri.isValid() and tri.isConnected() and not tri.hasBoundaryFacets()
            and tri.isOrientable(), "move output fails structural checks")


def legal(tri, kind, dimension):
    guard = {"P": tri.hasPachner, "Z": tri.has20,
             "C": tri.hasCollapseEdge, "F": tri.has44}[kind]
    return [[kind, dimension, f.index()] for f in faces(tri, dimension) if guard(f)]


def replay(certificate, source_file=None):
    pinned()
    start = certificate["start"]
    state_format = certificate.get("state_format", "regina-isosig4")
    if source_file is not None:
        source_file = Path(source_file)
        require(hashlib.sha256(source_file.read_bytes()).hexdigest() == certificate["source_sha256"],
                "source file hash mismatch")
        require(encode(read4(source_file), state_format) == start,
                "certificate starts at another triangulation")
    source = decode(start, state_format)
    admission = admit(source)
    source_size = source.size()
    tri = source
    previous = start
    counts = Counter()
    for i, row in enumerate(certificate["moves"]):
        require(row["before"] == previous, f"discontinuous move chain at {i}")
        # Tight states preserve simplex labels. Keep the triangulation whose
        # previous guarded move has already updated its topology caches.
        # Isomorphism signatures instead prescribe a fresh canonical labelling.
        if state_format == "regina-isosig4":
            tri = decode(previous, state_format)
        perform(tri, row["operation"])
        actual = encode(tri, state_format)
        require(actual == row["after"], f"wrong endpoint at move {i}")
        previous = actual
        counts[str(row["operation"][:2])] += 1
    require(previous == certificate["end"], "wrong final endpoint")
    standard = is_standard(decode(previous, state_format))
    require(standard == certificate["standard_endpoint"], "incorrect standardness flag")
    return {"status": "PASS_STANDARDIZATION" if standard else "PASS_PARTIAL_PATH",
            "moves": len(certificate["moves"]), "move_counts": dict(counts),
            "start_pentachora": source_size,
            "end_pentachora": decode(previous, state_format).size(),
            "source_binding_checked": source_file is not None,
            "source_admission": admission, "standard_endpoint": standard,
            "trust": "Regina 7.4.1 guarded moves and 3-sphere recognition; explicit double-simplex target"}


def search(source_file, output, seconds=120, seed=0, growth=24):
    pinned()
    output = Path(output)
    require(not output.exists(), "choose a new output directory")
    output.mkdir(parents=True)
    source_file = Path(source_file)
    tri = read4(source_file)
    print(json.dumps({"stage": "source_admission", "pentachora": tri.size()}), flush=True)
    admission = admit(tri)
    save_json(output / "SOURCE_ADMISSION.json", admission)
    state_format = "regina-tight4-v1"
    start = encode(tri, state_format)
    target = standard_sphere().isoSig()
    rng = random.Random(seed)
    started = time.monotonic()
    states, moves, positions = [start], [], {start: 0}
    best_end, best_moves = start, []
    best_size = tri.size()
    accepted = 0
    max_size = tri.size()
    saved = started
    last_progress = started
    mode = "reduce"
    growth_left = 0
    neutral_left = 0
    interrupted = False

    def checkpoint():
        certificate = {
            "schema": "motoo.checked-4d-path.v1", "source_file": source_file.name,
            "source_sha256": hashlib.sha256(source_file.read_bytes()).hexdigest(),
            "start": start, "end": best_end, "moves": list(best_moves),
            "state_format": state_format,
            "standard_endpoint": is_standard(decode(best_end, state_format)),
            "reference": "Two 4-simplices glued by the identity on each of their five facets",
            "reference_isosig": target, "regina_distribution": "7.4.1",
            "seed": seed, "elapsed_seconds": time.monotonic() - started,
            "accepted_search_moves": accepted, "max_pentachora": max_size,
            "interrupted": interrupted, "exhaustive": False,
            "native_simplify_used": False,
        }
        save_json(output / "best_certificate.json", certificate)
        return certificate

    with gzip.open(output / "accepted_moves.jsonl.gz", "wt") as log:
        try:
            while time.monotonic() - started < seconds:
                if is_standard(tri):
                    best_end, best_moves, best_size = states[-1], list(moves), tri.size()
                    break
                choices = []
                if mode == "reduce":
                    for kind, dimension in (("C", 1), ("P", 0), ("P", 1),
                                             ("Z", 2), ("Z", 1), ("Z", 0)):
                        choices = legal(tri, kind, dimension)
                        if choices:
                            break
                    if not choices:
                        mode = "neutral"
                        neutral_left = 20
                if mode == "neutral" and not choices:
                    choices = legal(tri, "P", 2) + legal(tri, "F", 1)
                    neutral_left -= 1
                    if not choices or neutral_left < 0:
                        mode = "grow"
                        growth_left = rng.randint(1, 6)
                        choices = []
                if mode == "grow" and not choices:
                    if tri.size() >= best_size + growth:
                        states = [start] + [row["after"] for row in best_moves]
                        moves = list(best_moves)
                        positions = {s: i for i, s in enumerate(states)}
                        tri = decode(best_end, state_format)
                        mode = "reduce"
                        continue
                    choices = legal(tri, "P", 3)
                    if not choices and tri.size() + 4 <= best_size + growth:
                        choices = [["P", 4, rng.randrange(tri.size())]]
                    if not choices:
                        break
                    growth_left -= 1
                    if growth_left <= 0:
                        mode = "reduce"
                operation = rng.choice(choices)
                before = states[-1]
                perform(tri, operation)
                after = encode(tri, state_format)
                accepted += 1
                max_size = max(max_size, tri.size())
                row = {"before": before, "operation": operation, "after": after}
                log.write(json.dumps(row, separators=(",", ":")) + "\n")
                log.flush()
                if after in positions:
                    old = positions[after]
                    for removed in states[old + 1:]:
                        del positions[removed]
                    del states[old + 1:]
                    del moves[old:]
                else:
                    moves.append(row)
                    positions[after] = len(states)
                    states.append(after)
                # This search uses labelled tight encodings, so retaining tri
                # preserves the exact recorded state and Regina's valid caches.
                # Actual rollbacks above still decode the saved checkpoint.
                if tri.size() < best_size or is_standard(tri):
                    best_end, best_moves, best_size = after, list(moves), tri.size()
                    mode = "reduce"
                if time.monotonic() - saved > 10:
                    checkpoint()
                    saved = time.monotonic()
                if time.monotonic() - last_progress > 10:
                    print(json.dumps({"accepted": accepted, "current": tri.size(),
                                      "best": best_size, "seconds": time.monotonic() - started}), flush=True)
                    last_progress = time.monotonic()
        except KeyboardInterrupt:
            interrupted = True
    certificate = checkpoint()
    verification = replay(certificate, source_file)
    save_json(output / "REPLAY.json", verification)
    if verification["standard_endpoint"]:
        save_json(output / "SUCCESS.json", certificate)
    print(json.dumps({k: verification[k] for k in
                      ("status", "moves", "start_pentachora", "end_pentachora")}), flush=True)
    return verification


def selftest():
    pinned()
    original = standard_sphere()
    original_sig = original.isoSig()
    grown = regina.Triangulation4.fromIsoSig(original_sig)
    perform(grown, ["P", 4, 0])
    grown_sig = grown.isoSig()
    inverse = None
    for move in legal(regina.Triangulation4.fromIsoSig(grown_sig), "P", 0):
        test = regina.Triangulation4.fromIsoSig(grown_sig)
        perform(test, move)
        if test.isoSig() == original_sig:
            inverse = move
            break
    require(inverse is not None, "subdivision has no inverse move")
    certificate = {"start": grown_sig, "end": original_sig, "standard_endpoint": True,
                   "moves": [{"before": grown_sig, "operation": inverse, "after": original_sig}]}
    require(replay(certificate)["standard_endpoint"], "valid inverse proof rejected")
    bad = json.loads(json.dumps(certificate))
    bad["moves"][0]["operation"][2] = 10**9
    try:
        replay(bad)
    except ValueError:
        pass
    else:
        raise ValueError("corrupt face accepted")
    bad = json.loads(json.dumps(certificate))
    bad["standard_endpoint"] = False
    try:
        replay(bad)
    except ValueError:
        pass
    else:
        raise ValueError("false endpoint claim accepted")
    try:
        admit(regina.Example4.cp2())
    except ValueError:
        pass
    else:
        raise ValueError("CP2 passed homology-sphere filter")
    print("PASS: explicit sphere inverse path; corrupt index/endpoint flags and CP2 rejected")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    s = commands.add_parser("search")
    s.add_argument("source", type=Path)
    s.add_argument("--out", required=True, type=Path)
    s.add_argument("--seconds", type=float, default=120)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--growth", type=int, default=24)
    r = commands.add_parser("replay")
    r.add_argument("certificate", type=Path)
    r.add_argument("--source", type=Path)
    r.add_argument("--output", type=Path)
    commands.add_parser("selftest")
    args = parser.parse_args()
    if args.command == "search":
        search(args.source, args.out, args.seconds, args.seed, args.growth)
    elif args.command == "replay":
        result = replay(json.loads(args.certificate.read_text()), args.source)
        if args.output:
            save_json(args.output, result)
        print(json.dumps({k: result[k] for k in
                          ("status", "moves", "start_pentachora", "end_pentachora",
                           "source_binding_checked")}), flush=True)
    else:
        selftest()

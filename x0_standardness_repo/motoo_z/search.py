#!/usr/bin/env python3
"""Run retained guarded simplification producers with explicit input binding.

small: parallel excursions from a canonical signature, with fresh independent
       replay of the resulting small certificate. The default is the published
       eight-pentachoron endpoint, not the old script's redraw default.
large: sampled excursions on a labelled closed triangulation. Produces a
       chunked move journal that requires a separate independent replay.

No size threshold or user-provided list is treated as proof of standardness.
"""
from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
import sys
import time

from common import (ROOT, LEGACY, REFERENCE, check_integrity, fresh_output,
                    require, run_stage, sha256, write_json)
from verify import (START8, admit_small, helper_command, prepare_recognition_helper,
                    read_json, recognition_runtime_check, replay_terminal_certificate,
                    runtime)


def load_triangulation_file(path, regina):
    """Accept a root Triangulation4 packet or one unambiguous descendant."""
    path = Path(path).expanduser().resolve()
    require(path.is_file(), "Missing triangulation input: " + str(path))
    if path.suffix.lower() in (".sig", ".isosig"):
        return regina.Triangulation4.fromIsoSig(path.read_text().strip()), False
    packet = regina.open(str(path))
    require(packet is not None, "Regina could not open the input file")
    if isinstance(packet, regina.Triangulation4):
        return regina.Triangulation4(packet), True
    stack, candidates = [packet], []
    while stack:
        node = stack.pop()
        if isinstance(node, regina.Triangulation4):
            candidates.append(node)
        child = node.firstChild()
        while child is not None:
            stack.append(child)
            child = child.nextSibling()
    require(len(candidates) == 1,
            "Input container must have exactly one 4D triangulation; export the intended packet separately")
    return regina.Triangulation4(candidates[0]), False


def resolve_start(args, regina):
    if args.start_file:
        source = args.start_file.expanduser().resolve()
        digest = sha256(source)
        tri, root_packet = load_triangulation_file(source, regina)
        require(sha256(source) == digest, "Source file changed while loading")
        return tri, {"kind": "file", "path": str(source), "sha256": digest,
                     "root_triangulation_packet": root_packet}
    signature = args.start or START8
    tri = regina.Triangulation4.fromIsoSig(signature)
    return tri, {"kind": "isomorphism_signature", "supplied_signature": signature,
                 "default_published_endpoint": args.start is None}


def prepare_known(path, out, regina):
    if path is None:
        return None
    path = path.expanduser().resolve()
    before = sha256(path)
    signatures = set()
    for line in path.read_text().splitlines():
        for value in line.split("#", 1)[0].split():
            signatures.add(regina.Triangulation4.fromIsoSig(value).isoSig())
    require(sha256(path) == before, "Known-target list changed while loading")
    require(signatures, "Known-target list is empty")
    frozen = out / "KNOWN_TARGETS.txt"
    frozen.write_text("\n".join(sorted(signatures)) + "\n")
    return {"original_path": str(path), "original_sha256": before,
            "frozen_path": frozen.name, "frozen_sha256": sha256(frozen),
            "signatures": sorted(signatures),
            "standardness_of_listed_targets_assumed": False}


def new_search(args, out, regina, integrity):
    tri, binding = resolve_start(args, regina)
    meta = {"schema": "motoo-z-search-input-v1", "phase": args.phase,
            "input": binding, "release_manifest_sha256": integrity["manifest_sha256"],
            "regina_distribution": "7.4.1", "created_unix_time": time.time()}
    if args.phase == "small":
        meta["start"] = tri.isoSig()
        meta["source_admission"] = admit_small(tri)
        meta["producer_sha256"] = sha256(REFERENCE / "code/simplify_penta.py")
        meta["known_targets"] = prepare_known(args.known, out, regina)
    else:
        if binding["kind"] == "file" and binding["root_triangulation_packet"]:
            source = Path(binding["path"])
        else:
            source = out / "SOURCE.rga"
            require(regina.PacketOfTriangulation4(tri).save(str(source)), "Could not save selected source packet")
        meta["source"] = {"path": str(source), "sha256": sha256(source),
                          "tight_sha256": hashlib.sha256(tri.tightEncoding().encode("ascii")).hexdigest(),
                          "pentachora": tri.size()}
        meta["source_admission"] = "Performed afresh by the guarded large runner, optionally inside the helper"
        meta["producer_sha256"] = sha256(LEGACY / "reduce_closed_sampled.py")
        meta["guard_journal_sha256"] = sha256(LEGACY / "reduce_closed_bursts.py")
    write_json(out / "SEARCH_INPUT.json", meta)
    return meta


def small(args, out, meta, report):
    require(args.neutral >= 0 and args.pool > 0 and args.stuck > 0 and args.kmax >= 2
            and args.connect_at >= 0 and args.connect_extra >= 0 and args.connect_seconds > 0,
            "Invalid small-search policy")
    script = REFERENCE / "code/simplify_penta.py"
    require(meta["producer_sha256"] == sha256(script), "Small-search producer changed since this run began")
    walk = out / "walk"
    command = [sys.executable, "-B", script, "search",
               "--hours", str(args.hours), "--workers", str(args.workers),
               "--neutral", str(args.neutral), "--pool", str(args.pool),
               "--stuck", str(args.stuck), "--kmax", str(args.kmax),
               "--connect-at", str(args.connect_at), "--connect-extra", str(args.connect_extra),
               "--connect-seconds", str(args.connect_seconds)]
    if args.resume:
        state = read_json(walk / "STATE.json")
        require(state["start"] == meta["start"], "Resume state starts at a different triangulation")
        command.extend(["--resume", walk])
    else:
        command.extend(["--out", walk, "--start", meta["start"], "--seed", str(args.seed)])
    known = meta.get("known_targets")
    if known:
        frozen = out / known["frozen_path"]
        require(sha256(frozen) == known["frozen_sha256"], "Frozen known-target list changed")
        command.extend(["--known", frozen])
    run_stage("small_search_%03d" % report["attempt"], command, out, report)
    certificate = walk / "CERTIFICATE.json"
    if not certificate.is_file():
        certificate = walk / "BEST_CHAIN.json"
    require(certificate.is_file(), "Producer did not save a complete or partial path")
    replay = replay_terminal_certificate(certificate, meta["start"])
    replay_file = out / ("SEARCH_REPLAY_%03d.json" % report["attempt"])
    write_json(replay_file, replay)
    report["independent_replay"] = {"path": replay_file.name, "sha256": sha256(replay_file),
                                    "moves": replay["moves"], "start": replay["start"], "end": replay["end"]}
    report["standard_sphere_proved"] = replay["standard_sphere_proved"]
    if replay["standard_sphere_proved"]:
        report["status"] = "PASS_SEARCH_PATH_TO_EXPLICIT_STANDARD_SPHERE"
    elif known and replay["end"] in known["signatures"]:
        report["status"] = "PASS_SEARCH_PATH_TO_USER_LISTED_TARGET_ONLY"
    else:
        report["status"] = "PASS_PARTIAL_SEARCH_PATH"
    report["scope"] = "Only the explicitly bound search start; no new construction-to-start identification is inferred"


def large(args, out, meta, report):
    require(args.growth >= 4 and args.burst24 > 0 and args.mix >= 0
            and args.birth_every > 0 and args.birth_vertices >= 2
            and args.chunk_size > 0 and args.probes > 0, "Invalid large-search policy")
    source = Path(meta["source"]["path"])
    require(sha256(source) == meta["source"]["sha256"], "Selected large-search source changed")
    command = [LEGACY / "reduce_closed_sampled.py", source, "--out", out / "walk",
               "--seconds", str(3600 * args.hours), "--seed", str(args.seed),
               "--growth", str(args.growth), "--burst24", str(args.burst24),
               "--mix", str(args.mix), "--birth-every", str(args.birth_every),
               "--birth-vertices", str(args.birth_vertices), "--chunk-size", str(args.chunk_size),
               "--probes", str(args.probes)]
    if args.edge20_first:
        command.append("--edge20-first")
    if args.initial_burst:
        command.append("--initial-burst")
    if args.recognition_helper:
        bridge = prepare_recognition_helper(out, report)
        argv = helper_command(bridge, out / "recognition", command)
    else:
        argv = [sys.executable, "-B", *command]
    run_stage("large_search", argv, out, report)
    if args.recognition_helper:
        report["recognition"] = recognition_runtime_check(out / "recognition")
    checkpoint = read_json(out / "walk/CHECKPOINT.json")
    recorded_source = read_json(out / "walk/SOURCE.json")
    require(recorded_source["source_sha256"] == meta["source"]["sha256"]
            and recorded_source["start_tight_sha256"] == meta["source"]["tight_sha256"],
            "Large-search journal starts at a different source")
    report["status"] = "SEARCH_JOURNAL_RECORDED_REPLAY_REQUIRED"
    report["standard_sphere_proved"] = False
    report["producer_checkpoint"] = {key: checkpoint[key] for key in
                                      ("status", "moves", "current_pentachora", "best_pentachora",
                                       "standard_endpoint", "native_4d_simplify_used")}
    report["fresh_independent_replay_performed"] = False
    report["replay_command"] = [str(sys.executable), "-B", str(LEGACY / "verify_closed_bursts.py"),
                                str(out / "walk"), "--source", str(source),
                                "--out", str(out / "INDEPENDENT_REPLAY.json")]
    report["scope"] = "A discovery journal, including any producer target flag, requires independent move replay"


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    sub = result.add_subparsers(dest="phase", required=True)
    for phase in ("small", "large"):
        p = sub.add_parser(phase)
        source = p.add_mutually_exclusive_group(required=phase == "large")
        source.add_argument("--start", help="Starting isomorphism signature")
        source.add_argument("--start-file", type=Path, help="Root/container RGA or .sig/.isosig text file")
        output = p.add_mutually_exclusive_group(required=True)
        output.add_argument("--output", type=Path)
        if phase == "small":
            output.add_argument("--resume", type=Path, help="Resume a prior wrapper output directory")
        else:
            p.set_defaults(resume=None)
        p.add_argument("--hours", type=float, default=8)
        p.add_argument("--seed", type=int, default=1)
        if phase == "small":
            p.add_argument("--workers", type=int, default=0)
            p.add_argument("--known", type=Path, help="Optional target signatures; membership is not a standardness proof")
            p.add_argument("--neutral", type=int, default=40)
            p.add_argument("--pool", type=int, default=400)
            p.add_argument("--stuck", type=int, default=20000)
            p.add_argument("--kmax", type=int, default=14)
            p.add_argument("--connect-at", type=int, default=4)
            p.add_argument("--connect-extra", type=int, default=6)
            p.add_argument("--connect-seconds", type=int, default=900)
        else:
            p.add_argument("--recognition-helper", action="store_true")
            p.add_argument("--growth", type=int, default=40)
            p.add_argument("--burst24", type=int, default=12)
            p.add_argument("--mix", type=int, default=80)
            p.add_argument("--birth-every", type=int, default=3)
            p.add_argument("--birth-vertices", type=int, default=5)
            p.add_argument("--chunk-size", type=int, default=1000)
            p.add_argument("--probes", type=int, default=64)
            p.add_argument("--edge20-first", action=argparse.BooleanOptionalAction, default=True)
            p.add_argument("--initial-burst", action="store_true")
    return result


def main():
    args = parser().parse_args()
    require(math.isfinite(args.hours) and args.hours > 0, "--hours must be finite and positive")
    if args.phase == "small":
        require(args.workers >= 0, "--workers must be nonnegative")
    integrity = check_integrity()
    regina = runtime()
    if args.resume:
        require(not args.start and not args.start_file and not args.known,
                "A resumed search retains its original start and target list")
        out = args.resume.expanduser().resolve()
        meta = read_json(out / "SEARCH_INPUT.json")
        require(meta["schema"] == "motoo-z-search-input-v1" and meta["phase"] == "small",
                "Resume requires a previous small-search wrapper directory")
        require(meta["release_manifest_sha256"] == integrity["manifest_sha256"],
                "Resume must use the same released module")
        prior = read_json(out / "SEARCH.json")
        require(prior.get("input_metadata_sha256") == sha256(out / "SEARCH_INPUT.json"),
                "Resume input metadata differs from the previous wrapper receipt")
        binding = meta["input"]
        if binding["kind"] == "file":
            require(sha256(binding["path"]) == binding["sha256"],
                    "Original resume input file changed")
            original, _ = load_triangulation_file(binding["path"], regina)
        else:
            require(binding["kind"] == "isomorphism_signature", "Unknown resume input kind")
            original = regina.Triangulation4.fromIsoSig(binding["supplied_signature"])
        require(original.isoSig() == meta["start"],
                "Resume start is not the canonical triangulation of its bound original input")
        attempt = prior.get("attempt", 0) + 1
        report = {"schema": "motoo-z-search-v1", "phase": "small", "status": "RUNNING",
                  "attempt": attempt, "stages": prior.get("stages", [])}
    else:
        out = fresh_output(args.output)
        meta = new_search(args, out, regina, integrity)
        report = {"schema": "motoo-z-search-v1", "phase": args.phase, "status": "RUNNING",
                  "attempt": 1, "stages": []}
    report["release_integrity"] = integrity
    report["requested_additional_hours"] = args.hours
    report["input_metadata_sha256"] = sha256(out / "SEARCH_INPUT.json")
    write_json(out / "SEARCH.json", report)
    started = time.monotonic()
    try:
        if args.phase == "small":
            small(args, out, meta, report)
        else:
            large(args, out, meta, report)
        require(check_integrity() == integrity, "Released inputs changed during search")
        require(sha256(out / "SEARCH_INPUT.json") == report["input_metadata_sha256"],
                "Search input metadata changed")
        if meta["input"]["kind"] == "file":
            require(sha256(meta["input"]["path"]) == meta["input"]["sha256"], "Original source file changed")
        report["seconds_this_attempt"] = time.monotonic() - started
        write_json(out / "SEARCH.json", report)
        print(report["status"], flush=True)
        print("Receipt: " + str(out / "SEARCH.json"), flush=True)
    except BaseException as error:
        report["status"] = "FAILED_OR_INTERRUPTED"
        report["error"] = type(error).__name__ + ": " + str(error)
        report["seconds_this_attempt"] = time.monotonic() - started
        write_json(out / "SEARCH.json", report)
        raise


if __name__ == "__main__":
    main()

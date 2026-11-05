#!/usr/bin/env python3
"""Replay the released Motoo Z closed triangulation-to-sphere certificate.

Terminal scope proves only that the specified eight-pentachoron endpoint is
standard. Full scope also replays the 780,847-move path from the exact closed
24,170-pentachoron input. Neither scope reruns the Kirby-to-closed construction.
All essential checks are explicit runtime failures, including under python -O.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import sys
import tarfile
import time

from common import (ROOT, LEGACY, REFERENCE, check_integrity, fresh_output,
                    require, run_stage, sha256, write_json)

START8 = "iLAALMQkabcccfefhghhhcaaavava2aaaoacbZb6aEa8a8a"
STANDARD = "cPkbbbbaaaaaaaa"
SOURCE_SHA256 = "e6281897c358e315757a77f851433db1527f48d71b97ec0056174e0a7eb2c467"
INITIAL_TIGHT_SHA256 = "1f698ed1061bffaf7f380ec5b4ed8075306e6cffd8f6c3cc281873f67003f92d"
FINAL_TIGHT_SHA256 = "43fdb32848802ce7f49df04f4c100410f4a411275f76661426243681fa686bdb"
BEST_RGA_SHA256 = "1e170b5e297d64e95a1438dbbbfee8de9dc5f792280eb6a6c2ec747cc8cd815e"
TERMINAL_SHA256 = "54819a4b09ece5267e08d00e6168b674642c5d2c4b379259db8ff1b62a6d8a91"
SOURCE_METADATA_SHA256 = "1f562bc3bc3f0cd35fdc20784d9f4c7dfe49d739a4cc0f438146276985269da3"
CHECKPOINT_SHA256 = "4b1256a2defd20c77cca2cf8d7119c4c5307aefc89ba339de73fe03dd7e47de2"
MOVES_ARCHIVE_SHA256 = "cf0a077bc3680e331a444b8b4fe1564e5077d00d20ca51d81c610717cb70e4ab"
LONG_VERIFIER_SHA256 = "5fc138e9cd563b8f323814f90a22067d9120ff6432a2db143c916625309f2461"
HELPER_ENGINE_SHA256 = "42e0c5a630d71f98c6e5ff3c70c520d1c3d136151aad4ac37fe807004fbaca1d"
LONG_COUNTS = {"['C', 1]": 1287, "['Z', 1]": 71713,
               "['Z', 2]": 15485, "['P', 1]": 15036,
               "['P', 4]": 8160, "['P', 3]": 82355,
               "['F', 1]": 286182, "['P', 2]": 293738,
               "['Z', 0]": 6891}


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "Duplicate JSON key: " + str(key))
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs)


def runtime():
    require(importlib.metadata.version("regina") == "7.4.1", "Requires regina==7.4.1")
    import regina
    return regina


def standard_double(regina):
    tri = regina.Triangulation4()
    first, second = tri.newSimplices(2)
    for facet in range(5):
        first.join(facet, second, regina.Perm5())
    require(tri.isoSig() == STANDARD, "Explicit standard-double signature differs")
    return tri


def admit_small(tri):
    require(tri.size() > 0 and tri.isValid() and tri.isConnected()
            and tri.isOrientable() and not tri.isIdeal() and tri.isClosed(),
            "Triangulation is not a finite closed connected orientable PL manifold")
    return {"pentachora": tri.size(), "vertices": tri.countVertices(),
            "euler_characteristic": tri.eulerCharTri(),
            "finite_closed_connected_orientable_valid": True}


def checked_move(tri, move):
    require(isinstance(move, list) and len(move) == 3, "Malformed move")
    kind, dimension, index = move
    require(kind in ("P", "20", "44", "C"), "Unsupported terminal move kind")
    require(type(dimension) is int and 0 <= dimension <= 4, "Invalid face dimension")
    count = tri.size() if dimension == 4 else tri.countFaces(dimension)
    require(type(index) is int and 0 <= index < count, "Invalid face index")
    face = tri.simplex(index) if dimension == 4 else tri.face(dimension, index)
    before, vertices = tri.size(), tri.countVertices()
    if kind == "P":
        require(tri.hasPachner(face), "Pachner precondition rejected")
        require(tri.pachner(face), "Pachner application rejected")
        require(tri.size() == before + 2 * dimension - 4, "Wrong Pachner count change")
        if dimension in (0, 4):
            require(tri.countVertices() == vertices + (1 if dimension == 4 else -1),
                    "Wrong extreme-Pachner vertex change")
    elif kind == "20":
        require(dimension in (0, 1, 2), "Invalid 2-0 face dimension")
        require(tri.has20(face), "2-0 precondition rejected")
        require(tri.move20(face) and tri.size() == before - 2, "2-0 application/count rejected")
    elif kind == "44":
        require(dimension == 1 and tri.has44(face), "4-4 precondition rejected")
        require(tri.move44(face) and tri.size() == before, "4-4 application/count rejected")
    else:
        require(dimension == 1 and face.vertex(0) != face.vertex(1),
                "Collapse must use an edge with distinct endpoints")
        require(tri.hasCollapseEdge(face), "Collapse precondition rejected")
        require(tri.collapseEdge(face), "Collapse application rejected")
        require(tri.size() < before and tri.countVertices() == vertices - 1,
                "Wrong collapse count or vertex change")


def replay_terminal_certificate(path, expected_start, require_standard=False):
    """Independent guarded replay, also usable for a newly discovered small path.

An arbitrary complete/target label in the input never proves standardness.
Only an actual isomorphism to the explicitly constructed double does so.
"""
    regina = runtime()
    path = Path(path)
    before_hash = sha256(path)
    cert = read_json(path)
    require(isinstance(cert, dict) and isinstance(cert.get("steps"), list),
            "Malformed terminal certificate")
    require(isinstance(cert.get("start"), str) and cert["start"] == expected_start,
            "Unexpected terminal start")
    require(isinstance(cert.get("end"), str), "Missing terminal endpoint")
    if "complete" in cert:
        require(type(cert["complete"]) is bool, "Invalid completeness flag")
    tri = regina.Triangulation4.fromIsoSig(cert["start"])
    require(tri.isoSig() == cert["start"], "Start signature is not canonical")
    initial = admit_small(tri)
    cur = cert["start"]
    counts, directions = Counter(), Counter()
    moves = 0
    maximum = tri.size()
    segments = []
    for number, step in enumerate(cert["steps"]):
        require(isinstance(step, dict), "Malformed terminal segment")
        require(set(step) <= {"from", "to", "moves", "direction"}, "Unknown segment field")
        require(step.get("from") == cur and isinstance(step.get("to"), str),
                "Broken terminal segment chain")
        direction = step.get("direction", "forward")
        require(direction in ("forward", "backward"), "Unknown segment direction")
        require(isinstance(step.get("moves"), list) and step["moves"], "Empty/malformed segment")
        source, target = ((step["from"], step["to"]) if direction == "forward"
                          else (step["to"], step["from"]))
        tri = regina.Triangulation4.fromIsoSig(source)
        require(tri.isoSig() == source, "Segment source is not canonical")
        source_size = tri.size()
        source_admission = admit_small(tri)
        require(source_admission["euler_characteristic"] == initial["euler_characteristic"],
                "Segment source Euler characteristic differs")
        maximum = max(maximum, source_size)
        for move in step["moves"]:
            checked_move(tri, move)
            counts[str(move[:2])] += 1
            directions[direction] += 1
            moves += 1
            maximum = max(maximum, tri.size())
        require(tri.isoSig() == target, "Segment does not produce the recorded exact endpoint")
        state = admit_small(tri)
        require(state["euler_characteristic"] == initial["euler_characteristic"],
                "Euler characteristic changed")
        segments.append({"index": number, "direction": direction,
                         "moves": len(step["moves"]), "applied_source_pentachora": source_size,
                         "applied_target_pentachora": tri.size()})
        cur = step["to"]
    require(cur == cert["end"], "Wrong final endpoint")
    endpoint = regina.Triangulation4.fromIsoSig(cur)
    require(endpoint.isoSig() == cur, "Final signature is not canonical")
    final = admit_small(endpoint)
    require(final["euler_characteristic"] == initial["euler_characteristic"],
            "Endpoint Euler characteristic changed")
    if "end_size" in cert:
        require(type(cert["end_size"]) is int and cert["end_size"] == endpoint.size(),
                "Wrong recorded endpoint size")
    if "total_moves" in cert:
        require(type(cert["total_moves"]) is int and cert["total_moves"] == moves,
                "Wrong total move count")
    standard = standard_double(regina)
    proved_standard = endpoint.size() == 2 and endpoint.isIsomorphicTo(standard) is not None
    if require_standard:
        require(proved_standard and cur == STANDARD, "The path does not reach the explicit standard sphere")
    require(sha256(path) == before_hash, "Certificate changed during replay")
    return {"schema": "motoo-z-independent-terminal-replay-v1",
            "status": "PASS_STANDARD_S4" if proved_standard else "PASS_PARTIAL_PATH",
            "certificate_sha256": before_hash, "start": cert["start"], "end": cur,
            "moves": moves, "move_counts": dict(counts), "directions": dict(directions),
            "segments": segments, "maximum_pentachora": maximum,
            "source_admission": initial, "final_admission": final,
            "standard_reference": STANDARD, "explicit_identity_glued_double_constructed": True,
            "standard_sphere_proved": proved_standard,
            "all_moves_guarded": True, "discovery_code_imported": False,
            "regina_distribution": importlib.metadata.version("regina"),
            "trust": "Regina 7.4.1 local-move guards, skeleton, signature codec and isomorphism checks"}


def recognition_runtime_check(proofs):
    result = read_json(Path(proofs) / "RUNTIME.json")
    require(result["status"] == "TARGET_COMPLETED", "Recognition target did not complete normally")
    require(result["original_engine_sha256"] == HELPER_ENGINE_SHA256,
            "Recognition helper used an unreviewed engine")
    require(result["four_dimensional_guard_kernel_modified"] is False,
            "Recognition helper changed the 4D guard kernel")
    require(result["pinned_files_unchanged_at_exit"] and
            all(result["pinned_files_unchanged_at_exit"].values()), "Recognition pins changed")
    return {"runtime_sha256": sha256(Path(proofs) / "RUNTIME.json"),
            "original_engine_sha256": result["original_engine_sha256"],
            "calls": len(result["calls"]), "counters": result["bridge_counters"],
            "four_dimensional_guard_kernel_modified": False}


def prepare_recognition_helper(out, report):
    """Build and test the unchanged helper for its exact admitted Linux wheel."""
    regina = runtime()
    require(sys.platform == "linux" and platform.machine() == "x86_64",
            "The optional recognition helper requires Linux x86-64")
    engines = list(Path(regina.__file__).parent.glob("engine*.so"))
    require(len(engines) == 1 and sha256(engines[0]) == HELPER_ENGINE_SHA256,
            "Optional helper requires its pinned CPython 3.12 Regina wheel; do not change the pin")
    compiler = shutil.which("c++") or shutil.which("g++")
    require(compiler is not None, "A C++ compiler is required for the optional helper")
    tools = Path(out) / "helper"
    tools.mkdir()
    bridge = tools / "primal_recognition_bridge.so"
    run_stage("compile_recognition_helper",
              [compiler, "-std=c++17", "-O2", "-fPIC", "-shared",
               LEGACY / "primal_recognition_bridge.cpp",
               "-Wl,--version-script=" + str(LEGACY / "primal_recognition_bridge.map"),
               "-o", bridge, "-ldl"], out, report)
    controls = tools / "controls"
    run_stage("recognition_controls",
              helper_command(bridge, controls, [LEGACY / "test_primal_recognition_bridge.py"]),
              out, report)
    result = read_json(controls / "CONTROLS.json")
    require(result["status"] == "PASS_EXACT_SERIALIZATION_POSITIVE_AND_ORIGINAL_FALLBACK_CONTROLS",
            "Recognition controls did not pass")
    require(result["native_4d_internal_calls_verified"] == 5,
            "Recognition controls omitted the internal 4D link calls")
    report["recognition_controls"] = recognition_runtime_check(controls)
    return bridge


def helper_command(bridge, proofs, command):
    return [sys.executable, "-B", LEGACY / "run_with_primal_recognition.py",
            "--bridge", bridge, "--proofs", proofs, "--", *command]


def extract_journal(archive, out):
    destination = Path(out) / "move_archive"
    destination.mkdir()
    seen, total = set(), 0
    with tarfile.open(archive, "r:xz") as stream:
        for item in stream:
            name = PurePosixPath(item.name)
            require(not name.is_absolute() and ".." not in name.parts and name.parts
                    and name.parts[0] == "best8_prefix", "Unsafe move-archive path")
            require(item.isdir() or item.isfile(), "Move archive contains a link or special file")
            require(str(name) not in seen, "Duplicate move-archive member")
            seen.add(str(name))
            target = destination.joinpath(*name.parts)
            if item.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            total += item.size
            require(total <= 300 * 1024 * 1024, "Unexpected move-archive expansion")
            target.parent.mkdir(parents=True, exist_ok=True)
            source = stream.extractfile(item)
            require(source is not None, "Cannot read archive member")
            with source, target.open("xb") as writer:
                shutil.copyfileobj(source, writer)
    journal = destination / "best8_prefix"
    require(sha256(journal / "SOURCE.json") == SOURCE_METADATA_SHA256,
            "Wrong extracted source metadata")
    require(sha256(journal / "CHECKPOINT.json") == CHECKPOINT_SHA256,
            "Wrong extracted checkpoint")
    return journal


def bind_selected_source(selected):
    """Bridge compression-only RGA differences without changing the certificate.

Both the complete decompressed XML and the exact labelled tight encoding must
equal the released source. A mere isomorphism signature match is insufficient.
"""
    regina = runtime()
    selected = Path(selected).expanduser().resolve()
    released = REFERENCE / "closed/closed_candidate.rga"
    require(sha256(released) == SOURCE_SHA256, "Released closed source changed")
    selected_hash = sha256(selected)

    def payload(path):
        raw = Path(path).read_bytes()
        return gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw

    expected_xml, selected_xml = payload(released), payload(selected)
    require(selected_xml == expected_xml, "Selected source XML differs from the released construction")
    original = regina.open(str(released))
    candidate = regina.open(str(selected))
    require(isinstance(original, regina.Triangulation4)
            and isinstance(candidate, regina.Triangulation4), "Source is not a root 4D triangulation")
    expected_encoding, selected_encoding = original.tightEncoding(), candidate.tightEncoding()
    require(selected_encoding == expected_encoding, "Selected labelled triangulation differs")
    require(hashlib.sha256(expected_encoding.encode("ascii")).hexdigest() == INITIAL_TIGHT_SHA256,
            "Released source does not have the certified initial numbered state")
    require(sha256(selected) == selected_hash and sha256(released) == SOURCE_SHA256,
            "Source changed during identity comparison")
    replay_source = selected if selected_hash == SOURCE_SHA256 else released
    bridge = {"status": "PASS_EXACT_XML_AND_LABELLED_STATE_IDENTITY",
              "selected_path": str(selected), "selected_sha256": selected_hash,
              "released_path": str(released), "released_sha256": SOURCE_SHA256,
              "replay_source_path": str(replay_source), "replay_source_sha256": SOURCE_SHA256,
              "compressed_bytes_equal": selected_hash == SOURCE_SHA256,
              "complete_xml_payload_equal": True,
              "xml_payload_sha256": hashlib.sha256(expected_xml).hexdigest(),
              "exact_tight_encoding_equal": True, "initial_tight_sha256": INITIAL_TIGHT_SHA256,
              "selected_source_file_modified": False, "certificate_rewritten": False}
    return replay_source, bridge


def check_long_receipt(receipt, source, journal):
    require(receipt["status"] == "PASS_PARTIAL_CLOSED_PATH" and
            receipt["standard_sphere_proved"] is False, "Wrong long-path receipt scope")
    require(receipt["source_sha256"] == sha256(source) == SOURCE_SHA256,
            "Long replay used a different source")
    require(receipt["source_metadata_sha256"] == SOURCE_METADATA_SHA256
            and receipt["checkpoint_sha256"] == CHECKPOINT_SHA256
            and receipt["verifier_sha256"] == LONG_VERIFIER_SHA256, "Long replay bindings differ")
    require(receipt["moves"] == 780847 and receipt["counts"] == LONG_COUNTS,
            "Long replay move counts differ")
    require(receipt["best_pentachora"] == 8 and receipt["best_at_move"] == 780847
            and receipt["maximum_pentachora"] == 24170
            and receipt["final_tight_sha256"] == FINAL_TIGHT_SHA256, "Long replay endpoint differs")
    for key, size, vertices in (("source_admission", 24170, 19), ("final_admission", 8, 1)):
        require(receipt[key]["finite_closed_connected_orientable_valid"] is True
                and receipt[key]["pentachora"] == size and receipt[key]["vertices"] == vertices,
                "Incomplete manifold admission: " + key)
    require(receipt["all_bound_files_unchanged_through_replay"] is True
            and receipt["runner_imported"] is False
            and receipt["native_4d_simplify_used"] is False, "Long replay omitted required guarantees")
    regina = runtime()
    best = journal / "states/best.rga"
    require(sha256(best) == BEST_RGA_SHA256, "Wrong saved best triangulation")
    packet = regina.open(str(best))
    require(isinstance(packet, regina.Triangulation4), "Best state is not a 4D triangulation")
    require(hashlib.sha256(packet.tightEncoding().encode("ascii")).hexdigest() == FINAL_TIGHT_SHA256,
            "Best labelled state differs")
    require(packet.isoSig() == START8, "Long and terminal paths do not meet")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scope", choices=("full", "terminal"), default="full")
    parser.add_argument("--source", type=Path,
                        help="Fresh closed_candidate.rga; exact XML and labelled state are required; compression may differ")
    parser.add_argument("--recognition-helper", action="store_true",
                        help="Use the unchanged audited helper with its exact Linux CPython 3.12 wheel")
    args = parser.parse_args()
    require(args.scope == "full" or (not args.source and not args.recognition_helper),
            "--source and --recognition-helper apply only to full scope")
    integrity = check_integrity()
    out = fresh_output(args.output)
    report = {"schema": "motoo-z-composed-replay-v1", "status": "RUNNING", "scope": args.scope,
              "release_integrity": integrity, "construction_replayed": False,
              "diagram_identification_reproved": False, "long_path_replayed": False,
              "source_selection": "explicit" if args.source else "released_reference",
              "stages": []}
    started = time.monotonic()
    write_json(out / "VERIFICATION.json", report)
    try:
        runtime()
        certificate = REFERENCE / "final/CERTIFICATE.json"
        require(sha256(certificate) == TERMINAL_SHA256, "Wrong released terminal certificate")
        terminal = replay_terminal_certificate(certificate, START8, require_standard=True)
        require(terminal["moves"] == 430 and terminal["directions"] == {"forward": 423, "backward": 7},
                "Wrong released terminal path counts")
        write_json(out / "TERMINAL_REPLAY.json", terminal)
        report["terminal"] = terminal
        if args.scope == "terminal":
            report["status"] = "PASS_TERMINAL_STANDARDIZATION_ONLY"
            report["claim"] = "The specified eight-pentachoron endpoint is PL homeomorphic to S4."
        else:
            source = (args.source or REFERENCE / "closed/closed_candidate.rga").expanduser().resolve()
            replay_source, identity_bridge = bind_selected_source(source)
            archive = ROOT / "certificates/published_fig2_closed_path_moves.tar.xz"
            require(sha256(archive) == MOVES_ARCHIVE_SHA256, "Wrong long-move archive")
            require(sha256(LEGACY / "verify_closed_bursts.py") == LONG_VERIFIER_SHA256,
                    "The retained independent long-path verifier changed")
            report["source"] = {"path": str(source), "sha256": identity_bridge["selected_sha256"],
                                "initial_tight_sha256": INITIAL_TIGHT_SHA256}
            report["replay_source"] = {"path": str(replay_source), "sha256": SOURCE_SHA256}
            report["source_identity_bridge"] = identity_bridge
            journal = extract_journal(archive, out)
            command = [LEGACY / "verify_closed_bursts.py", journal,
                       "--source", replay_source, "--out", out / "LONG_REPLAY.json"]
            if args.recognition_helper:
                bridge = prepare_recognition_helper(out, report)
                argv = helper_command(bridge, out / "recognition", command)
            else:
                argv = [sys.executable, "-B", *command]
            run_stage("long_path", argv, out, report)
            long_receipt = read_json(out / "LONG_REPLAY.json")
            check_long_receipt(long_receipt, replay_source, journal)
            require(sha256(source) == identity_bridge["selected_sha256"],
                    "Selected fresh source changed during replay")
            if args.recognition_helper:
                report["recognition"] = recognition_runtime_check(out / "recognition")
            report["long_path_replayed"] = True
            report["long_path"] = long_receipt
            report["composed_moves"] = 781277
            report["interface"] = {"tight_sha256": FINAL_TIGHT_SHA256, "isomorphism_signature": START8,
                                   "long_endpoint_equals_terminal_start": True}
            report["status"] = "PASS_CLOSED_SOURCE_STANDARDIZATION"
            report["claim"] = "The exact selected 24170-pentachoron closed source is PL homeomorphic to S4."
        after = check_integrity()
        require(after == integrity, "Released inputs changed during verification")
        report["released_inputs_unchanged"] = True
        report["seconds"] = time.monotonic() - started
        write_json(out / "VERIFICATION.json", report)
        print(report["status"], flush=True)
        print("Receipt: " + str(out / "VERIFICATION.json"), flush=True)
    except BaseException as error:
        report["status"] = "FAILED_OR_INTERRUPTED"
        report["error"] = type(error).__name__ + ": " + str(error)
        report["seconds"] = time.monotonic() - started
        write_json(out / "VERIFICATION.json", report)
        raise


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit the published Z-circle file bindings and journal structure.

This deliberately DOES NOT replay any four-dimensional move or recognise the
source's vertex links.  PASS_STRUCTURAL_BINDINGS_ONLY certifies file identity,
record continuity, and the interfaces between the supplied certificates; it
does not itself establish a PL homeomorphism or standardness.

The large tar archive is read once as a stream.  Only its tiny best.rga member
is temporarily written, so Regina can read it.  All other members remain in
memory one at a time.  Input paths are relative to this file, not the caller's
working directory.  Requires the pinned Regina 7.4.1 Python distribution.

Usage:
    python audit_bindings.py --output runs/bindings.json
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import re
import sys
import tarfile
import tempfile
import time

from common import ROOT, REFERENCE, require, sha256


ARCHIVE = ROOT / "certificates/published_fig2_closed_path_moves.tar.xz"
ARCHIVE_SHA256 = "cf0a077bc3680e331a444b8b4fe1564e5077d00d20ca51d81c610717cb70e4ab"
PAPER_FILES = {
    "diagram/fig2_published_katie.json":
        "5aeaea4ffa70b3d4084fbae11cd144c256ae3af361cbebb10da07166cad3427f",
    "katie/raw_katie.rga":
        "f1590ec6de3e7a9d9fdeb2175f98a9cd70688027776d274de0f2da4718a119bb",
    "closed/closed_candidate.rga":
        "e6281897c358e315757a77f851433db1527f48d71b97ec0056174e0a7eb2c467",
    "final/CERTIFICATE.json":
        "54819a4b09ece5267e08d00e6168b674642c5d2c4b379259db8ff1b62a6d8a91",
}
START_TIGHT = "1f698ed1061bffaf7f380ec5b4ed8075306e6cffd8f6c3cc281873f67003f92d"
END_TIGHT = "43fdb32848802ce7f49df04f4c100410f4a411275f76661426243681fa686bdb"
Z8 = "iLAALMQkabcccfefhghhhcaaavava2aaaoacbZb6aEa8a8a"
W4 = "eLAQcabbcdddcaaaga2amamaya"
S2 = "cPkbbbbaaaaaaaa"
LONG_MOVES = 780847
LONG_CHUNKS = 781
TERMINAL_MOVES = 430
METADATA_NAMES = ("SOURCE.json", "CHECKPOINT.json", "REPLAY.json", "EXTRACTION.json")
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
LONG_COUNTS = {
    "['C', 1]": 1287, "['Z', 0]": 6891, "['Z', 1]": 71713,
    "['Z', 2]": 15485, "['P', 1]": 15036, "['P', 2]": 293738,
    "['P', 3]": 82355, "['P', 4]": 8160, "['F', 1]": 286182,
}
FORWARD_COUNTS = {
    "['P', 1]": 7, "['P', 2]": 369, "['P', 3]": 20,
    "['P', 4]": 3, "['20', 1]": 13, "['20', 2]": 7,
    "['44', 1]": 3, "['C', 1]": 1,
}
BACKWARD_COUNTS = {"['C', 1]": 2, "['P', 2]": 1, "['P', 3]": 4}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: " + key)
        result[key] = value
    return result


def reject_constant(value):
    raise RuntimeError("Non-finite JSON numeric constant: " + value)


def decode(data, label):
    try:
        return json.loads(data, object_pairs_hook=unique_pairs,
                          parse_constant=reject_constant)
    except (ValueError, UnicodeError, RuntimeError) as error:
        raise RuntimeError("Invalid JSON in " + label + ": " + str(error)) from error


def exact_integer(value, minimum=0):
    return type(value) is int and value >= minimum


def checked_input(path):
    path = Path(path)
    require(path.is_file(), "Missing input: " + str(path))
    require(not path.is_symlink() and path.resolve().is_relative_to(ROOT),
            "Input must be a regular file inside this module: " + str(path))
    return path


def operation(value, terminal=False):
    require(isinstance(value, list) and len(value) == 3,
            "A move must be [kind, dimension, index]")
    kind, dimension, index = value
    require(type(kind) is str and exact_integer(dimension) and exact_integer(index),
            "Invalid move field types or negative dimension/index")
    zero = "20" if terminal else "Z"
    neutral = "44" if terminal else "F"
    require((kind == "P" and dimension <= 4) or
            (kind == zero and dimension <= 2) or
            (kind in (neutral, "C") and dimension == 1),
            "Unsupported move family/dimension: " + repr(value))
    return kind, dimension, index


def inspect_chunk(data, item):
    """Check one journal chunk without applying its moves."""
    label = item["file"]
    require(digest(data) == item["sha256"], "Chunk SHA-256 mismatch: " + label)
    rows = data.splitlines()
    require(len(rows) == item["moves"], "Chunk record count mismatch: " + label)
    require(rows, "Empty move chunk: " + label)
    first_hash = first_size = previous_hash = previous_size = None
    counts = Counter()
    minimum = maximum = None
    for offset, line in enumerate(rows):
        record = decode(line, label + ":" + str(offset + 1))
        require(isinstance(record, dict), "Move record must be an object")
        index = record.get("index")
        require(exact_integer(index) and index == item["first_index"] + offset,
                "Nonconsecutive move index in " + label)
        before_hash, after_hash = record.get("before_sha256"), record.get("after_sha256")
        require(type(before_hash) is str and HEX64.fullmatch(before_hash) is not None and
                type(after_hash) is str and HEX64.fullmatch(after_hash) is not None,
                "Invalid tight-encoding hash at move " + str(index))
        before, after = record.get("before_pentachora"), record.get("pentachora")
        require(exact_integer(before, 1) and exact_integer(after, 1) and
                before % 2 == 0 and after % 2 == 0,
                "Invalid closed size at move " + str(index))
        if offset == 0:
            first_hash, first_size = before_hash, before
            minimum, maximum = min(before, after), max(before, after)
        else:
            require(before_hash == previous_hash and before == previous_size,
                    "Hash or size discontinuity at move " + str(index))
        kind, dimension, _ = operation(record.get("operation"))
        delta = after - before
        if kind == "P":
            require(delta == 2 * dimension - 4, "Wrong Pachner size change")
        elif kind == "Z":
            require(delta == -2, "Wrong 2-0 size change")
        elif kind == "F":
            require(delta == 0, "Wrong 4-4 size change")
        else:
            require(delta < 0, "Edge-collapse record does not decrease size")
        counts[str([kind, dimension])] += 1
        previous_hash, previous_size = after_hash, after
        minimum, maximum = min(minimum, after), max(maximum, after)
    return {"first_hash": first_hash, "last_hash": previous_hash,
            "first_size": first_size, "last_size": previous_size,
            "minimum_size": minimum, "maximum_size": maximum,
            "records": len(rows), "counts": counts}


def member_name(member):
    """Reject traversal, duplicates handled by caller, links and special files."""
    name = member.name.rstrip("/") if member.isdir() else member.name
    require(type(name) is str and name and "\\" not in name and "\x00" not in name,
            "Unsafe archive member name")
    require(not PurePosixPath(name).is_absolute() and
            all(part not in ("", ".", "..") for part in name.split("/")),
            "Unsafe archive member path: " + name)
    require(member.isdir() or member.isreg(), "Archive contains a link or special file: " + name)
    require(not getattr(member, "sparse", None), "Sparse archive member is not supported: " + name)
    return name


def inspect_archive(metadata_bytes, checkpoint):
    chunks = checkpoint.get("chunks")
    require(isinstance(chunks, list) and len(chunks) == LONG_CHUNKS,
            "Checkpoint must name exactly 781 chunks")
    expected_chunks = {}
    for number, item in enumerate(chunks):
        require(isinstance(item, dict), "Invalid chunk description")
        expected_file = "moves/" + f"{number:06d}" + ".jsonl"
        expected_count = 1000 if number < LONG_CHUNKS - 1 else 847
        require(item.get("file") == expected_file and
                type(item.get("first_index")) is int and item["first_index"] == 1000 * number and
                type(item.get("moves")) is int and item["moves"] == expected_count and
                type(item.get("sha256")) is str and HEX64.fullmatch(item["sha256"]) is not None,
                "Invalid checkpoint chunk: " + expected_file)
        expected_chunks["best8_prefix/" + expected_file] = item
    expected_metadata = {"best8_prefix/" + name: data for name, data in metadata_bytes.items()}
    best_name = "best8_prefix/states/best.rga"
    expected_files = set(expected_chunks) | set(expected_metadata) | {best_name}
    allowed_directories = {"best8_prefix", "best8_prefix/moves", "best8_prefix/states"}
    seen = set()
    seen_files = set()
    summaries = {}
    best_bytes = None
    uncompressed = 0
    with tarfile.open(ARCHIVE, mode="r|xz") as archive:
        for member in archive:
            name = member_name(member)
            require(name not in seen, "Duplicate archive member: " + name)
            seen.add(name)
            if member.isdir():
                require(name in allowed_directories, "Unexpected archive directory: " + name)
                continue
            require(name in expected_files, "Unexpected archive member: " + name)
            # Every member of this fixed release is less than 300 KiB.  Bound
            # reads defensively without buffering or extracting the archive.
            require(exact_integer(member.size) and member.size <= 2 * 1024 * 1024,
                    "Unexpected archive member size: " + name)
            stream = archive.extractfile(member)
            require(stream is not None, "Cannot read archive member: " + name)
            data = stream.read()
            require(len(data) == member.size, "Short archive member: " + name)
            seen_files.add(name)
            uncompressed += len(data)
            if name in expected_chunks:
                summaries[name] = inspect_chunk(data, expected_chunks[name])
            elif name in expected_metadata:
                require(data == expected_metadata[name],
                        "Archive/reference metadata differ: " + name)
            else:
                best_bytes = data
    require(seen_files == expected_files, "Missing archive members: " + repr(sorted(expected_files - seen_files)))
    require(best_bytes is not None, "No endpoint triangulation in archive")

    previous_hash, previous_size = START_TIGHT, 24170
    minimum = maximum = previous_size
    total = 0
    counts = Counter()
    for item in chunks:
        summary = summaries["best8_prefix/" + item["file"]]
        require(summary["first_hash"] == previous_hash and summary["first_size"] == previous_size,
                "Cross-chunk hash/size discontinuity: " + item["file"])
        previous_hash, previous_size = summary["last_hash"], summary["last_size"]
        total += summary["records"]
        counts.update(summary["counts"])
        minimum = min(minimum, summary["minimum_size"])
        maximum = max(maximum, summary["maximum_size"])
    require(total == LONG_MOVES and previous_hash == END_TIGHT and previous_size == 8,
            "Journal has the wrong final count, hash or size")
    require(dict(counts) == LONG_COUNTS, "Journal move-family totals differ")
    require(minimum == 8 and maximum == 24170, "Journal size extrema differ")
    return best_bytes, {"archive_sha256": ARCHIVE_SHA256,
                        "regular_files": len(seen_files), "chunks": len(chunks),
                        "uncompressed_bytes": uncompressed, "records_scanned": total,
                        "member_safety_checked": True, "chunk_sha256_checked": True,
                        "indices_hashes_sizes_continuous": True,
                        "duplicate_metadata_byte_identical": True,
                        "minimum_pentachora": minimum, "maximum_pentachora": maximum,
                        "operation_counts": dict(counts)}


def check_metadata(metadata, metadata_bytes, best_bytes):
    source, checkpoint = metadata["SOURCE.json"], metadata["CHECKPOINT.json"]
    replay, extraction = metadata["REPLAY.json"], metadata["EXTRACTION.json"]
    for name, schema in (("SOURCE.json", "motoo.closed-bursts.source.v1"),
                         ("CHECKPOINT.json", "motoo.closed-bursts.checkpoint.v1"),
                         ("REPLAY.json", "motoo.closed-bursts.independent-replay.v1")):
        require(metadata[name].get("schema") == schema, "Wrong schema in " + name)
    require(source.get("source_sha256") == replay.get("source_sha256") ==
            PAPER_FILES["closed/closed_candidate.rga"], "Source file binding differs")
    require(source.get("start_tight_sha256") == START_TIGHT, "Wrong start tight hash")
    require(source.get("source_admission", {}).get("pentachora") == 24170,
            "Wrong source pentachoron count")
    require(replay.get("source_metadata_sha256") == extraction.get("source_metadata_sha256") ==
            digest(metadata_bytes["SOURCE.json"]), "Source-metadata hash differs")
    require(replay.get("checkpoint_sha256") == extraction.get("extracted_checkpoint_sha256") ==
            digest(metadata_bytes["CHECKPOINT.json"]), "Checkpoint hash differs")
    require(checkpoint.get("moves") == checkpoint.get("best_moves") ==
            replay.get("moves") == replay.get("best_at_move") == extraction.get("extracted_moves") ==
            LONG_MOVES, "Metadata move count differs")
    require(checkpoint.get("current_tight_sha256") == checkpoint.get("best_tight_sha256") ==
            replay.get("final_tight_sha256") == extraction.get("best_tight_sha256") ==
            END_TIGHT, "Endpoint tight-hash binding differs")
    require(checkpoint.get("current_pentachora") == checkpoint.get("best_pentachora") ==
            replay.get("best_pentachora") == extraction.get("best_pentachora") == 8,
            "Endpoint size metadata differs")
    require(checkpoint.get("current_file") == checkpoint.get("best_file") == "states/best.rga",
            "Checkpoint endpoint filename differs")
    require(checkpoint.get("current_rga_sha256") == checkpoint.get("best_rga_sha256") ==
            extraction.get("best_rga_sha256") == digest(best_bytes),
            "Endpoint RGA byte hash differs")
    require(checkpoint.get("counts") == replay.get("counts") == LONG_COUNTS,
            "Recorded move-family totals differ")
    require(checkpoint.get("max_pentachora") == replay.get("maximum_pentachora") == 24170,
            "Recorded maximum size differs")
    require(replay.get("status") == "PASS_PARTIAL_CLOSED_PATH" and
            replay.get("standard_sphere_proved") is False and
            checkpoint.get("standard_endpoint") is False,
            "Prefix receipt must identify a partial path, not standardness")


def load_root_tri4(regina, path):
    root = regina.open(str(path))
    # These released files are root Triangulation4 packets, not containers.
    require(root is not None and isinstance(root, regina.Triangulation4),
            "Expected a root Triangulation4 packet in " + str(path))
    return root


def check_triangulation_interfaces(regina, best_bytes, certificate):
    closed = load_root_tri4(regina, REFERENCE / "closed/closed_candidate.rga")
    closed_hash = digest(closed.tightEncoding().encode())
    require(closed.size() == 24170 and closed_hash == START_TIGHT,
            "Actual closed triangulation does not match the journal source")
    with tempfile.TemporaryDirectory(prefix="motoo_z_bindings_") as temporary:
        endpoint_path = Path(temporary) / "best.rga"
        endpoint_path.write_bytes(best_bytes)
        endpoint = load_root_tri4(regina, endpoint_path)
        endpoint_hash = digest(endpoint.tightEncoding().encode())
        endpoint_signature = endpoint.isoSig()
        require(endpoint.size() == 8 and endpoint_hash == END_TIGHT and endpoint_signature == Z8,
                "Actual journal endpoint does not match the terminal certificate start")
    require(certificate.get("start") == endpoint_signature, "Terminal start is disconnected")

    standard = regina.Triangulation4()
    first, second = standard.newPentachora(2)
    for facet in range(5):
        first.join(facet, second, regina.Perm5())
    require(standard.isoSig() == S2, "Explicit double has an unexpected signature")
    target = regina.Triangulation4.fromIsoSig(certificate["end"])
    require(target.isIsomorphicTo(standard) is not None,
            "Certificate target is not the explicit double of a pentachoron")
    return {"closed_pentachora": closed.size(), "closed_tight_sha256": closed_hash,
            "journal_endpoint_pentachora": endpoint.size(),
            "journal_endpoint_tight_sha256": endpoint_hash,
            "journal_endpoint_isosig": endpoint_signature,
            "standard_anchor_isosig": standard.isoSig(),
            "standard_anchor_construction": "two pentachora, all five facet gluings identity",
            "standard_anchor_isomorphism_checked": True,
            "manifold_recognition_performed": False}


def inspect_terminal(regina, certificate):
    require(isinstance(certificate, dict), "Terminal certificate must be an object")
    require(certificate.get("start") == Z8 and certificate.get("end") == S2 and
            certificate.get("complete") is True and
            type(certificate.get("end_size")) is int and certificate["end_size"] == 2 and
            type(certificate.get("total_moves")) is int and
            certificate["total_moves"] == TERMINAL_MOVES,
            "Terminal certificate endpoints or declared result differ")
    steps = certificate.get("steps")
    require(isinstance(steps, list) and len(steps) == 26,
            "Terminal certificate must contain 26 blocks")
    current = Z8
    forward, backward = Counter(), Counter()
    signatures = set()
    for number, step in enumerate(steps):
        require(isinstance(step, dict) and step.get("from") == current,
                "Discontinuous terminal block " + str(number))
        direction = step.get("direction", "forward")
        require(direction == ("forward" if number < 19 else "backward"),
                "Unexpected terminal move direction")
        require(type(step.get("to")) is str and step["to"], "Missing terminal target signature")
        signatures.update((current, step["to"]))
        moves = step.get("moves")
        require(isinstance(moves, list) and moves, "Empty or invalid terminal block")
        if direction == "backward":
            require(len(moves) == 1, "Backward terminal blocks must contain one move")
        for move in moves:
            kind, dimension, _ = operation(move, terminal=True)
            (forward if direction == "forward" else backward)[str([kind, dimension])] += 1
        current = step["to"]
        if number == 18:
            require(current == W4, "Unexpected common four-pentachoron triangulation")
    require(current == S2 and dict(forward) == FORWARD_COUNTS and dict(backward) == BACKWARD_COUNTS,
            "Terminal count or endpoint mismatch")
    for signature in signatures:
        triangulation = regina.Triangulation4.fromIsoSig(signature)
        require(triangulation.isoSig() == signature,
                "Terminal signature is not canonical: " + signature)
    return {"blocks": len(steps), "moves_described": sum(forward.values()) + sum(backward.values()),
            "forward_blocks": 19, "forward_moves": sum(forward.values()),
            "backward_blocks": 7, "backward_moves": sum(backward.values()),
            "forward_operation_counts": dict(forward), "backward_operation_counts": dict(backward),
            "common_four_pentachoron_isosig": W4,
            "canonical_signatures_checked": len(signatures), "move_legality_checked": False}


def audit_bindings():
    """Return a fresh STRUCTURAL_ONLY report; raise on any failed binding."""
    started = time.monotonic()
    import regina
    version = importlib.metadata.version("regina")
    require(version == "7.4.1", "Install the pinned Python distribution regina==7.4.1")
    auditor_hash, common_hash = sha256(Path(__file__)), sha256(ROOT / "common.py")
    file_bindings = {}
    for name, expected in PAPER_FILES.items():
        path = checked_input(REFERENCE / name)
        actual = sha256(path)
        require(actual == expected, "Published file hash differs: " + name)
        file_bindings[name] = {"sha256": actual, "bytes": path.stat().st_size}
    checked_input(ARCHIVE)
    require(sha256(ARCHIVE) == ARCHIVE_SHA256, "Published move archive hash differs")
    metadata_bytes = {name: checked_input(REFERENCE / "closed_path" / name).read_bytes()
                      for name in METADATA_NAMES}
    metadata = {name: decode(data, name) for name, data in metadata_bytes.items()}
    best_bytes, journal = inspect_archive(metadata_bytes, metadata["CHECKPOINT.json"])
    check_metadata(metadata, metadata_bytes, best_bytes)
    certificate = decode((REFERENCE / "final/CERTIFICATE.json").read_bytes(), "terminal certificate")
    terminal = inspect_terminal(regina, certificate)
    interfaces = check_triangulation_interfaces(regina, best_bytes, certificate)
    for name, expected in PAPER_FILES.items():
        require(sha256(REFERENCE / name) == expected, "Input changed during audit: " + name)
    for name, data in metadata_bytes.items():
        require((REFERENCE / "closed_path" / name).read_bytes() == data,
                "Metadata changed during audit: " + name)
    require(sha256(ARCHIVE) == ARCHIVE_SHA256, "Archive changed during audit")
    require(sha256(Path(__file__)) == auditor_hash and sha256(ROOT / "common.py") == common_hash,
            "Auditor code changed during audit")
    return {"schema": "motoo-z-structural-bindings-v1", "status": "PASS_STRUCTURAL_BINDINGS_ONLY",
            "scope": "STRUCTURAL_ONLY", "standardness_proved": False,
            "four_dimensional_moves_replayed": 0, "move_legality_checked": False,
            "manifold_recognition_performed": False,
            "paper_bound_files": file_bindings, "journal": journal, "terminal": terminal,
            "triangulation_interfaces": interfaces,
            "total_moves_described": LONG_MOVES + TERMINAL_MOVES,
            "all_bound_inputs_and_code_unchanged": True,
            "auditor_sha256": auditor_hash, "common_sha256": common_hash,
            "python": sys.version, "regina_distribution": version,
            "regina_version_string": regina.versionString(),
            "seconds": time.monotonic() - started,
            "limitations": ["No journal move was applied or checked for admissibility.",
                            "No source vertex link or geometric construction was recognised.",
                            "Historical replay receipts were checked for binding, not substituted for a fresh replay."]}


def output_path(value):
    path = Path(value).expanduser()
    require(not path.is_symlink(), "Output cannot be a symlink")
    path = path.resolve()
    require(not path.exists(), "Output already exists: " + str(path))
    for name in ("legacy", "reference", "certificates", "provenance", "verification", "environment"):
        require(not path.is_relative_to(ROOT / name), "Output would modify released inputs: " + str(path))
    require(not ROOT.is_relative_to(path), "Output path cannot contain the module")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True, type=Path,
                        help="new JSON report filename; existing files are never overwritten")
    args = parser.parse_args()
    try:
        destination = output_path(args.output)
        report = audit_bindings()
        destination.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation also prevents a concurrent process from having its
        # report overwritten between the initial path check and this write.
        with destination.open("x", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
        print(report["status"] + ": 781 chunks; 780847 journal records; 430 terminal moves described.")
        print("No four-dimensional moves replayed; no standardness conclusion from this audit alone.")
        print("Report: " + str(destination))
    except (OSError, ValueError, KeyError, RuntimeError, ImportError, tarfile.TarError,
            importlib.metadata.PackageNotFoundError) as error:
        print("FAIL_STRUCTURAL_BINDINGS: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

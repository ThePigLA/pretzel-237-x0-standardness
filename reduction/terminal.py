#!/usr/bin/env python3
"""Bind checked terminal paths to each other and to a supplied capped endpoint.

This helper checks certificate continuity, receipt/code identities, and the
explicit standard anchor. The independent replay programs execute the moves;
this program does not substitute receipt inspection for those executions.
The caller is responsible for running them into fresh output files first.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
from typing import Any

import regina


CERTIFICATE_SHA256 = "e9ba9f14831bb0da32043d7d93eb0148f0683e2a6e4f67d5156e8e9fc94adcdf"
REPLAY_SHA256 = "f55c0618ab3728dfa458f94e676d452e007d0c74b875942240ca64416d833829"
TWO_ZERO_SHA256 = "81a626d8b3e452ecbb411923a77bf13ed7c169297301cfe85141ed0fb34a98ba"
G6 = "gAMLMQaabcdeffff2a+a2avbmaoaxaOb0aUa"
MEETING = "gLAMPQaccddefefffaaaubaaMafapboapboa"
STANDARD6 = "gLAPPQacbedeeffffaaawaaaobcacaQbPboa"
SPHERE2 = "cPkbbbbaaaaaaaa"
PATHS = {
    "G6_to_X": (G6, MEETING, 6814, 6, {"P2", "P3", "Z1", "Z2"}),
    "std6_to_X": (STANDARD6, MEETING, 1519, 6, {"P2", "P3", "Z1", "Z2"}),
    "sphere_to_std6": (SPHERE2, STANDARD6, 103, 2,
                       {"P0", "P1", "P2", "P3", "P4"}),
}


class TerminalBindingError(ValueError):
    """A supplied certificate, receipt, endpoint, or code identity is invalid."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise TerminalBindingError(message)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _input_path(path: Path, label: str) -> Path:
    try:
        result = Path(path).expanduser().resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise TerminalBindingError(f"{label}: cannot resolve input {path}: {error}") from error
    require(result.is_file(), f"{label}: input is not a regular file: {result}")
    return result


def _json_file(path: Path, label: str) -> tuple[Path, bytes, dict[str, Any]]:
    path = _input_path(path, label)
    raw = path.read_bytes()

    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            require(key not in value, f"{label}: duplicate JSON key {key!r}")
            value[key] = item
        return value

    def bad_constant(value: str) -> None:
        raise TerminalBindingError(f"{label}: non-finite JSON constant {value}")

    try:
        data = json.loads(raw, object_pairs_hook=pairs_hook, parse_constant=bad_constant)
    except (UnicodeError, json.JSONDecodeError) as error:
        raise TerminalBindingError(f"{label}: invalid JSON: {error}") from error
    require(type(data) is dict, f"{label}: top-level JSON value must be an object")
    return path, raw, data


def _mapping(value: Any, label: str) -> dict[str, Any]:
    require(type(value) is dict, f"{label}: expected an object")
    return value


def _counts(value: Any, label: str) -> dict[str, int]:
    result = _mapping(value, label)
    require(all(type(k) is str and type(v) is int and v > 0
                for k, v in result.items()), f"{label}: invalid positive integer counts")
    return result


def _properties(tri: Any) -> dict[str, Any]:
    require(tri.isValid() and tri.isConnected() and tri.isClosed()
            and tri.isOrientable() and tri.eulerCharTri() == 2,
            "terminal endpoint is not valid, connected, closed, orientable with chi=2")
    return {
        "pentachora": tri.size(),
        "vertices": tri.countVertices(),
        "f_vector": list(tri.fVector()),
        "vertex_degrees": sorted(v.degree() for v in tri.vertices()),
    }


def _standard_double() -> Any:
    sphere = regina.Triangulation4()
    first, second = sphere.newPentachora(2)
    for facet in range(5):
        first.join(facet, second, regina.Perm5())
    require(sphere.isoSig() == SPHERE2, "explicit double has unexpected signature")
    require(regina.Example4.sphere().isoSig() == SPHERE2,
            "Regina's sphere example differs from the explicit double")
    return sphere


def _derive_paths(certificate: dict[str, Any]) -> tuple[dict[str, Any], set[str], dict[str, int]]:
    require(set(certificate) == set(PATHS) | {"X", "std6"},
            "terminal certificate has unexpected top-level keys")
    require(certificate["X"] == MEETING and certificate["std6"] == STANDARD6,
            "terminal certificate has different meeting or standard-six metadata")
    derived: dict[str, Any] = {}
    all_states: set[str] = set()
    known_sizes: dict[str, int] = {}
    total_moves: Counter[str] = Counter()
    for label, (start, end, length, initial_size, allowed) in PATHS.items():
        steps = certificate[label]
        require(type(steps) is list and len(steps) == length,
                f"{label}: wrong number or type of steps")
        previous = start
        size = initial_size
        states = {start}
        sizes: Counter[str] = Counter({str(size): 1})
        counts: Counter[str] = Counter()
        degrees: Counter[str] = Counter()
        zero_sizes: list[int] = []
        first_above_ten = None
        require(start not in known_sizes or known_sizes[start] == size,
                f"{label}: repeated state has inconsistent size")
        known_sizes[start] = size
        for index, step in enumerate(steps):
            where = f"{label} step {index}"
            require(type(step) is list and len(step) == 4, f"{where}: malformed step")
            source, code, face_index, destination = step
            require(type(source) is str and source == previous,
                    f"{where}: discontinuity or wrong initial state")
            require(type(destination) is str and bool(destination),
                    f"{where}: invalid destination signature")
            require(type(code) is str and code in allowed, f"{where}: forbidden move code")
            require(type(face_index) is int and face_index >= 0,
                    f"{where}: invalid face index")
            dimension = int(code[1])
            if code.startswith("P"):
                delta, degree = 2 * dimension - 4, 5 - dimension
            else:
                delta, degree = -2, 2
                zero_sizes.append(size)
            counts[code] += 1
            degrees[f"{code}:{degree}"] += 1
            size += delta
            require(size > 0 and size % 2 == 0, f"{where}: impossible closed size")
            require(destination not in known_sizes or known_sizes[destination] == size,
                    f"{where}: repeated state has inconsistent derived size")
            known_sizes[destination] = size
            if first_above_ten is None and size > 10:
                first_above_ten = {"step_0_based": index, "size": size,
                                   "signature": destination}
            sizes[str(size)] += 1
            states.add(destination)
            previous = destination
        require(previous == end, f"{label}: wrong final state")
        require(size == 6, f"{label}: derived endpoint size differs from six")
        derived[label] = {
            "start": start, "end": end, "moves": len(steps),
            "move_counts": dict(sorted(counts.items())),
            "move_face_degree_counts": dict(sorted(degrees.items())),
            "state_visits_including_endpoints": len(steps) + 1,
            "unique_states": len(states),
            "pentachora_min": min(map(int, sizes)),
            "pentachora_max": max(map(int, sizes)),
            "state_visit_size_histogram": dict(sorted(sizes.items(), key=lambda p: int(p[0]))),
            "first_state_above_ten": first_above_ten,
            "two_zero_counts": {code: n for code, n in sorted(counts.items())
                                if code.startswith("Z")},
            "largest_two_zero_source_size": max(zero_sizes, default=0),
        }
        all_states.update(states)
        total_moves.update(counts)
    require(derived["G6_to_X"]["end"] == derived["std6_to_X"]["end"]
            and derived["std6_to_X"]["start"] == derived["sphere_to_std6"]["end"],
            "terminal paths do not meet exactly")
    return derived, all_states, dict(sorted(total_moves.items()))


def _check_isomorphism(source: Any, destination: Any) -> dict[str, Any]:
    iso = source.isIsomorphicTo(destination)
    require(iso is not None, "capped pairing has no direct isomorphism to terminal start")
    require(source.size() == destination.size(), "isomorphism sizes differ")
    images = [iso.simpImage(i) for i in range(source.size())]
    permutations = [[iso.facetPerm(i)[k] for k in range(5)]
                    for i in range(source.size())]
    require(sorted(images) == list(range(destination.size())),
            "isomorphism simplex map is not bijective")
    require(all(sorted(p) == list(range(5)) for p in permutations),
            "isomorphism has invalid vertex permutations")
    equations = 0
    for i, permutation in enumerate(permutations):
        first, second = source.pentachoron(i), destination.pentachoron(images[i])
        for facet in range(5):
            a, b = first.adjacentSimplex(facet), second.adjacentSimplex(permutation[facet])
            require(a is not None and b is not None, "isomorphism has unexpected boundary")
            require(b.index() == images[a.index()], "isomorphism adjacency mismatch")
            p, q = first.adjacentGluing(facet), second.adjacentGluing(permutation[facet])
            require(all(permutations[a.index()][p[k]] == q[permutation[k]] for k in range(5)),
                    "isomorphism gluing conjugacy equation failed")
            equations += 1
    require(equations == 30, "wrong number of capped isomorphism equations")
    return {"simplex_images": images, "vertex_permutations": permutations,
            "paired_half_facet_equations": equations}


def _check_capped(path: Path, canonical: Any) -> dict[str, Any]:
    path, raw, capped = _json_file(path, "capped endpoint")
    require(set(capped) == {"simplex_count", "gluings", "isoSig", "f_vector",
                            "numbered_sha256", "numbered_rows"},
            "capped endpoint has unexpected schema")
    require(type(capped["simplex_count"]) is int and capped["simplex_count"] == 6,
            "capped endpoint must contain six pentachora")
    require(capped["isoSig"] == G6, "capped endpoint metadata differs from terminal start")
    rows = capped["numbered_rows"]
    require(type(rows) is list and len(rows) == 6
            and all(type(row) is list and len(row) == 5 for row in rows),
            "capped endpoint numbered-pairing shape is invalid")
    for i, row in enumerate(rows):
        for facet, entry in enumerate(row):
            require(type(entry) is list and len(entry) == 2,
                    "capped endpoint has malformed paired half-facet")
            neighbor, encoded = entry
            require(type(neighbor) is int and 0 <= neighbor < 6
                    and type(encoded) is str and sorted(encoded) == list("01234"),
                    "capped endpoint has invalid neighbor or permutation")
            require((i, facet) != (neighbor, int(encoded[facet])),
                    "a capped half-facet is paired to itself")
    tri = regina.Triangulation4()
    tri.newPentachora(6)
    equations = 0
    for i, row in enumerate(rows):
        for facet, (neighbor, encoded) in enumerate(row):
            permutation = tuple(map(int, encoded))
            back, reverse_encoded = rows[neighbor][permutation[facet]]
            reverse = tuple(map(int, reverse_encoded))
            require(back == i and all(reverse[permutation[k]] == k for k in range(5)),
                    "capped endpoint gluing is not reciprocal")
            if tri.pentachoron(i).adjacentSimplex(facet) is None:
                tri.pentachoron(i).join(facet, tri.pentachoron(neighbor),
                                       regina.Perm5(*permutation))
            equations += 1
    require(equations == 30, "wrong number of reciprocal capped equations")
    actual_rows = [[[p.adjacentSimplex(f).index(), str(p.adjacentGluing(f))]
                    for f in range(5)] for p in tri.pentachora()]
    require(actual_rows == rows, "constructed capped numbering differs from supplied rows")
    require(capped["numbered_sha256"] == _sha(json.dumps(actual_rows).encode()),
            "capped numbered-pairing hash mismatch")
    expected_gluings = []
    for i, row in enumerate(actual_rows):
        for facet, (neighbor, encoded) in enumerate(row):
            permutation = list(map(int, encoded))
            if (i, facet) < (neighbor, permutation[facet]):
                expected_gluings.append([i, facet, neighbor, permutation])
    require(capped["gluings"] == expected_gluings and len(expected_gluings) == 15,
            "capped redundant gluing list disagrees with its numbered pairing")
    require(tri.isoSig() == G6, "actual capped pairing differs from terminal start")
    properties = _properties(tri)
    require(capped["f_vector"] == properties["f_vector"], "capped f-vector mismatch")
    isomorphism = _check_isomorphism(tri, canonical)
    return {"path": str(path), "sha256": _sha(raw), "isoSig": tri.isoSig(),
            "properties": properties, "reciprocal_half_facet_equations": equations,
            "redundant_gluing_encodings_agree": True,
            "numbered_sha256": capped["numbered_sha256"],
            "direct_isomorphism_to_terminal_start": isomorphism}


def verify_terminal_bindings(
    certificate: Path,
    capped_path: Path | None,
    independent_receipt: Path,
    two_zero_receipt: Path,
) -> dict[str, Any]:
    """Bind supplied replay receipts and optional actual capped endpoint.

    The caller must produce the two receipts by fresh subprocess executions.
    This function never reads an old acceptance record or default receipt.
    """
    certificate, raw, data = _json_file(certificate, "terminal certificate")
    digest = _sha(raw)
    require(digest == CERTIFICATE_SHA256, "terminal certificate bytes differ from reviewed input")
    independent_receipt, replay_raw, replay = _json_file(independent_receipt, "terminal replay receipt")
    two_zero_receipt, zero_raw, zero = _json_file(two_zero_receipt, "two-zero receipt")
    require(certificate != independent_receipt and certificate != two_zero_receipt
            and independent_receipt != two_zero_receipt, "terminal input paths must be distinct")
    require(replay.get("status") == "PASS" and zero.get("status") == "PASS",
            "one or both terminal replay receipts did not pass")
    require(replay.get("schema") == "x0-standardness-independent-replay-v1",
            "unexpected full replay receipt schema")
    for label, receipt in (("full replay", replay), ("two-zero", zero)):
        require(receipt.get("certificate_sha256") == digest,
                f"{label} receipt is bound to different certificate bytes")
    source_root = Path(__file__).resolve().parent / "source" / "terminal"
    replay_source = _input_path(source_root / "independent_replay.py", "independent replay source")
    zero_source = _input_path(source_root / "independent_two_zero.py", "two-zero source")
    require(_sha(replay_source.read_bytes()) == replay.get("verifier_sha256") == REPLAY_SHA256,
            "full replay receipt/checker identity mismatch")
    require(_sha(zero_source.read_bytes()) == zero.get("checker_sha256") == TWO_ZERO_SHA256,
            "two-zero receipt/checker identity mismatch")
    package_version = importlib.metadata.version("regina")
    require(package_version == "7.4.1" and replay.get("regina_package_version") == package_version,
            "terminal checks require the recorded Regina package 7.4.1")
    require(replay.get("regina_engine_version") == zero.get("regina_version") == regina.versionString(),
            "terminal receipts and current Regina engine versions disagree")
    scope = _mapping(replay.get("scope"), "full replay scope")
    for flag in ("independent_of_supplied_checker", "canonical_states",
                 "all_unique_states_checked_valid_closed_connected_orientable_chi2",
                 "all_moves_have_separate_has_guard_and_checked_execution",
                 "independent_basic_face_degree_and_distinct_incident_pentachora_checks"):
        require(scope.get(flag) is True, f"full replay omitted required check: {flag}")
    require(zero.get("regina_move_functions_called") is False,
            "second two-zero receipt does not attest separate move implementation")
    require(zero.get("guard_method") == "closed-input doubled-link check plus forests of identified faces"
            and zero.get("output_method") == "new facet pairing by explicit traversal through flattened pillow",
            "unexpected second two-zero checking method")
    derived, all_states, total_counts = _derive_paths(data)
    replay_paths = _mapping(replay.get("paths"), "replay paths")
    zero_by_path = _mapping(zero.get("by_path"), "two-zero path counts")
    require(set(replay_paths) == set(zero_by_path) == set(PATHS),
            "receipt path set differs from the certificate")
    max_vertices = 0
    for label, expected in derived.items():
        actual = _mapping(replay_paths[label], f"{label} replay result")
        require(actual.get("status") == "PASS", f"{label}: replay did not pass")
        for key, value in expected.items():
            if key not in ("two_zero_counts", "largest_two_zero_source_size"):
                require(actual.get(key) == value, f"{label}: derived {key} differs from replay receipt")
        require(_counts(zero_by_path[label], f"{label} two-zero counts") == expected["two_zero_counts"],
                f"{label}: second-checker move counts differ from certificate")
        vertex_histogram = _counts(actual.get("state_visit_vertex_histogram"), f"{label} vertex histogram")
        require(all(k.isdecimal() and str(int(k)) == k and int(k) > 0 for k in vertex_histogram),
                f"{label}: invalid vertex histogram labels")
        require(sum(vertex_histogram.values()) == expected["moves"] + 1,
                f"{label}: inconsistent vertex-histogram total")
        require(actual.get("vertices_min") == min(map(int, vertex_histogram))
                and actual.get("vertices_max") == max(map(int, vertex_histogram)),
                f"{label}: inconsistent vertex range")
        max_vertices = max(max_vertices, actual["vertices_max"])
    all_moves = sum(total_counts.values())
    two_zero_counts = {code: count for code, count in total_counts.items() if code.startswith("Z")}
    require(replay.get("all_moves") == all_moves
            and _counts(replay.get("all_move_counts"), "aggregate replay move counts") == total_counts,
            "aggregate replay move counts differ from the certificate")
    require(replay.get("distinct_states_all_paths") == len(all_states),
            "distinct replay state count differs from the certificate")
    require(replay.get("global_max_pentachora") == max(p["pentachora_max"] for p in derived.values())
            and replay.get("global_max_vertices") == max_vertices,
            "global replay size range is inconsistent")
    require(_counts(zero.get("counts"), "aggregate two-zero counts") == two_zero_counts
            and zero.get("two_zero_moves_checked") == sum(two_zero_counts.values()),
            "aggregate second-checker counts differ from the certificate")
    require(zero.get("largest_source_size") == max(p["largest_two_zero_source_size"] for p in derived.values()),
            "second-checker largest source size differs from the certificate")
    sphere = _standard_double()
    anchors = {SPHERE2: sphere}
    for signature in (G6, MEETING, STANDARD6):
        anchors[signature] = regina.Triangulation4.fromIsoSig(signature)
        require(anchors[signature].isoSig() == signature, "noncanonical terminal endpoint")
    anchor_properties = {signature: _properties(tri) for signature, tri in anchors.items()}
    require(replay.get("anchor_properties") == anchor_properties,
            "replay endpoint properties differ from independently reconstructed anchors")
    standard_anchor = _mapping(replay.get("standard_anchor"), "standard anchor receipt")
    require(standard_anchor.get("isoSig") == SPHERE2
            and standard_anchor.get("matches_Example4_sphere") is True
            and standard_anchor.get("properties") == anchor_properties[SPHERE2],
            "replay standard anchor differs from the explicit simplex double")
    cap = _check_capped(capped_path, anchors[G6]) if capped_path is not None else None
    return {
        "schema": "x0-reduction-terminal-binding-v1", "status": "PASS",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "checker_sha256": _sha(Path(__file__).read_bytes()),
        "regina_package_version": package_version, "regina_engine_version": regina.versionString(),
        "certificate": {"path": str(certificate), "sha256": digest},
        "receipts": {
            "independent_replay": {"path": str(independent_receipt), "sha256": _sha(replay_raw),
                                   "checker_sha256": REPLAY_SHA256},
            "independent_two_zero": {"path": str(two_zero_receipt), "sha256": _sha(zero_raw),
                                     "checker_sha256": TWO_ZERO_SHA256},
        },
        "paths": derived, "all_moves": all_moves, "all_move_counts": total_counts,
        "pachner_moves": sum(n for code, n in total_counts.items() if code.startswith("P")),
        "two_zero_moves": sum(two_zero_counts.values()), "distinct_states": len(all_states),
        "exact_path_intersections_checked": True,
        "standard_anchor": {"isoSig": SPHERE2, "properties": anchor_properties[SPHERE2],
                            "construction": "Double of one abstract four-simplex: two distinct pentachora with corresponding facets joined by identity"},
        "capped_endpoint": cap,
        "capped_endpoint_bound": cap is not None,
        "scope": {
            "role": "Composition and input/code/receipt binding after separate move replay",
            "move_replay_performed_by_this_helper": False,
            "receipt_freshness": "Caller responsibility; this helper binds the supplied receipts but does not establish when their move replays occurred",
            "whole_terminal_certificate_pachner_only": False,
            "geometric_construction_identification_checked_here": False,
            "earlier_reduction_checked_here": False,
            "remaining_native_dependencies": ["Regina signature codec", "Regina skeleton and validity", "Regina simplex construction and combinatorial isomorphism"],
        },
        "conclusion": "The checked terminal paths identify their six-pentachoron start with the explicit double of a four-simplex by Pachner and admissible two-zero moves. The supplied capped endpoint is included in this conclusion only when capped_endpoint_bound is true.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--capped", type=Path)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--two-zero", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="new JSON output path")
    args = parser.parse_args(argv)
    output = args.output.expanduser().resolve()
    if output.exists():
        print(f"FAIL: refusing to overwrite an existing output: {output}", file=sys.stderr)
        return 2
    inputs = [args.certificate, args.replay, args.two_zero, Path(__file__)]
    if args.capped is not None:
        inputs.append(args.capped)
    if any(output == p.expanduser().resolve() for p in inputs):
        print("FAIL: output path aliases an input or this checker", file=sys.stderr)
        return 2
    try:
        result = verify_terminal_bindings(args.certificate, args.capped, args.replay, args.two_zero)
        status = 0
    except Exception as error:
        result = {"schema": "x0-reduction-terminal-binding-v1", "status": "FAIL",
                  "checked_utc": datetime.now(timezone.utc).isoformat(),
                  "error_type": type(error).__name__, "error": str(error)}
        status = 1
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "receipt": str(output),
                      "error": result.get("error")}))
    return status


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Independent, assertion-free replay of the supplied X0 standardness paths.

This script imports no reviewer/search module. It independently checks the
certificate structure, canonical state chain, face indices, degree and
incident-simplex requirements, Regina's modern guarded APIs, move size changes,
and exact resulting signatures. Regina remains the shared move/isoSig engine.

Default input is the unchanged extracted reviewer certificate. Run with
Python -O if desired: essential verification does not depend on assert.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time

import regina


EXPECTED_CERTIFICATE_SHA256 = "e9ba9f14831bb0da32043d7d93eb0148f0683e2a6e4f67d5156e8e9fc94adcdf"
G6 = "gAMLMQaabcdeffff2a+a2avbmaoaxaOb0aUa"
SPHERE2 = "cPkbbbbaaaaaaaa"
STD6 = "gLAPPQacbedeeffffaaawaaaobcacaQbPboa"
MEETING = "gLAMPQaccddefefffaaaubaaMafapboapboa"
FACE_METHOD = {0: "vertex", 1: "edge", 2: "triangle", 3: "tetrahedron", 4: "pentachoron"}
COUNT_METHOD = {0: "countVertices", 1: "countEdges", 2: "countTriangles", 3: "countTetrahedra", 4: "size"}
EXPECTED_PATHS = {
    "G6_to_X": (G6, MEETING, 6814, {"P2", "P3", "Z1", "Z2"}),
    "std6_to_X": (STD6, MEETING, 1519, {"P2", "P3", "Z1", "Z2"}),
    "sphere_to_std6": (SPHERE2, STD6, 103, {"P0", "P1", "P2", "P3", "P4"}),
}


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def load_json_no_duplicate_keys(data):
    def pairs_hook(pairs):
        out = {}
        for key, value in pairs:
            require(key not in out, "duplicate JSON key: " + key)
            out[key] = value
        return out
    return json.loads(data, object_pairs_hook=pairs_hook)


def construct_standard_double():
    """Glue each facet of two disjoint abstract 4-simplices by the identity."""
    t = regina.Triangulation4()
    a, b = t.newPentachora(2)
    for facet in range(5):
        a.join(facet, b, regina.Perm5())
    require(t.isoSig() == SPHERE2, "explicit simplex double has wrong signature")
    require(regina.Example4.sphere().isoSig() == SPHERE2,
            "Regina sphere example disagrees with explicit simplex double")
    return t


def check_state(t, signature, state_cache):
    if signature in state_cache:
        return state_cache[signature]
    require(isinstance(signature, str) and signature, "empty/non-string signature")
    require(t.isoSig() == signature, "state is not a canonical classic isoSig: " + signature)
    require(t.size() > 0, "empty state")
    require(t.isConnected(), "disconnected state: " + signature)
    require(t.isValid(), "invalid state: " + signature)
    require(t.isClosed(), "nonclosed or ideal state: " + signature)
    require(t.isOrientable(), "nonorientable state: " + signature)
    require(t.eulerCharTri() == 2, "Euler characteristic is not two: " + signature)
    stat = {
        "pentachora": t.size(),
        "vertices": t.countVertices(),
        "f_vector": [t.countVertices(), t.countEdges(), t.countTriangles(),
                     t.countTetrahedra(), t.size()],
        "vertex_degrees": sorted(v.degree() for v in t.vertices()),
    }
    state_cache[signature] = stat
    return stat


def verify_path(label, steps, specification, state_cache):
    start, end, length, allowed = specification
    require(isinstance(steps, list) and len(steps) == length,
            label + ": unexpected path length/type")
    expected_source = start
    move_counts = Counter()
    state_sigs = {start}
    sizes = Counter()
    vertex_counts = Counter()
    degrees = Counter()
    first_exceed_ten = None
    start_time = time.perf_counter()
    for index, step in enumerate(steps):
        where = f"{label} step {index}"
        require(isinstance(step, list) and len(step) == 4, where + ": invalid record shape")
        source, code, face_index, destination = step
        require(isinstance(source, str) and isinstance(destination, str),
                where + ": state is not a string")
        require(source == expected_source, where + ": chain discontinuity or wrong start")
        require(isinstance(code, str) and code in allowed, where + ": forbidden move code")
        require(type(face_index) is int and face_index >= 0, where + ": invalid face index")
        t = regina.Triangulation4.fromIsoSig(source)
        source_info = check_state(t, source, state_cache)
        dimension = int(code[1])
        require(face_index < getattr(t, COUNT_METHOD[dimension])(), where + ": face out of range")
        face = getattr(t, FACE_METHOD[dimension])(face_index)
        expected_degree = 5 - dimension if code.startswith("P") else 2
        if dimension < 4:
            require(face.isValid(), where + ": invalid move face")
            require(not face.isBoundary(), where + ": boundary move face")
            require(face.degree() == expected_degree, where + ": wrong move-face degree")
            incident = [embedding.simplex().index() for embedding in face.embeddings()]
            require(len(incident) == expected_degree and len(set(incident)) == expected_degree,
                    where + ": incident pentachora are not distinct")
        degrees[f"{code}:{expected_degree}"] += 1
        size_before = t.size()
        if code.startswith("P"):
            require(t.hasPachner(face), where + ": Pachner legality guard rejected")
            require(t.pachner(face), where + ": guarded Pachner move failed")
            expected_delta = 2 * dimension - 4
        else:
            require(t.has20(face), where + ": 2-0 legality guard rejected")
            require(t.move20(face), where + ": guarded 2-0 move failed")
            expected_delta = -2
        require(t.size() == size_before + expected_delta, where + ": incorrect size change")
        require(t.isoSig() == destination, where + ": exact target signature mismatch")
        destination_info = check_state(t, destination, state_cache)
        sizes[source_info["pentachora"]] += 1
        vertex_counts[source_info["vertices"]] += 1
        if first_exceed_ten is None and destination_info["pentachora"] > 10:
            first_exceed_ten = {"step_0_based": index, "size": destination_info["pentachora"],
                                "signature": destination}
        move_counts[code] += 1
        state_sigs.add(destination)
        expected_source = destination
    require(expected_source == end, label + ": wrong final endpoint")
    final_info = state_cache[end]
    sizes[final_info["pentachora"]] += 1
    vertex_counts[final_info["vertices"]] += 1
    result = {
        "status": "PASS", "start": start, "end": end,
        "moves": length, "move_counts": dict(sorted(move_counts.items())),
        "state_visits_including_endpoints": length + 1,
        "unique_states": len(state_sigs),
        "pentachora_min": min(sizes), "pentachora_max": max(sizes),
        "vertices_min": min(vertex_counts), "vertices_max": max(vertex_counts),
        "state_visit_size_histogram": dict(sorted(sizes.items())),
        "state_visit_vertex_histogram": dict(sorted(vertex_counts.items())),
        "move_face_degree_counts": dict(sorted(degrees.items())),
        "first_state_above_ten": first_exceed_ten,
        "elapsed_seconds": round(time.perf_counter() - start_time, 6),
    }
    print(label + ": PASS " + str(length) + " guarded moves", flush=True)
    return result


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path,
                        default=root / "source/x0_is_standard/x0_standard_certificate.json")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("independent_replay_receipt.json"))
    args = parser.parse_args()
    require(not args.output.exists(), "refusing to overwrite an existing audit receipt")
    receipt = {
        "schema": "x0-standardness-independent-replay-v1",
        "status": "RUNNING",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(args.certificate.resolve()),
        "verifier_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "python": platform.python_version(),
        "python_optimization_level": sys.flags.optimize,
        "regina_package_version": importlib.metadata.version("regina"),
        "regina_engine_version": regina.versionString(),
        "scope": {
            "independent_of_supplied_checker": True,
            "shared_engine": "Regina guarded Pachner/2-0 operations and isomorphism signatures",
            "canonical_states": True,
            "all_unique_states_checked_valid_closed_connected_orientable_chi2": True,
            "all_moves_have_separate_has_guard_and_checked_execution": True,
            "independent_basic_face_degree_and_distinct_incident_pentachora_checks": True,
            "full_non_Regina_reimplementation_of_move_guards": False,
            "new_explicit_sphere_recognition_of_all_vertex_links": False,
            "accepted_X0_to_G6_certificate_replayed": False,
            "random_walk_exploration_counts_verified": False,
        },
    }
    t0 = time.perf_counter()
    try:
        raw = args.certificate.read_bytes()
        receipt["certificate_sha256"] = sha256_bytes(raw)
        require(receipt["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256,
                "input bytes do not match the reviewed certificate")
        c = load_json_no_duplicate_keys(raw)
        require(type(c) is dict and set(c) == set(EXPECTED_PATHS) | {"X", "std6"},
                "unexpected certificate top-level keys")
        require(c["X"] == MEETING and c["std6"] == STD6, "endpoint metadata differs")
        state_cache = {}
        double = construct_standard_double()
        check_state(double, SPHERE2, state_cache)
        receipt["standard_anchor"] = {
            "construction": "Two distinct 4-simplices; each of five corresponding facets glued by the identity permutation",
            "isoSig": SPHERE2,
            "matches_Example4_sphere": True,
            "properties": state_cache[SPHERE2],
        }
        receipt["paths"] = {
            label: verify_path(label, c[label], spec, state_cache)
            for label, spec in EXPECTED_PATHS.items()
        }
        receipt["all_moves"] = sum(x["moves"] for x in receipt["paths"].values())
        total_counts = Counter()
        for p in receipt["paths"].values():
            total_counts.update(p["move_counts"])
        receipt["all_move_counts"] = dict(sorted(total_counts.items()))
        receipt["distinct_states_all_paths"] = len(state_cache)
        receipt["global_max_pentachora"] = max(x["pentachora"] for x in state_cache.values())
        receipt["global_max_vertices"] = max(x["vertices"] for x in state_cache.values())
        receipt["anchor_properties"] = {s: state_cache[s] for s in (G6, MEETING, STD6, SPHERE2)}
        receipt["conclusion"] = (
            "The three supplied paths connect G6 and the explicitly constructed double of a "
            "4-simplex by checked Regina Pachner and 2-0 PL-preserving moves, joined at exactly "
            "matching canonical states. This independently verifies the new G6 standardness "
            "certificate at the Regina move/isoSig trust level. X0 standardness additionally "
            "inherits the previously accepted identification and X0-to-G6 certificate."
        )
        receipt["status"] = "PASS"
    except Exception as error:
        receipt["status"] = "FAIL"
        receipt["error_type"] = type(error).__name__
        receipt["error"] = str(error)
    receipt["elapsed_seconds"] = round(time.perf_counter() - t0, 6)
    receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": receipt["status"], "elapsed_seconds": receipt["elapsed_seconds"],
                      "receipt": str(args.output), "error": receipt.get("error")}), flush=True)
    return 0 if receipt["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

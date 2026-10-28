#!/usr/bin/env python3
"""Small subprocess adapters for two construction audit stages.

These adapters use the archived independent checker for source coarsening and
the archived native source producer followed by the independent endpoint
checker. They do not implement a new simplification search.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
import time
from itertools import combinations
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def source_coarsening(work, receipt):
    """Replay all source replacements and reconstruct the exact descendant."""
    started = time.monotonic()
    root = work / "frozen_release"
    paths = {"master": root / "coarsened/master_source.json.gz",
             "certificate": root / "coarsened/coarsening_certificate.json.gz",
             "descendant": root / "coarsened/source.json.gz",
             "checker": root / "lib/verify_compact3.py"}
    master, certificate, descendant = (load(paths[k]) for k in
                                       ("master", "certificate", "descendant"))
    checker = module(paths["checker"], "archived_independent_source_coarsening")
    deck = master["deck"]
    original = [tuple(sorted(tet, key=lambda v: min(v, deck[v])))
                for tet in master["facets"]]
    tags = sorted({v for tet in original for v in tet}, key=repr)
    deck_map = {v: deck[v] for v in tags}
    protected = {v for v in tags if deck[v] == v}
    endpoint, checked = checker.replay_part(certificate, original,
                                           deck_map, protected)
    require(len(certificate["steps"]) == 588, "Wrong source coarsening trace")
    active = sorted({v for tet in endpoint for v in tet})
    indices = {v: i for i, v in enumerate(active)}
    facets = sorted(tuple(sorted(indices[v] for v in tet)) for tet in endpoint)
    endpoint_deck = [indices[deck[v]] for v in active]
    rebuilt = {
        "schema": "marked-pretzel-source-descendant-v1",
        "facets": facets, "closed_cover_facets": facets,
        "deck": endpoint_deck,
        "vertex_ranks": [min(i, endpoint_deck[i]) for i in range(len(active))],
        "parent_vertices": active,
        "parent_source_sha256": digest(master),
        "parent_file_sha256": sha(paths["master"]),
        "ancestry": "equivariant pairs of independently replayed relative 3-ball replacements; fixed circle pointwise preserved",
        "ramification_edges": sorted({edge for tet in facets
            for edge in combinations(tet, 2)
            if all(endpoint_deck[v] == v for v in edge)}),
        "source_interpretation": master["source_interpretation"],
    }
    require(digest(rebuilt) == digest(descendant),
            "Replayed source descendant differs from the complete stored object")
    require(len(facets) == 628, "Wrong source descendant size")
    result = {"schema": "fresh-source-coarsening-audit-v1", "status": "PASS",
              "source_tetrahedra": len(original),
              "descendant_tetrahedra": len(facets),
              "complete_descendant_canonical_sha256": digest(rebuilt),
              "all_descendant_fields_reconstructed": True,
              "checked": checked,
              "input_file_sha256": {k: sha(p) for k, p in paths.items()},
              "producer_search_rerun": False,
              "assertions_enabled": not sys.flags.optimize,
              "wall_seconds": time.monotonic() - started}
    receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


def native_transport(work, receipt):
    """Fresh native peripheral transport, followed by exact endpoint checks."""
    started = time.monotonic()
    root = work / "frozen_release"
    # The archived producer imports its sibling prepare_source module, which
    # in turn imports the small lib/ helpers. Import only from this stage copy.
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(root / "lib"))
    expected_coned = load(root / "peripheral/marked_coned_exterior.json.gz")
    producer = module(root / "marked_exterior.py", "archived_marked_exterior")
    transport = producer.run()
    require(digest(load(root / "peripheral/marked_coned_exterior.json.gz")) ==
            digest(expected_coned),
            "Native transport regenerated a different marked compact source")
    checker = module(root / "diagram/verify_marked_isomorphism.py",
                     "fresh_independent_marked_endpoint")
    exact = checker.build_certificate(root / "peripheral")
    selected = exact["all_isomorphisms"][
        exact["selected_identity_peripheral_isomorphism"]]
    require(exact["verified"] is True, "Fresh marked endpoint check failed")
    require(selected["orientation_sign"] == 1,
            "Fresh endpoint reverses the specified ambient orientation")
    require(selected["cusp_matrix_columns_images_of_source_M_L"] == [[1, 0], [0, 1]],
            "Fresh endpoint changes a specified peripheral class")
    endpoint_path = receipt.with_name("native_transport_endpoint.json")
    endpoint_path.write_text(json.dumps(exact, indent=2, sort_keys=True) + "\n")
    result = {"schema": "fresh-native-source-transport-audit-v1",
              "status": "PASS",
              "compact_source_canonical_sha256": digest(expected_coned),
              "native_transport": transport,
              "selected_exact_endpoint": selected,
              "endpoint_certificate": endpoint_path.name,
              "endpoint_certificate_sha256": sha(endpoint_path),
              "scope": "Fresh native import, finite-vertex removal and peripheral transport, followed by independent exact marked endpoint verification.",
              "trust_boundary": "SnapPy/Spherogram native constructor and move/peripheral-transport semantics remain trusted; this is not an elementary 3D move transcript.",
              "assertions_enabled": not sys.flags.optimize,
              "wall_seconds": time.monotonic() - started}
    receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("source-coarsening", "native-transport"))
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    require(not sys.flags.optimize, "Run with Python assertions enabled")
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    action = source_coarsening if args.stage == "source-coarsening" else native_transport
    action(args.work.resolve(), args.receipt.resolve())


if __name__ == "__main__":
    main()

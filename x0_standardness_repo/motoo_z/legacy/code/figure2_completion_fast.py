#!/usr/bin/env python3
"""Explicit small truncation for the colour-preserving Katie triangulation.

When every pentachoron has exactly one corner at the chosen vertex, and
the remaining corners have a face-compatible total order, removing a
small open neighbourhood of that vertex leaves a tetrahedral prism in
each pentachoron.  The standard staircase triangulates each such prism
with FOUR pentachora.  This is an explicit PL truncation, not a
recognition heuristic and not a singular cone filling.

Here the raw Katie triangulation has identity gluing permutations, so
local vertex colours supply the compatible ordering.  The implementation
checks this precondition on every original gluing.  All local labels and
all gluing equations are reproducible from the input and the algorithm.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import regina


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def save4(t, path):
    require(regina.PacketOfTriangulation4(t).save(str(path)), "save failed")


def save3(t, path):
    require(regina.PacketOfTriangulation3(t).save(str(path)), "save failed")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def staircase_truncate(source, vertex_index):
    """Truncate one vertex using a verified globally compatible ordering.

    A symbolic vertex (colour,0) is the point on the radial edge just
    outside the deleted neighbourhood.  (colour,1) is the old finite
    vertex. The staircase simplex with split i contains low vertices
    through i and high vertices from i onwards.
    """
    count = source.size()
    apex = []
    orders = []
    for old in source.simplices():
        candidates = [v for v in range(5) if old.vertex(v).index() == vertex_index]
        require(len(candidates) == 1, "every old simplex must contain exactly one apex corner")
        a = candidates[0]
        apex.append(a)
        require(len({old.vertex(v).index() for v in range(5)}) == 5,
                "each old simplex must have five distinct global vertices")
        # Global vertex indices supply a face-compatible ordering even
        # after guarded moves have destroyed Katie's local colouring.
        orders.append(sorted((v for v in range(5) if v != a),
                             key=lambda v: old.vertex(v).index()))
    for old in source.simplices():
        i = old.index()
        for facet in range(5):
            other = old.adjacentSimplex(facet)
            if other is None:
                continue
            j = other.index()
            gluing = old.adjacentGluing(facet)
            mapped = [gluing[v] for v in orders[i] if v != facet]
            wanted = [v for v in orders[j] if v != gluing[facet]]
            require(mapped == wanted, "face pairing does not preserve the staircase order")
            if facet != apex[i]:
                require(gluing[apex[i]] == apex[j], "face gluing does not preserve apex")

    result = regina.Triangulation4()
    result.newSimplices(4 * count)
    labels = []
    outer = {}
    radial_boundary = []
    internal_pairs = 0
    old_pairs = 0

    for old in source.simplices():
        i = old.index()
        order = orders[i]
        unmatched = {}
        for split in range(4):
            row = [(v, 0) for v in order[:split+1]] + [(v, 1) for v in order[split:]]
            require(len(row) == len(set(row)) == 5, "degenerate staircase simplex")
            simplex = 4*i + split
            labels.append(row)
            for facet in range(5):
                face = tuple(sorted(row[:facet] + row[facet+1:]))
                if face in unmatched:
                    other_simplex, other_facet = unmatched.pop(face)
                    other_labels = labels[other_simplex]
                    permutation = [other_facet if v == facet else other_labels.index(row[v])
                                   for v in range(5)]
                    result.simplex(simplex).join(facet, result.simplex(other_simplex),
                                                regina.Perm5(*permutation))
                    internal_pairs += 1
                else:
                    unmatched[face] = (simplex, facet)
        for face, location in unmatched.items():
            colours = {v for v, level in face}
            levels = {level for v, level in face}
            if levels == {0}:
                require(len(colours) == 4, "bad new truncation boundary")
                radial_boundary.append({"old_pentachoron": i,
                                        "new_pentachoron": location[0],
                                        "new_facet": location[1]})
                continue
            if len(colours) == 4:
                require(levels == {1}, "unrecognised old opposite facet")
                old_facet = apex[i]
            else:
                require(len(colours) == 3, "unrecognised prism side")
                old_facet = next(v for v in order if v not in colours)
            outer[(i, old_facet, face)] = location

    visited = set()
    for key, location in outer.items():
        i, old_facet, face = key
        if key in visited:
            continue
        old = source.simplex(i)
        other = old.adjacentSimplex(old_facet)
        if other is None:
            continue
        gluing = old.adjacentGluing(old_facet)
        mapped = tuple(sorted((gluing[v], level) for v, level in face))
        other_key = (other.index(), gluing[old_facet], mapped)
        require(other_key in outer, "old facet subdivision has no matching mate")
        other_location = outer[other_key]
        a, af = location
        b, bf = other_location
        permutation = []
        for v, (colour, level) in enumerate(labels[a]):
            permutation.append(bf if v == af else labels[b].index((gluing[colour], level)))
        require(sorted(permutation) == list(range(5)), "invalid prism gluing")
        result.simplex(a).join(af, result.simplex(b), regina.Perm5(*permutation))
        visited.update((key, other_key))
        old_pairs += 1
    require(len(radial_boundary) == count, "wrong truncation boundary count")
    require(result.size() == 4 * count, "wrong prism simplex count")
    return result, {
        "construction": "staircase triangulation of each truncated pentachoron as Delta3 x I",
        "chosen_vertex": vertex_index,
        "input_pentachora": count,
        "output_pentachora": result.size(),
        "new_simplices_per_old": 4,
        "local_apex_colours": apex,
        "local_remaining_colours_in_order": orders,
        "new_simplex_symbolic_vertices": labels,
        "internal_prism_facet_pairs": internal_pairs,
        "facet_pairs_across_old_simplices": old_pairs,
        "truncation_boundary_locations": radial_boundary,
        "all_old_gluings_preserve_local_order": True,
    }


def audit_source(source, out):
    records = []
    ideal_candidates = []
    for vertex in source.vertices():
        link = vertex.buildLink()
        initial = link.size()
        then = time.monotonic()
        link.simplify()
        sphere = link.isSphere()
        record = {"vertex": vertex.index(), "original_link_tetrahedra": initial,
                  "simplified_link_tetrahedra": link.size(), "simplified_link_isosig": link.isoSig(),
                  "simplified_link_is_sphere": bool(sphere), "H1": str(link.homology()),
                  "seconds": time.monotonic() - then,
                  "trust_boundary": "native Regina3D simplify and sphere recognition on a link copy"}
        if not sphere:
            summands = link.summands()
            ref = regina.Example3.s2xs1()
            expected = len(summands) == 2 and all(s.isIsomorphicTo(ref) is not None for s in summands)
            record["summand_isosigs"] = [s.isoSig() for s in summands]
            record["recognised_two_S1xS2_summands"] = bool(expected)
            require(expected, "non-spherical link is not independently recognised as #2(S1xS2)")
            ideal_candidates.append(vertex.index())
        records.append(record)
        print(json.dumps({"link_audit": record}), flush=True)
    require(len(ideal_candidates) == 1, "expected one ideal vertex")
    report = {"status": "PASS_NATIVE_LINK_RECOGNITION", "vertices": records,
              "ideal_vertex": ideal_candidates[0],
              "all_other_links_are_3_spheres": True,
              "ideal_link_is_two_S1xS2": True}
    save_json(out / "SOURCE_VERTEX_LINKS.json", report)
    return ideal_candidates[0]


def construct(source_path, out_path, vertex_index=None):
    source_path, out = Path(source_path), Path(out_path)
    require(not out.exists(), "output already exists")
    out.mkdir(parents=True)
    source = regina.Triangulation4(regina.open(str(source_path)))
    original = regina.Triangulation4(source)
    started = time.monotonic()
    require(source.isValid() and source.isConnected() and source.isOrientable(), "bad raw input")
    recognised_vertex = audit_source(source, out)
    if vertex_index is None:
        vertex_index = recognised_vertex
    require(vertex_index == recognised_vertex, "requested vertex differs from independently recognised ideal vertex")
    result, proof = staircase_truncate(source, vertex_index)
    save4(result, out / "finite_staircase.rga")
    save_json(out / "STAIRCASE_CONSTRUCTION.json", proof)
    require(result.isValid() and result.isConnected() and result.isOrientable(), "bad staircase output")
    require(source == original, "source was modified")
    # Every original vertex link was certified above. The explicit local
    # prism construction establishes that the new piece is finite without
    # requiring expensive fresh 4D isIdeal() recognition on every move.
    real = [b for b in result.boundaryComponents() if b.isReal()]
    require(len(real) == 1, "expected one real boundary component")
    boundary = regina.Triangulation3(real[0].build())
    link = source.vertex(vertex_index).buildLink()
    iso = boundary.isIsomorphicTo(link)
    require(iso is not None, "new boundary is not isomorphic to original ideal link")
    save3(boundary, out / "finite_boundary.rga")
    save3(link, out / "original_ideal_link.rga")
    boundary_iso = {"source_boundary_tetrahedron": list(range(boundary.size())),
                    "target_original_link_tetrahedron": [iso.simpImage(i) for i in range(boundary.size())],
                    "vertex_permutation": [[iso.facetPerm(i)[j] for j in range(4)] for i in range(boundary.size())]}
    save_json(out / "BOUNDARY_ISOMORPHISM.json", boundary_iso)
    report = {"status": "PASS_EXPLICIT_FINITE_TRUNCATION", "source": str(source_path),
              "source_sha256": digest(source_path), "finite_sha256": digest(out / "finite_staircase.rga"),
              "source_pentachora": source.size(), "finite_pentachora": result.size(),
              "finite_vertices": result.countVertices(), "boundary_tetrahedra": boundary.size(),
              "boundary_is_exact_original_vertex_link_up_to_explicit_isomorphism": True,
              "input_vertex_links_audited": True,
              "finite_manifold_follows_from_explicit_truncation": True,
              "native_regina_truncateIdeal_used": False,
              "closed_manifold_constructed": False, "standardness_asserted": False,
              "seconds": time.monotonic()-started,
              "code_sha256": digest(__file__)}
    save_json(out / "FINITE_TRUNCATION_VERIFICATION.json", report)
    print(json.dumps(report), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--vertex", type=int)
    args = parser.parse_args()
    construct(args.source, args.out, args.vertex)

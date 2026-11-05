#!/usr/bin/env python3
"""Explicit reference upper handlebody and 3D-Pachner-to-4D layering.

The genus-two upper piece is the boundary connected sum of two explicit
triangulations of Delta^3 x S^1.  A logged 3-dimensional Pachner path on
its real boundary is realised by gluing one pentachoron per move.
This file never fills an S1xS2 boundary by a singular cone.
"""

from __future__ import annotations

import argparse
import hashlib
import gzip
import json
from pathlib import Path
import time

import regina


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def save4(tri, path):
    require(regina.PacketOfTriangulation4(tri).save(str(path)), "failed to save 4D triangulation")


def save3(tri, path):
    require(regina.PacketOfTriangulation3(tri).save(str(path)), "failed to save 3D triangulation")


def get_face(tri, dimension, index):
    return tri.simplex(index) if dimension == 3 else tri.face(dimension, index)


def all_faces(tri, dimension):
    return tri.simplices() if dimension == 3 else tri.faces(dimension)


def d3_times_s1():
    """Three-edge circle times one tetrahedron, with staircase triangulation."""
    tri = regina.Triangulation4()
    unpaired = {}
    simplex_vertices = []
    for t in range(3):
        for split in range(4):
            vertices = [(t, j) for j in range(split + 1)]
            vertices += [((t + 1) % 3, j) for j in range(split, 4)]
            require(len(vertices) == len(set(vertices)) == 5, "degenerate product simplex")
            simplex = tri.newSimplex()
            simplex_vertices.append(vertices)
            for facet in range(5):
                key = tuple(sorted(vertices[:facet] + vertices[facet + 1:]))
                if key not in unpaired:
                    unpaired[key] = (simplex.index(), facet, vertices)
                else:
                    old_index, old_facet, old_vertices = unpaired.pop(key)
                    permutation = [
                        old_facet if v == facet else old_vertices.index(vertices[v])
                        for v in range(5)
                    ]
                    simplex.join(facet, tri.simplex(old_index), regina.Perm5(*permutation))
    require(tri.size() == 12 and tri.isValid() and not tri.isIdeal(),
            "product did not produce a finite 4-manifold")
    require(tri.countBoundaryComponents() == 1 and tri.countBoundaryFacets() == 36,
            "unexpected product boundary")
    return tri, simplex_vertices


def boundary_data(tri):
    require(tri.countBoundaryComponents() == 1, "expected exactly one boundary component")
    boundary = tri.boundaryComponent(0)
    require(boundary.isReal(), "expected a real, rather than ideal, boundary")
    triangulation = regina.Triangulation3(boundary.build())
    locations = []
    for i in range(boundary.size()):
        embedding = boundary.facet(i).embedding(0)
        p = embedding.vertices()
        locations.append({
            "pentachoron": embedding.simplex().index(),
            "facet": p[4], "vertices": [p[j] for j in range(4)],
        })
    return triangulation, locations


def reference_cap(simplify=True):
    product, coordinates = d3_times_s1()
    boundary, locations = boundary_data(product)
    face = locations[0]
    p = product.simplex(face["pentachoron"])
    require(len({p.vertex(v).index() for v in face["vertices"]}) == 4,
            "connected-sum gluing tetrahedron is not embedded")
    cap = regina.Triangulation4(product)
    offset = cap.size()
    cap.insertTriangulation(product)
    # The copies have identical labelled attaching tetrahedra.
    cap.simplex(face["pentachoron"]).join(
        face["facet"], cap.simplex(offset + face["pentachoron"]), regina.Perm5()
    )
    require(cap.isValid() and not cap.isIdeal() and cap.isConnected(), "bad boundary sum")
    require(cap.countBoundaryComponents() == 1, "boundary sum is not connected")
    proof = {
        "construction": "boundary connected sum of two staircase Delta3 x a three-edge S1",
        "product_simplex_coordinates": coordinates,
        "product_pentachora": product.size(), "product_boundary_tetrahedra": boundary.size(),
        "embedded_connecting_tetrahedron": face,
        "raw_cap_isosig": cap.isoSig(), "raw_cap_pentachora": cap.size(),
        "raw_cap_boundary_tetrahedra": cap.countBoundaryFacets(),
    }
    if simplify:
        cap.simplify()
    proof["native_regina_simplification_used"] = bool(simplify)
    proof["cap_isosig"] = cap.isoSig()
    proof["cap_pentachora"] = cap.size()
    proof["boundary_isosig"] = cap.boundaryComponent(0).build().isoSig()
    return cap, proof


class DSU:
    def __init__(self, elements):
        self.parent = {x: x for x in elements}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def local_pachner_labels(tri, dimension, index):
    """Five abstract vertices of the elementary 4-simplex trace cobordism."""
    require(dimension in range(4), "bad 3D move dimension")
    face = get_face(tri, dimension, index)
    require(tri.hasPachner(face), "3D Pachner precondition fails")
    if dimension == 3:
        return {index: [0, 1, 2, 3]}
    embeddings = list(face.embeddings())
    tetrahedra = [e.simplex().index() for e in embeddings]
    require(len(set(tetrahedra)) == len(tetrahedra) == 4 - dimension,
            "Pachner star does not contain the expected distinct tetrahedra")
    dsu = DSU((tet, v) for tet in tetrahedra for v in range(4))
    for embedding in embeddings:
        tet = embedding.simplex()
        centre = {embedding.vertices()[v] for v in range(dimension + 1)}
        for facet in set(range(4)) - centre:
            other = tet.adjacentSimplex(facet)
            require(other is not None and other.index() in tetrahedra,
                    "Pachner star is not internally closed around its centre")
            permutation = tet.adjacentGluing(facet)
            for v in set(range(4)) - {facet}:
                dsu.union((tet.index(), v), (other.index(), permutation[v]))
    classes = sorted({dsu.find(x) for x in dsu.parent})
    require(len(classes) == 5, "Pachner cluster does not have five abstract vertices")
    labels = {c: i for i, c in enumerate(classes)}
    result = {tet: [labels[dsu.find((tet, v))] for v in range(4)] for tet in tetrahedra}
    require(all(len(set(vs)) == 4 for vs in result.values()), "degenerate local old tetrahedron")
    require(len({tuple(sorted(vs)) for vs in result.values()}) == len(result),
            "repeated local old tetrahedron")
    return result


def layer_pachner(cap, source_isosig, dimension, index, expected_after=None, audit=True):
    """Attach one 4-simplex to realise a checked move on the cap boundary."""
    canonical = regina.Triangulation3.fromIsoSig(source_isosig)
    boundary, locations = boundary_data(cap)
    iso = canonical.isIsomorphicTo(boundary)
    require(iso is not None, "cap boundary does not match move source")
    labels = local_pachner_labels(canonical, dimension, index)
    gluing_records = []
    for old_tet, old_labels in sorted(labels.items()):
        boundary_tet = iso.simpImage(old_tet)
        vertex_map = iso.facetPerm(old_tet)
        location = locations[boundary_tet]
        missing = (set(range(5)) - set(old_labels)).pop()
        permutation = [None] * 5
        for v, label in enumerate(old_labels):
            permutation[label] = location["vertices"][vertex_map[v]]
        permutation[missing] = location["facet"]
        require(sorted(permutation) == list(range(5)), "invalid 4D gluing permutation")
        gluing_records.append({
            "new_facet": missing, "old_pentachoron": location["pentachoron"],
            "old_facet": location["facet"], "permutation": permutation,
            "source_tetrahedron": old_tet, "source_local_labels": old_labels,
        })
    require(canonical.pachner(get_face(canonical, dimension, index)), "move guard changed unexpectedly")
    target = canonical.isoSig()
    if expected_after is not None:
        require(target == expected_after, "3D move does not yield the expected state")
    simplex = cap.newSimplex()
    new_index = simplex.index()
    for item in gluing_records:
        simplex.join(item["new_facet"], cap.simplex(item["old_pentachoron"]),
                     regina.Perm5(*item["permutation"]))
    after, _ = boundary_data(cap)
    require(after.isoSig() == target, "4D layering and 3D Pachner move disagree")
    require(cap.hasShellBoundary(cap.simplex(new_index)),
            "layering cannot be undone by a guarded 4D boundary shelling")
    if audit:
        require(cap.isValid() and not cap.isIdeal() and cap.isConnected() and cap.isOrientable(),
                "layering produced an invalid/nonfinite/disconnected/nonorientable cap")
    return {
        "before": source_isosig, "dimension": dimension, "face": index,
        "after": target, "added_pentachoron": new_index, "gluings": gluing_records,
        "boundary_tetrahedra_after": after.size(),
        "inverse_4d_boundary_shelling_guard_passed": True,
    }


def cap_audit(cap):
    require(cap.isValid() and not cap.isIdeal() and cap.isConnected() and cap.isOrientable(),
            "cap is not a finite connected oriented 4-manifold")
    vertices = []
    for vertex in cap.vertices():
        link = regina.Triangulation3(vertex.buildLink())
        original_link_size = link.size()
        link.simplify()
        kind = "ball" if vertex.isBoundary() else "sphere"
        recognised = link.isBall() if vertex.isBoundary() else link.isSphere()
        require(recognised, "cap vertex link fails manifold recognition")
        vertices.append({"index": vertex.index(), "kind": kind,
                         "tetrahedra": original_link_size,
                         "simplified_tetrahedra": link.size(),
                         "native_link_simplification_used": True,
                         "recognised": bool(recognised)})
    boundary, _ = boundary_data(cap)
    summands = boundary.summands()
    reference = regina.Example3.s2xs1()
    require(len(summands) == 2 and all(s.isIsomorphicTo(reference) is not None for s in summands),
            "reference cap boundary is not #2(S1xS2)")
    return {
        "status": "PASS", "pentachora": cap.size(), "vertices": cap.countVertices(),
        "isosig": cap.isoSig(), "H1": str(cap.homology()), "H2": str(cap.homology(2)),
        "fundamental_group": str(cap.fundamentalGroup()),
        "boundary_tetrahedra": boundary.size(), "boundary_isosig": boundary.isoSig(),
        "boundary_summands": [s.isoSig() for s in summands], "vertex_links": vertices,
    }


def self_test(output):
    output = Path(output)
    require(not output.exists(), "output directory already exists")
    output.mkdir(parents=True)
    cap, proof = reference_cap(simplify=True)
    save4(cap, output / "reference_cap.rga")
    boundary, locations = boundary_data(cap)
    save3(boundary, output / "reference_boundary.rga")
    save_json(output / "reference_boundary_locations.json", locations)
    save_json(output / "REFERENCE_CONSTRUCTION.json", proof)
    before_audit = cap_audit(cap)
    # Exercise a 1-4/4-1 pair and a 2-3/3-2 pair by actual 4D layering.
    events = []
    for dimension in (3, 2):
        start = boundary_data(cap)[0].isoSig()
        tri = regina.Triangulation3.fromIsoSig(start)
        candidate = next((f.index() for f in all_faces(tri, dimension) if tri.hasPachner(f)), None)
        require(candidate is not None, "fixture has no required growth move")
        forward = layer_pachner(cap, start, dimension, candidate)
        events.append(forward)
        inverse_dimension = 3 - dimension
        following = regina.Triangulation3.fromIsoSig(forward["after"])
        inverse = None
        for face in all_faces(following, inverse_dimension):
            if not following.hasPachner(face):
                continue
            test = regina.Triangulation3(following)
            require(test.pachner(get_face(test, inverse_dimension, face.index())), "inverse move disappeared")
            if test.isoSig() == start:
                inverse = face.index()
                break
        require(inverse is not None, "could not find inverse fixture move")
        events.append(layer_pachner(cap, forward["after"], inverse_dimension, inverse, start))
    after_audit = cap_audit(cap)
    require(before_audit["boundary_isosig"] == after_audit["boundary_isosig"],
            "fixture inverse pairs did not restore the boundary")
    save4(cap, output / "layered_test_cap.rga")
    save_json(output / "layering_moves.json", events)
    report = {"status": "PASS", "construction": proof, "before": before_audit,
              "after": after_audit, "moves": len(events),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    save_json(output / "UPPER_CAP_SELFTEST.json", report)
    print(json.dumps({"status": "PASS", "moves": len(events),
                      "reference_pentachora": before_audit["pentachora"],
                      "reference_boundary_tetrahedra": before_audit["boundary_tetrahedra"]}), flush=True)


def layer_path(cap_file, path_file, output):
    output = Path(output)
    require(not output.exists(), "output directory already exists")
    output.mkdir(parents=True)
    cap_file, path_file = Path(cap_file), Path(path_file)
    packet = regina.open(str(cap_file))
    require(isinstance(packet, regina.Triangulation4), "source file is not a 4D triangulation")
    cap = regina.Triangulation4(packet)
    original = regina.Triangulation4(cap)
    path = json.loads(path_file.read_text())
    require(boundary_data(cap)[0].isoSig() == path["start"], "path starts at another boundary")
    initial_pentachora = cap.size()
    now = time.monotonic()
    before = path["start"]
    with gzip.open(output / "layering_events.jsonl.gz", "wt") as stream:
        for i, move in enumerate(path["moves"]):
            require(move["before"] == before, f"path discontinuity at {i}")
            event = layer_pachner(cap, before, move["dimension"], move["face"],
                                 move["after"], audit=False)
            event["step"] = i
            stream.write(json.dumps(event, sort_keys=True) + "\n")
            before = event["after"]
            if (i + 1) % 100 == 0:
                print(json.dumps({"layered_moves": i + 1, "pentachora": cap.size(),
                                  "seconds": time.monotonic() - now}), flush=True)
    require(before == path["end"], "wrong final boundary")
    shellback = regina.Triangulation4(cap)
    shelling_indices = []
    while shellback.size() > original.size():
        index = shellback.size() - 1
        require(shellback.shellBoundary(shellback.simplex(index)),
                "guarded reverse 4D boundary shelling failed")
        shelling_indices.append(index)
    require(shellback == original, "reverse shellings did not recover exact source piece")
    audit = cap_audit(cap)
    require(cap.size() == initial_pentachora + len(path["moves"]),
            "unexpected number of layering pentachora")
    save4(cap, output / "cap_with_matched_boundary.rga")
    boundary, locations = boundary_data(cap)
    save3(boundary, output / "matched_boundary.rga")
    save_json(output / "boundary_locations.json", locations)
    report = {"status": "PASS_EXPLICIT_4D_LAYERING", "source_cap_sha256":
              hashlib.sha256(cap_file.read_bytes()).hexdigest(),
              "source_path_sha256": hashlib.sha256(path_file.read_bytes()).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "moves": len(path["moves"]), "start_boundary": path["start"],
              "end_boundary": path["end"], "initial_pentachora": initial_pentachora,
              "elapsed_seconds": time.monotonic() - now, "final_cap_audit": audit,
              "guarded_4d_shellback_indices": shelling_indices,
              "shellback_recovers_exact_source_piece": True}
    save_json(output / "LAYERING_VERIFICATION.json", report)
    print(json.dumps({"status": report["status"], "moves": len(path["moves"]),
                      "pentachora": cap.size(), "boundary_tetrahedra": boundary.size(),
                      "seconds": report["elapsed_seconds"]}), flush=True)
    return report


def glue_pieces(first_file, second_file, output):
    """Close two finite pieces along an explicit boundary isomorphism."""
    output = Path(output)
    require(not output.exists(), "output directory already exists")
    output.mkdir(parents=True)
    first_file, second_file = Path(first_file), Path(second_file)
    first_packet, second_packet = regina.open(str(first_file)), regina.open(str(second_file))
    require(isinstance(first_packet, regina.Triangulation4) and
            isinstance(second_packet, regina.Triangulation4), "inputs must be 4D triangulations")
    first, second = regina.Triangulation4(first_packet), regina.Triangulation4(second_packet)
    require(first.isValid() and second.isValid() and not first.isIdeal() and not second.isIdeal(),
            "gluing inputs are not valid finite pieces")
    first_boundary, first_locations = boundary_data(first)
    second_boundary, second_locations = boundary_data(second)
    isomorphism = first_boundary.isIsomorphicTo(second_boundary)
    require(isomorphism is not None, "piece boundaries are not combinatorially isomorphic")
    maps = []
    for i, place in enumerate(first_locations):
        j = isomorphism.simpImage(i)
        other = second_locations[j]
        vertices = isomorphism.facetPerm(i)
        permutation = [None] * 5
        for v in range(4):
            permutation[place["vertices"][v]] = other["vertices"][vertices[v]]
        permutation[place["facet"]] = other["facet"]
        require(sorted(permutation) == list(range(5)), "invalid boundary gluing permutation")
        maps.append({"source_boundary_tetrahedron": i, "target_boundary_tetrahedron": j,
                     "source_pentachoron": place["pentachoron"], "source_facet": place["facet"],
                     "target_pentachoron": other["pentachoron"], "target_facet": other["facet"],
                     "permutation": permutation})
    result = regina.Triangulation4(first)
    offset = result.size()
    result.insertTriangulation(second)
    result.unlockAll()
    for row in maps:
        result.simplex(row["source_pentachoron"]).join(
            row["source_facet"], result.simplex(offset + row["target_pentachoron"]),
            regina.Perm5(*row["permutation"])
        )
    require(result.isValid() and result.isClosed() and not result.isIdeal() and
            result.isConnected() and result.isOrientable(), "gluing is not a closed oriented 4-manifold")
    links = []
    for vertex in result.vertices():
        link = regina.Triangulation3(vertex.buildLink())
        original_link_size = link.size()
        link.simplify()
        sphere = link.isSphere()
        require(sphere, "closed gluing has a vertex link that is not a sphere")
        links.append({"vertex": vertex.index(), "tetrahedra": original_link_size,
                      "simplified_tetrahedra": link.size(),
                      "native_link_simplification_used": True,
                      "is_sphere": bool(sphere)})
    save4(result, output / "closed_triangulation.rga")
    save_json(output / "boundary_gluing.json", maps)
    report = {
        "status": "PASS_EXPLICIT_CLOSED_4_MANIFOLD",
        "source_first_sha256": hashlib.sha256(first_file.read_bytes()).hexdigest(),
        "source_second_sha256": hashlib.sha256(second_file.read_bytes()).hexdigest(),
        "source_first_pentachora": first.size(), "source_second_pentachora": second.size(),
        "boundary_tetrahedra_glued": len(maps), "pentachora": result.size(),
        "isosig": result.isoSig(), "H1": str(result.homology()), "H2": str(result.homology(2)),
        "fundamental_group": str(result.fundamentalGroup()),
        "euler_characteristic": result.eulerCharManifold(), "vertex_links": links,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "standardness_asserted": False,
    }
    save_json(output / "CLOSED_GLUING_VERIFICATION.json", report)
    print(json.dumps({"status": report["status"], "pentachora": result.size(),
                      "vertices": result.countVertices(), "H1": report["H1"], "H2": report["H2"],
                      "euler_characteristic": report["euler_characteristic"]}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--cap")
    parser.add_argument("--path")
    parser.add_argument("--glue-first")
    parser.add_argument("--glue-second")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if args.self_test:
        self_test(args.out)
    elif args.glue_first or args.glue_second:
        require(args.glue_first and args.glue_second, "supply both gluing inputs")
        glue_pieces(args.glue_first, args.glue_second, args.out)
    else:
        require(args.cap and args.path, "supply --cap and --path, or --self-test")
        layer_path(args.cap, args.path, args.out)

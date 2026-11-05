#!/usr/bin/env python3
"""Independent local-cell verifier for a saved staircase truncation.

This checker does not invoke the constructor or native truncateIdeal().
It checks every saved simplex against the four standard staircase cells,
and checks every facet against either the other cell of the same prism,
an original facet pairing, or the new vertex-link boundary.
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


def verify(source_path, folder):
    now = time.monotonic()
    source_path, folder = Path(source_path), Path(folder)
    source = regina.Triangulation4(regina.open(str(source_path)))
    target = regina.Triangulation4(regina.open(str(folder / "finite_staircase.rga")))
    proof = json.loads((folder / "STAIRCASE_CONSTRUCTION.json").read_text())
    report = json.loads((folder / "FINITE_TRUNCATION_VERIFICATION.json").read_text())
    require(hashlib.sha256(source_path.read_bytes()).hexdigest() == report["source_sha256"], "wrong source digest")
    require(hashlib.sha256((folder / "finite_staircase.rga").read_bytes()).hexdigest() == report["finite_sha256"], "wrong target digest")
    require(target.size() == 4 * source.size(), "wrong cell count")
    rows = [[tuple(v) for v in row] for row in proof["new_simplex_symbolic_vertices"]]
    require(len(rows) == target.size(), "wrong symbolic cell count")
    vertex = proof["chosen_vertex"]
    facet_at = {}
    apex_by_simplex = []
    for old in source.simplices():
        oi = old.index()
        globals_ = [old.vertex(v).index() for v in range(5)]
        require(len(set(globals_)) == 5, "nonembedded original simplex vertices")
        apex = [v for v in range(5) if globals_[v] == vertex]
        require(len(apex) == 1, "not exactly one chosen corner")
        apex_by_simplex.append(apex[0])
        order = sorted((v for v in range(5) if v != apex[0]), key=lambda v: globals_[v])
        for split in range(4):
            si = 4*oi + split
            expected = [(v,0) for v in order[:split+1]] + [(v,1) for v in order[split:]]
            require(rows[si] == expected, "cell is not the corresponding staircase simplex")
            for f in range(5):
                face = frozenset(rows[si][:f] + rows[si][f+1:])
                facet_at.setdefault((oi,face), []).append((si,f))
    count = {"internal_directed_facets":0, "old_gluing_directed_facets":0, "new_boundary_facets":0,
             "old_boundary_facets":0}
    for si, row in enumerate(rows):
        oi = si // 4
        old = source.simplex(oi)
        for f in range(5):
            face = frozenset(row[:f] + row[f+1:])
            local = facet_at[oi,face]
            gluing_map = lambda v: v
            if len(local) == 2:
                mate = next(x for x in local if x != (si,f))
                count["internal_directed_facets"] += 1
            else:
                require(local == [(si,f)], "facet occurs in more than two local cells")
                levels = {lev for col,lev in face}
                colours = {col for col,lev in face}
                if levels == {0}:
                    require(len(face) == len(colours) == 4, "bad truncation facet")
                    require(target.simplex(si).adjacentSimplex(f) is None, "new truncation boundary was glued")
                    count["new_boundary_facets"] += 1
                    continue
                oldface = apex_by_simplex[oi] if levels == {1} else next(v for v in range(5) if v not in colours and v != apex_by_simplex[oi])
                oldmate = old.adjacentSimplex(oldface)
                if oldmate is None:
                    require(target.simplex(si).adjacentSimplex(f) is None, "old boundary facet was glued")
                    count["old_boundary_facets"] += 1
                    continue
                oldperm = old.adjacentGluing(oldface)
                mapped_face = frozenset((oldperm[col],lev) for col,lev in face)
                candidates = facet_at.get((oldmate.index(),mapped_face), [])
                require(len(candidates) == 1, "ambiguous or absent gluing target cell")
                mate = candidates[0]
                gluing_map = lambda v, p=oldperm: (p[v[0]],v[1])
                count["old_gluing_directed_facets"] += 1
            ms, mf = mate
            actual = target.simplex(si).adjacentSimplex(f)
            require(actual is not None and actual.index() == ms, "wrong adjacent simplex")
            perm = target.simplex(si).adjacentGluing(f)
            require(perm[f] == mf, "wrong adjacent facet")
            for j in range(5):
                if j != f:
                    require(rows[ms][perm[j]] == gluing_map(row[j]), "facet vertex map does not agree with prism gluing")
    require(count["new_boundary_facets"] == source.size(), "wrong number of new boundary facets")
    require(sum(count.values()) == 5*target.size(), "not all facets were accounted for")
    return {"status":"PASS_INDEPENDENT_PRISM_FACET_CHECK", "source_pentachora":source.size(),
            "target_pentachora":target.size(), "facet_counts":count,
            "seconds":time.monotonic()-now,
            "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "constructor_imported":False, "native_truncateIdeal_called":False}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--folder", required=True)
    p.add_argument("--out")
    a = p.parse_args()
    result = verify(a.source, a.folder)
    if a.out:
        Path(a.out).write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result), flush=True)

#!/usr/bin/env python3
"""Reconstruct the marked X0 triangulation from the fixed construction recipe.

The default rebuilds finite geometry and replays the recorded three-dimensional
compression choices. It does not rerun the searches that chose those moves or
the native source-identification transport. All writes go to a fresh --output
directory. The accepted source and reference files are read-only inputs.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
import gc
import gzip
import hashlib
import importlib.util
from itertools import combinations
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source"
REFERENCE = ROOT / "reference" / "x0"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    path = Path(path)
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def normalized(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): normalized(v) for k, v in value.items()}
    if isinstance(value, set):
        return sorted((normalized(v) for v in value), key=repr)
    if isinstance(value, (list, tuple)):
        return [normalized(v) for v in value]
    return value


def canonical_hash(value):
    raw = json.dumps(normalized(value), sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(normalized(data), indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "Cannot import " + str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def same_object(actual, expected, description):
    actual = normalized(actual)
    expected = normalized(expected)
    if actual != expected:
        keys = []
        if isinstance(actual, dict) and isinstance(expected, dict):
            keys = sorted(k for k in set(actual) | set(expected)
                          if actual.get(k) != expected.get(k))
        raise ValueError(description + " differs from the fixed recipe; top-level keys: " + repr(keys))
    return actual


def order_like(value, template):
    """Preserve serialization order without taking any data from the template."""
    if isinstance(value, dict) and isinstance(template, dict):
        require(set(value) == set(template), "Serialization template keys differ")
        return {k: order_like(value[k], template[k]) for k in template}
    if isinstance(value, list) and isinstance(template, list):
        require(len(value) == len(template), "Serialization template length differs")
        return [order_like(a, b) for a, b in zip(value, template)]
    return value


def serialize_like(value, template_path):
    """Serialize freshly computed data in the recorded JSON/gzip format.

    The template supplies whitespace, key order, and the gzip platform byte,
    never a mathematical value. Different zlib implementations may still give
    different compressed bytes; parsed-object equality remains checked.
    """
    template_path = Path(template_path)
    original = template_path.read_bytes()
    payload = gzip.decompress(original) if template_path.suffix == ".gz" else original
    template = json.loads(payload)
    value = order_like(normalized(value), template)
    stripped = payload.rstrip(b"\n")
    if b"\n" not in stripped:
        separators = (",", ":") if b'": ' not in stripped else (", ", ": ")
        result = json.dumps(value, separators=separators).encode()
    else:
        match = re.search(rb"\n( +)\S", stripped)
        indentation = len(match.group(1)) if match else 2
        result = json.dumps(value, indent=indentation).encode()
    if payload.endswith(b"\n"):
        result += b"\n"
    if template_path.suffix == ".gz":
        encoded = bytearray(gzip.compress(result, compresslevel=9, mtime=0))
        require(original[3] == 0, "Unsupported gzip header with optional fields")
        encoded[9] = original[9]
        return bytes(encoded)
    return result


def verify_inventory():
    path = ROOT / "provenance" / "source_files.json"
    inventory = load(path)
    files = inventory["files"]
    require(isinstance(files, list) and files, "Missing construction inventory")
    snapshot = {}
    for row in files:
        target = ROOT / row["path"]
        resolved = target.resolve()
        require(ROOT in resolved.parents and not target.is_symlink(), "Unsafe inventory path")
        require(target.is_file() and target.stat().st_size == row["bytes"],
                "Missing or wrong-sized input: " + row["path"])
        digest = file_hash(target)
        require(digest == row["sha256"], "Input hash mismatch: " + row["path"])
        snapshot[row["path"]] = digest
    snapshot[str(path.relative_to(ROOT))] = file_hash(path)
    return snapshot


class Build:
    def __init__(self, output, level, keep_work):
        self.output = output
        self.work = output / "work"
        self.frozen = self.work / "frozen_release"
        self.components = output / "components"
        self.level = level
        self.keep_work = keep_work
        self.report = {
            "schema": "x0-construction-build-v1",
            "status": "RUNNING",
            "level": level,
            "parameters": {"pretzel": [-2, 3, 7], "twist": 0, "roll": 1},
            "python": sys.version,
            "builder_sha256": file_hash(__file__),
            "stages": [],
            "components": [],
            "retained_recipe_inputs": [
                "Recorded source and adaptor compression move choices and endpoint branching orders",
                "Raw peripheral chart, checked independently as exact finite chart data",
                "Pinned native marked-exterior transport and diagram-identification records",
                "Puncture interface specification, with its geometry checked by the construction audit",
                "Historical X2 source-acceptance receipt and original master hash",
            ],
            "scope": (
                "Fresh finite source/component geometry, saved compression replay, and marked X0 assembly; "
                "native source transport is a fixed input, not a newly executed native path. "
                "This build alone does not prove the manifold identification or PL standardness."
                if level == "geometry" else
                "Fresh marked X0 assembly from the supplied fixed component recipe; upstream geometry "
                "is not regenerated in this level. No manifold-identification or standardness claim."
            ),
        }

    def checkpoint(self):
        atomic_json(self.output / "build_report.json", self.report)

    def stage(self, name, function):
        print("BUILD", name, flush=True)
        started = time.monotonic()
        log = self.output / "logs" / (name + ".log")
        row = {"name": name, "status": "RUNNING", "log": str(log.relative_to(self.output))}
        self.report["stages"].append(row)
        self.checkpoint()
        try:
            with log.open("w") as handle, redirect_stdout(handle), redirect_stderr(handle):
                result = function()
            row["status"] = "PASS"
            if result is not None:
                row["result"] = normalized(result)
        except BaseException as error:
            row["status"] = "FAILED"
            row["error"] = str(error)
            with log.open("a") as handle:
                traceback.print_exc(file=handle)
            raise
        finally:
            row["wall_seconds"] = time.monotonic() - started
            self.checkpoint()
            gc.collect()
        print("PASS", name, round(row["wall_seconds"], 2), "seconds", flush=True)

    def component(self, relative, generated, preserve=(), origin="generated_computation"):
        """Install generated values, retaining only explicitly listed metadata."""
        template = SOURCE / "frozen_release" / relative
        expected = load(template)
        actual = normalized(generated)
        preserved = {}
        for key in preserve:
            require(key in expected and key in actual, "Missing metadata field " + key)
            preserved[key] = {"fresh_value": actual[key], "historical_value": expected[key]}
            actual[key] = expected[key]
        actual = same_object(actual, expected, relative)
        encoded = serialize_like(actual, template)
        target = self.frozen / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(encoded)
        archived = self.components / relative
        archived.parent.mkdir(parents=True, exist_ok=True)
        archived.write_bytes(encoded)
        row = {"path": relative, "origin": origin, "canonical_sha256": canonical_hash(actual),
               "all_fields_equal_to_recipe": True, "file_sha256": file_hash(target),
               "reference_file_sha256": file_hash(template),
               "byte_equal_to_recipe": encoded == template.read_bytes()}
        if preserved:
            row["preserved_historical_metadata"] = preserved
        self.report["components"].append(row)
        self.checkpoint()
        return actual

    def source_geometry(self):
        module = load_module("construction_seifert_builder", self.frozen / "seifert/build_source.py")
        model, cover, quotient = module.build()
        oriented = self.component("seifert/data_oriented/source.json.gz", model)
        self.component("seifert/data_oriented/cover_regina.json.gz", cover)
        self.component("seifert/data_oriented/quotient_regina.json.gz", quotient)
        legacy = load(SOURCE / "frozen_release/coarsened/master_source.json.gz")
        projected = {key: oriented[key] for key in legacy}
        self.component("coarsened/master_source.json.gz", projected)
        self.component("seifert/data/source.json.gz", projected)
        orientation = load(SOURCE / "frozen_release/seifert/data/master_orientation.json")
        self.component("seifert/data/master_orientation.json", {key: oriented[key] for key in orientation})
        return {"vertices": len(oriented["vertices"]), "tetrahedra": len(oriented["facets"]),
                "legacy_source_projection": "Omit later-added orientation and gluing-path metadata only"}

    def source_compression(self):
        sys.path.insert(0, str(self.frozen / "lib"))
        replay = load_module("verify_compact3", self.frozen / "lib/verify_compact3.py")
        source = load(self.frozen / "coarsened/master_source.json.gz")
        certificate = load(self.frozen / "coarsened/coarsening_certificate.json.gz")
        tags = sorted({v for tet in source["facets"] for v in tet}, key=repr)
        initial = [tuple(sorted(tet, key=lambda v: min(v, source["deck"][v])))
                   for tet in source["facets"]]
        deck = {v: source["deck"][v] for v in tags}
        # replay_part checks every initial simplex and the full initial branching;
        # also bind that initial complex directly to the freshly generated source.
        require({tuple(sorted(t)) for t in initial} == {tuple(t) for t in source["facets"]},
                "Compression certificate is not bound to generated source")
        result, receipt = replay.replay_part(certificate, initial, deck,
                                             {v for v in tags if source["deck"][v] == v})
        active = sorted({v for tet in result for v in tet})
        index = {v: i for i, v in enumerate(active)}
        facets = sorted(tuple(sorted(index[v] for v in tet)) for tet in result)
        final_deck = [index[source["deck"][v]] for v in active]
        descendant = {
            "schema": "marked-pretzel-source-descendant-v1", "facets": facets,
            "closed_cover_facets": facets, "deck": final_deck,
            "vertex_ranks": [min(i, final_deck[i]) for i in range(len(final_deck))],
            "parent_vertices": active, "parent_source_sha256": canonical_hash(source),
            "parent_file_sha256": file_hash(self.frozen / "coarsened/master_source.json.gz"),
            "ancestry": "equivariant pairs of independently replayed relative 3-ball replacements; fixed circle pointwise preserved",
            "ramification_edges": sorted({e for tet in facets for e in combinations(tet, 2)
                                           if all(final_deck[v] == v for v in e)}),
            "source_interpretation": source["source_interpretation"],
        }
        self.component("coarsened/source.json.gz", descendant)
        return receipt

    def cut_geometry(self):
        module = load_module("prepare_source", self.frozen / "prepare_source.py")
        cut = module.cut_source(load(self.frozen / "coarsened/source.json.gz"))
        self.component("peripheral/cut_source.json.gz", cut)
        return {"N_tetrahedra": len(cut["N"]), "E_tetrahedra": len(cut["E"]),
                "torus_triangles": len(cut["T"])}

    def orientations_and_chart(self):
        module = load_module("construction_orientation_checker", self.frozen / "verify_source_geometry.py")
        paths = {"master": self.frozen / "coarsened/master_source.json.gz",
                 "mo": self.frozen / "seifert/data/master_orientation.json",
                 "cert": self.frozen / "coarsened/coarsening_certificate.json.gz",
                 "source": self.frozen / "coarsened/source.json.gz"}
        data = {key: load(path) for key, path in paths.items()}
        orientation, record = module.transport(data["master"], data["mo"], data["cert"], data["source"])
        record["input_file_sha256"] = {key: file_hash(path) for key, path in paths.items()}
        self.component("coarsened/source_orientation.json", record)
        cut = load(self.frozen / "peripheral/cut_source.json.gz")
        tori, record = module.clip_orient(data["source"], orientation, cut)
        record["source_orientation_file_sha256"] = file_hash(self.frozen / "coarsened/source_orientation.json")
        self.component("peripheral/cut_orientation.json", record)
        chart_report = module.chart_check(cut, load(self.frozen / "peripheral/peripheral_chart.json.gz"), tori)
        transport = load_module("construction_chart_transport", self.frozen / "transport_chart.py")
        chart = transport.build()
        self.component("peripheral/diagram_chart.json.gz", chart)
        self.component("peripheral/puncture_product_interface.json",
                       load(self.frozen / "peripheral/puncture_product_interface.json"),
                       origin="retained_recipe_with_bindings_refreshed; geometry_checked_by_audit")
        return {"raw_chart_freshly_checked": True, "native_transport_reexecuted": False,
                "chart_checks": chart_report}

    def adaptor_geometry(self):
        adaptor = load_module("torus_adaptor", self.frozen / "torus_adaptor.py")
        stages, deck, record, cut = adaptor.make_adaptor()
        self.component("adaptor/adaptor.json.gz", adaptor.convert(record))
        certificate = load(self.frozen / "adaptor/coarsening_certificate.json.gz")
        chart = load(self.frozen / "peripheral/diagram_chart.json.gz")
        from verify_compact3 import replay
        compact, final_deck, receipt = replay(adaptor.convert(certificate), stages, deck,
                                              adaptor.convert(record), chart)
        expected = load(SOURCE / "frozen_release/adaptor/model.json.gz")
        # Tagged deck maps are dictionaries serialized as pair lists. Restore the
        # fixed pair-list order only after proving equality of the entire map.
        expected_deck = {adaptor.tup(a): adaptor.tup(b) for a, b in expected["deck"]}
        require(final_deck == expected_deck, "Regenerated adaptor deck map differs")
        ordered_deck = [(adaptor.tup(a), final_deck[adaptor.tup(a)]) for a, _ in expected["deck"]]
        model = {
            "schema": "marked-pretzel-rectangular-adaptor-model-v1", "L": 3, "V": 4,
            "stages": compact, "deck": ordered_deck,
            "puncture_boundary": [[("N", "o", v) for v in sorted(face, key=lambda v: min(v, cut["deck"][v]))]
                                  for face in record["boundary"]],
            "source_cut_sha256": adaptor.digest(cut), "diagram_chart_sha256": adaptor.digest(chart),
            "adaptor_sha256": adaptor.digest(record), "coarsened": True,
            "coarsening_certificate_sha256": adaptor.digest(certificate),
            # This is the explicit corrected interpretation recorded by the
            # original orientation-metadata rebind, not the older producer note.
            "source_orientation": "J-oriented source agrees with canonical Seifert master; (u_J,v_J)=(u_old,-v_old), outward boundary N=-du_J wedge dv_J",
        }
        self.component("adaptor/model.json.gz", model)
        receipt["deck_pair_list_order"] = "Fixed recipe order after exact dictionary equality"
        return receipt

    def roll_geometry(self):
        module = load_module("construction_central_roll", self.frozen / "central_roll.py")
        output = self.output / "scratch" / "central_roll"
        certificate, model = module.build(output)
        self.component("central_roll/certificate.json.gz", module.conv(certificate))
        self.component("central_roll/model.json.gz", module.conv(model))
        return {"tetrahedra": len(model["tetrahedra"]), "roll": model["roll"]}

    def cap_geometry(self):
        sys.path.insert(0, str(self.frozen / "cap"))
        module = load_module("build_caps", self.frozen / "cap/build_caps.py")
        output = self.output / "scratch" / "caps"
        module.build(output, n=8, t=4)
        for filename in ("shared_boundary.json.gz", "cap_m0.json.gz", "cap_m2.json.gz"):
            self.component("cap/data/" + filename, load(output / filename))
        cap = load(output / "cap_m0.json.gz")
        return {"m0_pentachora": len(cap["facets"]), "shared_boundary_tetrahedra": len(cap["boundary"])}

    def bridge_geometry(self):
        module = load_module("construction_bridge", self.frozen / "cap/build_bridge.py")
        output = self.output / "scratch" / "bridge"
        module.build(self.frozen / "peripheral/puncture_product_interface.json",
                     self.frozen / "cap/data/shared_boundary.json.gz", output)
        bridge = load(output / "source_cap_bridge.json.gz")
        self.component("cap/data/source_cap_bridge.json.gz", bridge, preserve=("elapsed_seconds",))
        return {"pentachora": len(bridge["facets"]),
                "fresh_elapsed_seconds": bridge["elapsed_seconds"],
                "preserved_metadata": "Only elapsed_seconds retained from the historical recipe after all other fields match"}

    def assembly(self):
        module = load_module("construction_x0_assembler", self.work / "scripts/assemble_candidate.py")
        module.build()
        produced = self.work / "candidate"
        output = self.output / "candidate"
        output.mkdir()
        comparisons = {}
        for filename in ("marked_candidate.json.gz", "cover_regina.json.gz", "quotient_regina.json.gz"):
            actual = load(produced / filename)
            expected = load(REFERENCE / filename)
            same_object(actual, expected, "Fresh final " + filename)
            # Serialize generated data, rather than copying any reference output.
            (output / filename).write_bytes(serialize_like(actual, REFERENCE / filename))
            same_object(load(output / filename), actual, "Serialized fresh final " + filename)
            comparisons[filename] = {"all_fields_equal": True,
                                     "canonical_sha256": canonical_hash(actual),
                                     "file_sha256": file_hash(output / filename),
                                     "reference_file_sha256": file_hash(REFERENCE / filename),
                                     "byte_equal": (output / filename).read_bytes() == (REFERENCE / filename).read_bytes()}
        shutil.copy2(produced / "assembly_receipt.json", output / "assembly_receipt.json")
        master = load(output / "marked_candidate.json.gz")
        require(len(master["facets"]) == 103744 and len(master["vertices"]) == 4310,
                "Unexpected final triangulation size")
        self.report["outputs"] = comparisons
        return {"pentachora": len(master["facets"]), "vertices": len(master["vertices"]),
                "twist": master["twist"], "roll": master["roll"],
                "checked": "Every field of master, cover and quotient objects, including all markings, orientations, parts, ports and gluing permutations"}

    def run(self):
        before = verify_inventory()
        self.report["input_files_checked"] = len(before)
        self.report["input_manifest_sha256"] = file_hash(ROOT / "provenance/source_files.json")
        self.checkpoint()
        shutil.copytree(SOURCE, self.work, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        require(not (self.work / "candidate").exists(), "Source tree contains a candidate output")
        sys.path.insert(0, str(self.frozen))
        sys.path.insert(0, str(self.frozen / "lib"))
        if self.level == "geometry":
            for name, function in [
                ("montesinos_source", self.source_geometry),
                ("source_compression_replay", self.source_compression),
                ("source_cut", self.cut_geometry),
                ("orientations_and_fixed_chart", self.orientations_and_chart),
                ("adaptor_geometry_and_compression", self.adaptor_geometry),
                ("longitudinal_roll_collar", self.roll_geometry),
                ("caps", self.cap_geometry),
                ("source_cap_bridge", self.bridge_geometry),
            ]:
                self.stage(name, function)
        self.stage("marked_x0_assembly", self.assembly)
        after = verify_inventory()
        require(before == after, "Canonical construction inputs changed during the build")
        self.report["canonical_inputs_unchanged"] = True
        self.report["status"] = "PASS"
        if not self.keep_work:
            shutil.rmtree(self.work)
            self.report["disposable_work_tree_removed"] = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory for generated files and the build report")
    parser.add_argument("--level", choices=("geometry", "assembly"), default="geometry",
                        help="geometry rebuilds finite components and replays fixed choices; assembly uses supplied components")
    parser.add_argument("--keep-work", action="store_true", help="Retain the disposable staged source tree for inspection")
    args = parser.parse_args()
    require(__debug__, "Run build.py without Python -O; legacy source checkers require assertions")
    output = args.output.expanduser().resolve()
    require(not output.exists(), "Output already exists; choose a fresh directory")
    require(output != ROOT and output not in ROOT.parents, "Output cannot contain the construction source")
    for protected in (SOURCE, ROOT / "reference", ROOT / "provenance"):
        require(output != protected and protected not in output.parents, "Output is inside protected inputs")
    output.mkdir(parents=True)
    (output / "logs").mkdir()
    build = Build(output, args.level, args.keep_work)
    started = time.monotonic()
    try:
        build.run()
    except BaseException as error:
        build.report["status"] = "FAILED"
        build.report["error"] = str(error)
        raise
    finally:
        build.report["wall_seconds"] = time.monotonic() - started
        build.checkpoint()
    print("PASS: fresh marked X0 generated at", output / "candidate", flush=True)
    print("Build report:", output / "build_report.json", flush=True)


if __name__ == "__main__":
    main()

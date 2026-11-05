#!/usr/bin/env python3
"""Quick, deterministic replay of the three 2026-11-02 redraw transcriptions.

Requires Regina and Spherogram. No simplification, normal-surface recognition,
or exhaustive triangulation search is performed. This checks the encoded PDs,
necessary linking invariants, saved algebraic proofs, and Katie position. It
does not prove that a transcription agrees with the drawing or that either
redraw is isotopic to the earlier diagram.

Run from the repository root:
  python -O code/verify_redraw_transcription.py --out replay_fresh.json
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import regina
from spherogram import Link

from check_handlebody_group import replay, check_dotted_unlink
from check_redraw_markings import witnesses
from orient_katie_diagram import admission, orient


PDF_HASHES = {
    1: "cfc1eb4f163e7c504089dfb96bceaec1ec5ed3074ab97c466db78a751c2eede2",
    2: "2c1099d691aac6e68a7130f0a830eae78b7c223f564ebe9d7f2da62800cc2c11",
}
CASES = (
    ("figure1", 1, "assembled_f1_v1", 203, 0, True),
    ("figure2_literal_paint", 2, "assembled_f2_v1", 181, 0, False),
    ("figure2_dashed_back_arc", 2, "assembled_f2_dashed_v1", 181, 2, True),
)


def require(condition, message):
    # Deliberately not assert: every check remains active with python -O.
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def validate_pd(data):
    pd = data["pd_code"]
    n = data["component_count"]
    m = data["crossing_count"]
    require(len(pd) == m == len(data["crossings"]), "crossing count")
    counts = Counter(x for row in pd for x in row)
    require(all(len(row) == 4 for row in pd), "PD tuple length")
    require(set(counts) == set(range(2*m)), "zero-based consecutive PD labels")
    require(all(v == 2 for v in counts.values()), "each PD arc occurs twice")
    link = Link(pd, check_planarity=True)
    require(len(link.link_components) == n, "Spherogram component count")
    mapping = {}
    for i, c in enumerate(link.crossings):
        expected = data["crossings"][i]
        require(c.label == i == expected["id"], "crossing numbering")
        require(c.sign == expected["sign"], "Spherogram crossing sign")
        for imported, source in ((c.strand_components[0], expected["under_component"]),
                                 (c.strand_components[1], expected["over_component"])):
            require(imported not in mapping or mapping[imported] == source,
                    "Spherogram crossing ownership")
            mapping[imported] = source
    order = [mapping[i] for i in range(n)]
    require(sorted(order) == list(range(n)), "component ownership bijection")
    require(order == data["pd_component_to_source"], "stored component ownership")
    decorations = ["x" if data["components"][i]["dotted"]
                   else data["components"][i]["framing"] for i in order]
    require(decorations == data["decorations"], "component decorations")
    matrix = [[int(x) for x in row] for row in link.linking_matrix()]
    require(matrix == [[data["linking_matrix"][i][j] for j in order] for i in order],
            "saved linking matrix")

    rpd = [[a+1 for a in row] for row in pd]
    require(rpd == data["regina_pd_code"], "stored Regina PD")
    rlink = regina.Link.fromPD(rpd)
    require(rlink.size() == m and rlink.countComponents() == n, "Regina size")
    require(rlink.virtualGenus() == 0, "Regina planarity")
    owner = {}
    for i in range(n):
        start = rlink.component(i)
        at = start
        while True:
            key = (at.crossing().index(), at.strand())
            require(key not in owner, "Regina component repeated strand")
            owner[key] = i
            at = at.next()
            if at == start:
                break
    require(len(owner) == 2*m, "Regina covers all crossing strands")
    doubled = [[0]*n for _ in range(n)]
    for i, c in enumerate(data["crossings"]):
        sign = rlink.crossing(i).sign()
        require(sign == c["sign"], "Regina crossing sign")
        a, b = owner[i, 0], owner[i, 1]
        require(order[a] == c["under_component"] and order[b] == c["over_component"],
                "Regina crossing ownership")
        if a != b:
            doubled[a][b] += sign
            doubled[b][a] += sign
    require(all(x % 2 == 0 for row in doubled for x in row), "linking integrality")
    require([[x//2 for x in row] for row in doubled] == matrix,
            "independent Regina/Spherogram linking matrix agreement")
    return {"crossings": m, "components": n, "planar": True,
            "both_imports_match_every_crossing_sign_and_owner": True,
            "pd_component_to_source": order, "decorations": decorations,
            "linking_matrix_in_pd_order": matrix}


def verify_group(path, native_path, data):
    report = read(path)
    require(report["input_sha256"] == sha(native_path), "group/input hash binding")
    require(report["pd_code"] == data["pd_code"], "group PD binding")
    require(report["decorations"] == data["decorations"], "group marking binding")
    require(data["pd_component_to_source"] == list(range(data["component_count"])),
            "group checker expects identity source/component ordering")
    exterior = report["link_exterior"]
    added = []
    for j, decoration in enumerate(data["decorations"]):
        if decoration in ("x", "."):
            continue
        for kind, word in zip(("meridian", "longitude"), exterior["peripheral_curves"][j]):
            added.append({"component": j, "source_component": j,
                          "curve": kind, "word": word})
    require(added == report["added_relations"], "saved peripheral quotient relations")
    proof = report["tietze_proof"]
    require(proof["initial_generator_count"] == exterior["generator_count"],
            "Tietze generator binding")
    require(proof["initial_relators"] == exterior["relators"] + [a["word"] for a in added],
            "Tietze initial-relator binding")
    require(replay(proof), "Tietze proof did not reach trivial presentation")
    dotted = check_dotted_unlink(data)
    require(dotted["proved_unlink"], "dotted unlink check")
    require(dotted == report["dotted_sublink"], "fresh dotted-sublink replay agreement")
    require(report["status"] == "TRIVIAL_HANDLEBODY_PI1", "saved group status")
    return {"report_sha256": sha(path), "input_hash_binding": True,
            "tietze_steps_replayed": len(proof["steps"]),
            "trivial_quotient_proved_from_saved_presentation": True,
            "dotted_unlink_rechecked": True,
            "peripheral_extraction_recomputed": False,
            "qualification": "Saved presentation and peripheral extraction remain part of the input trust boundary."}


def verify_raw_conversion(root, folder, oriented_path):
    directory = root/"results/redraw_20261102"/folder
    receipt_path = directory/"RUN_RECEIPT.json"
    if not receipt_path.exists():
        return {"status": "NOT_INCLUDED_IN_THIS_COPY"}
    receipt = read(receipt_path)
    raw_path = directory/"raw_katie.rga"
    sig_path = directory/"triangulation.isosig"
    input_path = directory/"input.json"
    require(sha(input_path) == receipt["input_sha256"] == sha(oriented_path),
            "raw Katie/input hash binding")
    require(sha(raw_path) == receipt["raw_katie_sha256"], "raw Katie RGA hash")
    require(sha(sig_path) == receipt["triangulation_signature_sha256"], "raw signature file hash")
    require(receipt["returncode"] == 0, "Katie return code")
    raw = regina.open(str(raw_path))
    require(isinstance(raw, regina.Triangulation4), "raw RGA is a 4-triangulation")
    from_sig = regina.Triangulation4.fromIsoSig(sig_path.read_text().strip())
    require(raw == from_sig, "RGA and signature encode identical simplex gluings")
    basic = {"pentachora": raw.size(), "vertices": raw.countVertices(),
             "valid": raw.isValid(), "connected": raw.isConnected(),
             "orientable": raw.isOrientable(),
             "euler_characteristic_triangulation": raw.eulerCharTri()}
    require(basic == receipt["triangulation_basic"], "raw basic combinatorial checks")
    return {"status": "PASS_RAW_BYTES_AND_BASIC_COMBINATORICS",
            "receipt_sha256": sha(receipt_path), "raw_rga_sha256": sha(raw_path),
            "signature_file_sha256": sha(sig_path), "basic": basic,
            "vertex_link_recognition_performed": False,
            "closed_manifold_certified": False,
            "qualification": "No isIdeal(), isSphere(), normal-surface recognition, filling, or simplification is run."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out", "--output", dest="output", type=Path, default=Path("replay_fresh.json"))
    args = parser.parse_args()
    root = args.root.resolve()
    start = time.monotonic()
    report = {"schema": "motoo.redraw.transcription-replay.v1",
              "versions": {"regina_distribution": importlib.metadata.version("regina"),
                           "regina_engine": regina.versionString(),
                           "spherogram": importlib.metadata.version("spherogram")},
              "source_pdfs": {}, "cases": {},
              "source_diagram_isotopy_proved": False,
              "four_manifold_identification_proved": False,
              "qualification": "PASS means this bounded replay reproduced the stated transcription checks, including their failures. It does not certify source isotopy, upper-handle completion, or standardness."}
    originals = {}
    for figure in (1, 2):
        path = root/f"inputs/redraw_20261102/figure{figure}.pdf"
        require(sha(path) == PDF_HASHES[figure], "raw source PDF hash")
        report["source_pdfs"][str(figure)] = {"path": str(path.relative_to(root)), "sha256": sha(path)}
        old_path = root/f"inputs/figure{figure}_link.json"
        old = read(old_path)
        originals[figure] = (old, validate_pd(old), sha(old_path))
    for name, figure, folder, count, expected_matches, has_oriented in CASES:
        directory = root/"scratch/redraw_20261102"/folder
        path = directory/"native_checked_link.json"
        data = read(path)
        require(data["figure"] == figure and data["crossing_count"] == count,
                "expected figure/crossing count")
        checked = validate_pd(data)
        crossing_path = directory/"crossing_record.json"
        require(data["crossing_record_sha256"] == sha(crossing_path), "crossing record binding")
        old, old_check, old_hash = originals[figure]
        matches = witnesses(old_check["linking_matrix_in_pd_order"], old["decorations"],
                            checked["linking_matrix_in_pd_order"], data["decorations"])
        require(len(matches) == expected_matches, "expected necessary-condition witnesses")
        status = "PASS_LINKING_NECESSARY_CONDITION" if matches else "LINKING_MISMATCH"
        result = {"input_path": str(path.relative_to(root)), "input_sha256": sha(path),
                  "source_pdf_sha256": PDF_HASHES[figure],
                  "crossing_record_sha256": sha(crossing_path),
                  "pd_checks": checked,
                  "original_transcription_sha256": old_hash,
                  "marked_linking_status": status, "marked_linking_witnesses": matches,
                  "isotopy_proved": False,
                  "group": verify_group(directory/"group.json", path, data)}
        marking_file = directory/"marking_check.json"
        if marking_file.exists():
            saved = read(marking_file)
            require(saved["status"] == status and saved["redraw_sha256"] == sha(path)
                    and saved["original_sha256"] == old_hash, "saved marked-invariant bindings")
            stripped = [{k: v for k, v in w.items() if k != "new_pd_component_to_old_source_id"}
                        for w in saved["witnesses"]]
            require(stripped == matches, "fresh marked-invariant witness replay")
        if has_oriented:
            oriented_path = directory/"katie_oriented.json"
            saved = read(oriented_path)
            require(saved["source_record_sha256"] == sha(path), "Katie/source hash binding")
            minimal_path = directory/"minimal_katie_input.json"
            minimal = read(minimal_path)
            require(minimal["pd_code"] == data["pd_code"]
                    and minimal["decorations"] == data["decorations"], "Katie minimal input")
            require(saved["geometry_admission"]["input_sha256"] == sha(minimal_path),
                    "Katie orientation input hash binding")
            fresh = orient({"pd_code": data["pd_code"], "decorations": data["decorations"]})
            require(fresh["pd_code"] == saved["pd_code"], "exact Katie orientation replay PD")
            require(fresh["decorations"] == saved["decorations"], "Katie decorations")
            saved_geometry = {k: v for k, v in saved["geometry_admission"].items()
                              if k != "input_sha256"}
            require(fresh["geometry_admission"] == saved_geometry,
                    "exact Katie geometry admission replay")
            require(fresh["geometry_admission"]["dotted_position_admitted"], "Katie position")
            result["katie_position"] = {"status": "PASS_AFTER_AUXILIARY_ORIENTATION_ONLY",
                "oriented_file_sha256": sha(oriented_path),
                "reversed_components": fresh["geometry_admission"]["reversed_components"]}
            raw_folder = "katie_figure1_source_v1" if figure == 1 else "katie_figure2_dashed_v1"
            result["raw_katie_conversion"] = verify_raw_conversion(root, raw_folder, oriented_path)
        else:
            dotted = [i for i, d in enumerate(data["decorations"]) if d in ("x", ".")]
            try:
                diagnostic = admission(data["pd_code"], dotted)
                result["katie_position"] = {"status": diagnostic["status"]}
            except ValueError as error:
                result["katie_position"] = {"status": "NOT_ADMITTED", "reason": str(error)}
        report["cases"][name] = result
        print(json.dumps({"case": name, "crossings": count, "linking": status,
                          "witnesses": len(matches), "group_replay": "PASS",
                          "katie": result["katie_position"]["status"]}), flush=True)
    report["status"] = "PASS_BOUNDED_TRANSCRIPTION_REPLAY"
    report["wall_seconds"] = round(time.monotonic()-start, 3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "output": str(args.output),
                      "wall_seconds": report["wall_seconds"]}))


if __name__ == "__main__":
    main()

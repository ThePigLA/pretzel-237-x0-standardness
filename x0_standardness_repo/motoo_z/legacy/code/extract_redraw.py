#!/usr/bin/env python3
"""Inventory and separate the exact coloured geometry of a tablet Kirby redraw.

This conservatively exports PDF filled-path geometry and diagnostic layer
renders.  It does NOT infer centre lines, join drawing gaps, assign crossings,
produce a PD code, or certify isotopy.  The original clean-vector extractor
expects stroked centre lines and must not be applied to these pen outlines.

Usage:
  python code/extract_redraw.py inputs/redraw_20261102/figure1.pdf --out NEW_DIR

Only PyMuPDF 1.26.6 and the Python standard library are required.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path

import fitz


PALETTE = {
    (0.964706, 0.325490, 0.325490): "salmon",
    (0.011765, 0.274495, 0.560769): "navy",
    (0.819608, 0.827451, 0.000000): "yellow",
    (0.509804, 0.129412, 0.545098): "purple",
    (0.090196, 0.850980, 1.000000): "cyan",
    (0.000000, 0.478431, 1.000000): "blue",
    (0.823529, 0.000000, 0.000000): "red",
    (0.494102, 0.776471, 0.211749): "green",
    (0.960784, 0.603922, 0.137240): "orange",
}


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def colour_name(colour):
    key = tuple(round(float(x), 6) for x in colour)
    return PALETTE.get(key, "rgb_" + "_".join(f"{x:.6f}" for x in key))


def path_record(index, drawing):
    require(drawing["type"] == "f", "Unexpected non-filled path; inspect separately")
    items = []
    for item in drawing["items"]:
        require(item[0] in ("l", "c"), "Unexpected vector primitive")
        items.append({"kind": item[0], "points": [[float(p.x), float(p.y)] for p in item[1:]]})
    return {
        "drawing_index": index,
        "paint_sequence": drawing["seqno"],
        "kind": "filled_pen_outline",
        "colour": [float(x) for x in drawing["fill"]],
        "colour_name": colour_name(drawing["fill"]),
        "bbox": list(drawing["rect"]),
        "fill_opacity": drawing["fill_opacity"],
        "even_odd": drawing["even_odd"],
        "close_path": drawing["closePath"],
        "items": items,
    }


def paint(page, paths):
    """Repaint the extracted filled paths, retaining their original order."""
    for record in paths:
        shape = page.new_shape()
        for item in record["items"]:
            points = [fitz.Point(*p) for p in item["points"]]
            if item["kind"] == "c":
                shape.draw_bezier(*points)
            else:
                shape.draw_line(*points)
        shape.finish(color=None, fill=tuple(record["colour"]),
                     fill_opacity=record["fill_opacity"],
                     even_odd=record["even_odd"],
                     closePath=record["close_path"])
        shape.commit()


def inspect(source, output, scale):
    source, output = Path(source).resolve(), Path(output).resolve()
    require(not output.exists(), "Use a new output directory")
    require(importlib.metadata.version("PyMuPDF") == "1.26.6", "Requires PyMuPDF 1.26.6")
    document = fitz.open(source)
    require(len(document) == 1, "This intake expects one drawing page")
    output.mkdir(parents=True)
    page = document[0]
    drawings = page.get_drawings()
    records = [path_record(i, drawing) for i, drawing in enumerate(drawings)]
    colours = sorted({r["colour_name"] for r in records})
    geometry = {
        "schema": "motoo.coloured-tablet-geometry.v1",
        "source_sha256": digest(source),
        "page_rect": list(page.rect),
        "coordinates": "PDF page coordinates, x right and y down, in points",
        "paths": records,
        "centre_lines_extracted": False,
        "gap_pairings_inferred": False,
        "crossings_assigned": False,
        "pd_generated": False,
    }
    with gzip.open(output / "filled_paths.json.gz", "wt") as stream:
        json.dump(geometry, stream, separators=(",", ":"))
        stream.write("\n")

    rows = []
    for name in colours:
        selected = [r for r in records if r["colour_name"] == name]
        layer = fitz.open()
        target = layer.new_page(width=page.rect.width, height=page.rect.height)
        paint(target, selected)
        image_path = output / f"layer_{name}.png"
        target.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False).save(image_path)
        boxes = [r["bbox"] for r in selected]
        rows.append({
            "colour_name": name,
            "rgb": selected[0]["colour"],
            "filled_path_count": len(selected),
            "primitive_counts": dict(Counter(item["kind"] for r in selected for item in r["items"])),
            "bbox": [min(b[0] for b in boxes), min(b[1] for b in boxes),
                     max(b[2] for b in boxes), max(b[3] for b in boxes)],
            "diagnostic_render": image_path.name,
            "diagnostic_render_sha256": digest(image_path),
            "colour_is_not_an_automatically_verified_component_identity": True,
        })
        layer.close()

    combined = fitz.open()
    target = combined.new_page(width=page.rect.width, height=page.rect.height)
    paint(target, records)
    target.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False).save(output / "vectors_only.png")
    page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False).save(output / "original_page.png")
    report = {
        "schema": "motoo.redraw-geometry-inventory.v1",
        "status": "GEOMETRY_INVENTORIED_PD_TRANSCRIPTION_PENDING",
        "source_file": str(source), "source_sha256": digest(source),
        "source_bytes": source.stat().st_size,
        "pages": len(document), "page_rect": list(page.rect),
        "pdf_metadata": document.metadata,
        "pymupdf_version": importlib.metadata.version("PyMuPDF"),
        "filled_paths": len(records), "text_extraction": page.get_text(),
        "raster_images": [{"xref": i[0], "mask_xref": i[1], "width": i[2], "height": i[3],
                           "bits_per_component": i[4], "colour_space": i[5]}
                          for i in page.get_images(full=True)],
        "colour_layers": rows,
        "geometry_file": "filled_paths.json.gz",
        "geometry_sha256": digest(output / "filled_paths.json.gz"),
        "original_page_render_sha256": digest(output / "original_page.png"),
        "vectors_only_render_sha256": digest(output / "vectors_only.png"),
        "centre_lines_extracted": False, "crossings_assigned": False, "pd_generated": False,
        "isotopy_certified": False, "katie_position_certified": False,
        "limitations": [
            "Filled pen outlines are not centre-line link arcs.",
            "Colour layers also contain labels, dots and possible incidental marks.",
            "Erasures split some pen outlines into small fragments.",
            "No endpoint pairing or over/under convention is inferred by this tool.",
            "Diagnostic vector-only renders omit the embedded raster background.",
        ],
    }
    (output / "GEOMETRY_INVENTORY.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "source_sha256": report["source_sha256"],
                      "filled_paths": len(records),
                      "colour_layers": [{"name": r["colour_name"], "paths": r["filled_path_count"]} for r in rows],
                      "output": str(output)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scale", type=float, default=2.0)
    args = parser.parse_args()
    require(0 < args.scale <= 8, "Render scale must be positive and at most eight")
    inspect(args.source, args.out, args.scale)

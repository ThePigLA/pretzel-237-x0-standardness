"""Shared input-integrity checks for the public construction entry points."""
from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import json
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MASTER_SHA256 = "251d7ad38235869fa268e9ee6ec9da2f397cf1ca71db267d0bf0ca9b434c5e4c"
FACETS_SHA256 = "f40b43c09373dff94cde542a3076a8c2f482b1632726b63eed2c5610ddd1e92b"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def canonical_sha256(value) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def verify_sources(root: Path = ROOT) -> dict:
    """Check every unmodified imported file and the exact accepted master."""
    root = root.resolve()
    inventory_path = root / "provenance" / "source_files.json"
    inventory = read_json(inventory_path)
    if inventory.get("schema") != "construction-source-inventory-v1":
        raise ValueError("Unrecognized construction input inventory")
    rows = inventory["files"]
    seen = set()
    total_bytes = 0
    for row in rows:
        relative = row["path"]
        if relative in seen:
            raise ValueError("Duplicate inventory path: " + relative)
        seen.add(relative)
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or path == root:
            raise ValueError("Input path escapes the construction directory")
        if not path.is_file():
            raise ValueError("Missing construction input: " + relative)
        if path.stat().st_size != row["bytes"]:
            raise ValueError("Construction input size changed: " + relative)
        if file_sha256(path) != row["sha256"]:
            raise ValueError("Construction input hash changed: " + relative)
        total_bytes += row["bytes"]
    master = root / "reference" / "x0" / "marked_candidate.json.gz"
    if file_sha256(master) != MASTER_SHA256:
        raise ValueError("The accepted X0 reference master has changed")
    model = read_json(master)
    if canonical_sha256(model["facets"]) != FACETS_SHA256:
        raise ValueError("The reference facet list has changed")
    return {
        "status": "PASS",
        "checked_files": len(rows),
        "checked_bytes": total_bytes,
        "inventory_sha256": file_sha256(inventory_path),
        "reference_master_file_sha256": MASTER_SHA256,
        "reference_facets_canonical_sha256": FACETS_SHA256,
        "meaning": "Integrity of supplied inputs; this alone does not rerun their mathematical checks.",
    }


def environment_record() -> dict:
    packages = {}
    for name in ("regina", "snappy", "spherogram", "sympy", "psutil"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": packages,
        "python_optimization": sys.flags.optimize,
    }


if __name__ == "__main__":
    print(json.dumps(verify_sources(), indent=2))

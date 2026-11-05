#!/usr/bin/env python3
"""Shared integrity and logging helpers for the Motoo Z repository module."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parent
LEGACY = ROOT / "legacy/code"
REFERENCE = ROOT / "reference/published_fig2_S4"


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def check_integrity():
    manifest_path = ROOT / "MANIFEST_SHA256.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get("schema") == "motoo-z-files-v1", "Unknown manifest schema")
    entries = manifest["files"]
    total = 0
    for name, expected in entries.items():
        relative = Path(name)
        require(not relative.is_absolute() and ".." not in relative.parts,
                "Unsafe manifest path: " + name)
        path = ROOT / relative
        require(path.is_file(), "Missing released file: " + name)
        require(not path.is_symlink() and path.resolve().is_relative_to(ROOT),
                "Released file resolves outside module: " + name)
        require(sha256(path) == expected["sha256"], "Changed released file: " + name)
        require(path.stat().st_size == expected["bytes"], "Changed file size: " + name)
        total += expected["bytes"]
    return {"status": "PASS_RELEASED_FILE_INTEGRITY", "files": len(entries),
            "bytes": total, "manifest_sha256": sha256(manifest_path)}


def fresh_output(path):
    result = Path(path).expanduser().resolve()
    require(not result.exists(), "Output directory already exists: " + str(result))
    for protected in (ROOT / "legacy", ROOT / "reference", ROOT / "certificates",
                      ROOT / "provenance", ROOT / "verification", ROOT / "environment"):
        require(not result.is_relative_to(protected),
                "Choose an output directory outside the released inputs: " + str(result))
    require(not ROOT.is_relative_to(result), "Output would contain the module")
    result.mkdir(parents=True)
    return result


def run_stage(name, argv, out, report):
    """Run a supplied checker; retain its complete output and exit code."""
    out = Path(out)
    log = out / "logs" / (name + ".log")
    log.parent.mkdir(parents=True, exist_ok=True)
    command = [os.fspath(arg) for arg in argv]
    started = time.monotonic()
    print("Running " + name, flush=True)
    with log.open("w") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                cwd=out, check=False)
    record = {"name": name, "argv": command, "returncode": result.returncode,
              "seconds": time.monotonic() - started,
              "log": str(log.relative_to(out)), "log_sha256": sha256(log)}
    report.setdefault("stages", []).append(record)
    write_json(out / "STAGES.json", report)
    require(result.returncode == 0,
            "Stage failed: " + name + "; inspect " + str(log))
    return record


if __name__ == "__main__":
    print(json.dumps(check_integrity(), indent=2))

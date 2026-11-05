#!/usr/bin/env python3
"""Check every distributed file against the SHA256SUMS manifest."""
import hashlib
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = root / "SHA256SUMS"
    count = 0
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError("Manifest path escapes package")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError("SHA256 mismatch: " + relative)
        count += 1
    print(f"PASS: SHA256 checks for {count} distributed files")


if __name__ == "__main__":
    main()

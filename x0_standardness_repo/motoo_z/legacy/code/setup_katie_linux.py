#!/usr/bin/env python3
"""Extract pinned official Ubuntu 24.04 Katie packages locally, without apt/sudo.

For Ubuntu 24.04 x86_64 with the usual libc/libgmp/libxml/zlib dependencies.
On other platforms install Regina 7.4.1 and pass its Katie path to run_katie.py.
The package bytes are verified against the recorded SHA256 values before use.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import urllib.request

PACKAGES = [
    ("https://people.debian.org/~bab/regina/pool/main/r/regina-normal/"
     "regina-normal_7.4.1-1~ubuntu24.04.1~regina1_amd64.deb",
     "5e51f280a41437330cf37086e75b999cf7f8999b7db0f78884ea98502b2f0a6b"),
    ("https://archive.ubuntu.com/ubuntu/pool/main/g/gmp/"
     "libgmpxx4ldbl_6.3.0+dfsg-2ubuntu6.1_amd64.deb",
     "6f59344240b6dc139ed23ef236b8aa146e5d25e23aa90e79f814e2e4cc4b5752"),
    ("https://archive.ubuntu.com/ubuntu/pool/universe/t/tokyocabinet/"
     "libtokyocabinet9t64_1.4.48-15.1build3_amd64.deb",
     "21883c591e7f619eb6ae0acceff9d021ad35f6093daf5f6344f7bcecb5daecc1"),
]


def main():
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "amd64"):
        raise RuntimeError("These pinned packages target Linux x86_64 only")
    root = Path(__file__).resolve().parents[1] / "scratch/katie"
    downloads = root / "downloads"
    downloads.mkdir(parents=True, exist_ok=True)

    def fetch(item):
        url, expected = item
        path = downloads / url.rsplit("/", 1)[1]
        if not path.exists():
            with urllib.request.urlopen(url, timeout=60) as response:
                path.write_bytes(response.read())
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError("Checksum mismatch: " + str(path))
        return {"url": url, "file": str(path), "sha256": actual,
                "bytes": path.stat().st_size}

    receipts = list(ThreadPoolExecutor(3).map(fetch, PACKAGES))
    for item in receipts:
        subprocess.run(["dpkg-deb", "-x", item["file"], str(root / "prefix")], check=True)
    executable = root / "prefix/usr/bin/katie"
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = str(root / "prefix/usr/lib/x86_64-linux-gnu")
    version = subprocess.run([str(executable), "--version"], env=env, text=True,
                             capture_output=True, check=True).stdout.strip()
    if "7.4.1" not in version:
        raise RuntimeError("Unexpected Katie version: " + version)
    result = {"packages": receipts, "version": version,
              "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
              "system_packages_changed": False}
    (root / "SETUP_RECEIPT.json").write_text(json.dumps(result, indent=2) + "\n")
    print(version)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify the exact released files before running any receipt-refreshing check."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST_SHA256.json').read_text())
for name,row in manifest['files'].items():
    b=(root/name).read_bytes()
    if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:
        raise ValueError('Release file mismatch: '+name)
print('PASS',len(manifest['files']),'released files match their exact SHA-256 and byte counts.')

#!/usr/bin/env python3
"""Replay the saved Figure 2 handlebody paths from portable bundle paths.

This verifies partial reductions, not standardness of the closed manifold.
The closed-completion checks are separate and are described in the README.
Every invocation uses a new output directory. No original receipt is replaced.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def verify(root, out, manifest_only=False):
    root, out = Path(root).resolve(), Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    manifest = root / 'MANIFEST.json'
    manifest_receipt = None
    if manifest.exists():
        listing = json.loads(manifest.read_text())
        for item in listing['files']:
            candidate = (root / item['path']).resolve()
            require(candidate.is_relative_to(root), 'unsafe path in bundle manifest')
            require(candidate.is_file(), 'missing bundled file: ' + item['path'])
            require(candidate.stat().st_size == item['bytes'] and
                    file_sha(candidate) == item['sha256'],
                    'bundle checksum mismatch: ' + item['path'])
        manifest_receipt = {'status': 'PASS_EVERY_MANIFEST_FILE',
                            'files': len(listing['files']),
                            'manifest_sha256': file_sha(manifest)}
        print(json.dumps(manifest_receipt), flush=True)
    else:
        require(not manifest_only, 'MANIFEST.json is present only in the frozen bundle')
    if manifest_only:
        write_json(out / 'VERIFICATION.json', manifest_receipt)
        return manifest_receipt

    require(importlib.metadata.version('regina') == '7.4.1',
            'Install regina==7.4.1 before replaying these paths.')

    base = Path('results/figure2_working_20261102')
    raw = Path('results/redraw_20261102/katie_figure2_dashed_v1/raw_katie.rga')
    staircase = base / 'completion/raw_staircase'
    ideal = base / 'fast_ideal/retrace_270'
    finite = base / 'native_diagnosis/finite1691_verified/from_original'
    candidate = base / 'native_diagnosis/candidate2136_verified/from_original'
    jobs = [
        ('prism_construction', ['code/verify_staircase_truncation.py',
         '--source', str(raw), '--folder', str(staircase)]),
        ('ideal_270', ['code/audit_figure2_snaps.py', str(raw),
         str(ideal / 'moves.jsonl')]),
        ('finite_1691', ['code/reduce_finite_checked.py', 'replay',
         str(finite), '--source', str(staircase / 'finite_staircase.rga')]),
        ('completion_branch_2136', ['code/reduce_finite_checked.py', 'replay',
         str(candidate), '--source', str(staircase / 'finite_staircase.rga')]),
    ]
    results = []
    for name, command in jobs:
        receipt = out / (name + '.json')
        args = [sys.executable, '-O', *command, '--out', str(receipt)]
        print(json.dumps({'stage': name, 'status': 'RUNNING'}), flush=True)
        then = time.monotonic()
        completed = subprocess.run(args, cwd=root, capture_output=True, text=True)
        # Publish a closed file, rather than holding an append writer open
        # across the computation and its checkpoints.
        log = out / (name + '.log')
        log.write_text(completed.stdout + completed.stderr)
        require(completed.returncode == 0,
                name + ' failed; inspect ' + str(log))
        require(receipt.is_file(), name + ' returned without a verification receipt')
        evidence = json.loads(receipt.read_text())
        require(str(evidence.get('status', '')).startswith('PASS'),
                name + ' did not report PASS')
        record = {'stage': name, 'status': evidence['status'],
                  'receipt': receipt.name, 'receipt_sha256': file_sha(receipt),
                  'seconds': time.monotonic() - then}
        results.append(record)
        print(json.dumps(record), flush=True)

    # Bind the ideal replay to the archived endpoint, not just its advertised size.
    import regina
    ideal_receipt = json.loads((out / 'ideal_270.json').read_text())
    endpoint = regina.open(str(root / ideal / 'endpoint_270_974b67f26bf8.rga'))
    require(isinstance(endpoint, regina.Triangulation4), 'missing ideal endpoint')
    require(hashlib.sha256(endpoint.tightEncoding().encode()).hexdigest() ==
            ideal_receipt['final_tight_sha256'], 'ideal endpoint differs from replay')
    require(endpoint.size() == ideal_receipt['final_pentachora'] == 270,
            'wrong ideal endpoint size')
    finite_bindings = []
    for name, endpoint_relative, count in [
        ('finite_1691', base / 'native_diagnosis/finite1691_verified/lower_reference.rga', 1691),
        ('completion_branch_2136', base / 'completion/finite_candidate_phase3_a/lower_candidate.rga', 2136),
    ]:
        evidence = json.loads((out / (name + '.json')).read_text())
        endpoint_path = root / endpoint_relative
        selected = regina.open(str(endpoint_path))
        require(isinstance(selected, regina.Triangulation4), 'missing selected finite endpoint')
        selected_hash = hashlib.sha256(selected.tightEncoding().encode()).hexdigest()
        require(selected.size() == count and selected_hash == evidence['end_tight_sha256'],
                name + ': selected endpoint differs from replay')
        finite_bindings.append({'stage': name, 'file': str(endpoint_relative),
                                'file_sha256': file_sha(endpoint_path),
                                'tight_sha256': selected_hash, 'pentachora': count})
    final = {'status': 'PASS_SAVED_PARTIAL_HANDLEBODY_PATHS',
             'working_assumption': 'Figure 2 red passes over dashed green',
             'manifest_verification': manifest_receipt, 'checks': results,
             'ideal_endpoint_file_sha256': file_sha(root / ideal / 'endpoint_270_974b67f26bf8.rga'),
             'selected_finite_endpoint_bindings': finite_bindings,
             'closed_completion_verified_by_this_script': False,
             'standardness_asserted': False,
             'seconds': time.monotonic() - started}
    write_json(out / 'VERIFICATION.json', final)
    print(json.dumps(final), flush=True)
    return final


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--manifest-only', action='store_true')
    args = parser.parse_args()
    verify(args.root, args.out, args.manifest_only)

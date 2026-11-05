#!/usr/bin/env python3
"""Reproduce the Figure 2 reductions and closed PL construction.

No standard-sphere recognition is asserted by this command. The independent
closed-search path, if supplied, has its own replay command and success test.
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


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(root, out, stage='all'):
    root, out = Path(root).resolve(), Path(out).resolve()
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina 7.4.1')
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    checks = []

    def run(name, command, receipt, output_directory=False):
        receipt = Path(receipt)
        destination = receipt.parent if output_directory else receipt
        args = [sys.executable, '-O', *map(str, command), '--out', str(destination)]
        print(json.dumps({'stage': name, 'status': 'RUNNING'}), flush=True)
        then = time.monotonic()
        completed = subprocess.run(args, cwd=root, capture_output=True, text=True)
        log = out / (name + '.log')
        log.write_text(completed.stdout + completed.stderr)
        require(completed.returncode == 0 and receipt.is_file(), name + ' failed; inspect ' + str(log))
        report = json.loads(receipt.read_text())
        require(str(report.get('status', '')).startswith('PASS'), name + ' did not report PASS')
        result = {'stage': name, 'status': report['status'],
                  'receipt': str(receipt.relative_to(out)), 'receipt_sha256': sha(receipt),
                  'seconds': time.monotonic() - then}
        checks.append(result)
        print(json.dumps(result), flush=True)
        return report

    base = Path('results/figure2_working_20261102')
    branch = base / 'completion/finite_candidate_phase3_a'
    staircase = base / 'completion/raw_staircase'
    raw = Path('results/redraw_20261102/katie_figure2_dashed_v1/raw_katie.rga')
    upper9 = Path('results/boundary/upper_cap_selftest_v1/reference_cap.rga')
    upper474 = Path('results/boundary/standard_cap_six_boundary_v1/cap_with_matched_boundary.rga')
    upper483 = branch / 'matched_upper_cap/layered_piece.rga'
    lower = branch / 'lower_candidate.rga'
    layered = branch / 'layered_lower_fast/layered_piece.rga'
    closed = branch / 'closed_fast/closed_candidate.rga'
    if stage == 'all':
        run('saved_handlebody_paths', ['code/verify_saved_figure2_paths.py'],
            out / 'paths/VERIFICATION.json', True)
        prefix_receipt = out / 'paths/completion_branch_2136.json'
    else:
        if (root / 'MANIFEST.json').exists():
            run('archive_integrity', ['code/verify_saved_figure2_paths.py', '--manifest-only'],
                out / 'integrity/VERIFICATION.json', True)
        run('prism_construction', ['code/verify_staircase_truncation.py', '--source', raw,
            '--folder', staircase], out / 'prism.json')
        prefix_receipt = out / 'completion_prefix.json'
        run('completion_prefix', ['code/reduce_finite_checked.py', 'replay',
            base / 'native_diagnosis/candidate2136_verified/from_original',
            '--source', staircase / 'finite_staircase.rga'], prefix_receipt)

    run('upper_15_collapses', ['code/verify_upper_cap_guarded.py',
        '--certificate', base / 'upper_cap_guarded_v1', '--target', upper9], out / 'upper_collapses.json')
    run('upper_474_shellings', ['code/verify_upper_cap_shellings.py', '--upper', upper483,
        '--intermediate', upper474, '--reference', upper9], out / 'upper_shellings.json')
    finite_receipt = out / 'finite_admission/FINITE_MANIFOLD_ADMISSION.json'
    run('finite_manifold_admission', ['code/audit_raw_algebra4.py', layered,
        '--method', 'primal', '--finite'], finite_receipt, True)
    run('lower_22581_shellings', ['code/verify_layer_shellings_raw.py', '--source', lower,
        '--layered', layered, '--admission', finite_receipt], out / 'lower_shellings.json')
    run('independent_collar_attachments', ['code/audit_layering_collar.py', '--source', lower,
        '--path', branch / 'short_boundary_path.json', '--layered', layered,
        '--events', branch / 'layered_lower_fast/layering_events.jsonl.gz'], out / 'collar.json')
    run('exact_closed_gluing', ['code/audit_fast_closed_gluing.py', '--folder', branch / 'closed_fast'],
        out / 'closed_gluing.json')
    closed_receipt = out / 'closed_admission/CLOSED_MANIFOLD_ADMISSION.json'
    closed_report = run('closed_manifold_admission', ['code/audit_raw_algebra4.py', closed,
        '--method', 'primal'], closed_receipt, True)
    run('archived_construction_bindings', ['code/bind_collar_completion.py', '--root', root],
        out / 'construction_bindings.json')

    import regina
    prefix = json.loads(prefix_receipt.read_text())
    selected = regina.open(str(root / lower))
    require(isinstance(selected, regina.Triangulation4) and selected.size() == 2136 and
            hashlib.sha256(selected.tightEncoding().encode()).hexdigest() == prefix['end_tight_sha256'],
            'finite replay is not bound to the selected completion source')
    require(closed_report['source_sha256'] == sha(root / closed), 'closed input changed')
    result = {'status': 'PASS_REPRODUCED_FIGURE2_PL_CONSTRUCTION', 'stage_selection': stage,
              'working_assumption': 'Figure 2 red passes over dashed green; author confirmation pending',
              'checks': checks, 'closed_pentachora': 25200,
              'closed_sha256': sha(root / closed), 'closed_f_vector': closed_report['f_vector'],
              'closed_vertex_links_verified': len(closed_report['vertex_links']),
              'standardness_asserted': False,
              'seconds': time.monotonic() - started}
    destination = out / 'VERIFICATION.json'
    temporary = destination.with_suffix('.tmp')
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    temporary.replace(destination)
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--stage', choices=['all', 'construction'], default='all')
    args = parser.parse_args()
    verify(args.root, args.out, args.stage)

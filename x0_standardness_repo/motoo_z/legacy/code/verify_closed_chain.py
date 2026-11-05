#!/usr/bin/env python3
"""Replay successive closed-4D move certificates with exact join bindings.

Run from the extracted project root, or pass --root. The manifest uses paths
relative to that root. This imports only the separate verifier, not a search
implementation. It can itself be wrapped by run_with_primal_recognition.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from verify_closed_bursts import digest, load, require, state_digest, verify


def contained(root, relative):
    path = Path(relative)
    require(not path.is_absolute() and '..' not in path.parts,
            'chain paths must be contained project-relative names')
    result = root / path
    require(result.resolve().is_relative_to(root.resolve()), 'path escapes project root')
    return result


def run(manifest, root, output):
    manifest, root, output = Path(manifest).resolve(), Path(root).resolve(), Path(output).resolve()
    require(not output.exists(), 'chain output directory already exists')
    manifest_bytes = manifest.read_bytes()
    manifest_sha = hashlib.sha256(manifest_bytes).hexdigest()
    code_sha = digest(__file__)
    checker_path = Path(__file__).with_name('verify_closed_bursts.py')
    checker_sha = digest(checker_path)
    specification = json.loads(manifest_bytes)
    require(specification['schema'] == 'motoo.closed-move-chain.v1', 'unknown chain schema')
    require(specification['segments'], 'empty chain')
    output.mkdir(parents=True)
    started = time.monotonic()
    preceding_state = specification['initial_tight_sha256']
    total_moves = 0
    reports = []
    for number, segment in enumerate(specification['segments'], 1):
        directory = contained(root, segment['certificate'])
        source = contained(root, segment['source'])
        require(digest(source) == segment['source_sha256'], 'segment source file differs')
        source_state = state_digest(load(source))
        require(source_state == preceding_state == segment['start_tight_sha256'],
                'successive segments do not meet at the exact labelled triangulation')
        if number == 1:
            require(segment['source_sha256'] == specification['initial_source_sha256'],
                    'first source is not the claimed original input')
        report_file = output / ('segment_%02d.json' % number)
        result = verify(directory, source, report_file)
        require(result['source_sha256'] == segment['source_sha256'],
                'independently replayed source differs from the chain manifest')
        require(result['moves'] == segment['moves'], 'segment move count differs')
        require(result['final_tight_sha256'] == segment['end_tight_sha256'],
                'segment endpoint differs')
        require(result['final_admission']['pentachora'] == segment['end_pentachora'],
                'segment endpoint size differs')
        total_moves += result['moves']
        preceding_state = result['final_tight_sha256']
        reports.append({'certificate': segment['certificate'],
                        'source': segment['source'],
                        'result': str(report_file.relative_to(output)),
                        'result_sha256': digest(report_file),
                        'moves': result['moves'],
                        'end_pentachora': result['final_admission']['pentachora'],
                        'end_tight_sha256': preceding_state,
                        'standard_sphere_proved': result['standard_sphere_proved']})
    standard = reports[-1]['standard_sphere_proved']
    require(total_moves == specification['total_moves'], 'chain total move count differs')
    require(standard == specification['standard_sphere_proved'], 'chain standardness claim differs')
    require(digest(manifest) == manifest_sha, 'chain manifest changed during replay')
    require(digest(__file__) == code_sha and digest(checker_path) == checker_sha,
            'a checker file changed during replay')
    result = {'schema': 'motoo.closed-move-chain.replay.v1',
              'status': 'PASS_STANDARDIZATION' if standard else 'PASS_PARTIAL_CLOSED_CHAIN',
              'manifest_sha256': manifest_sha, 'segments': reports,
              'total_moves': total_moves, 'final_tight_sha256': preceding_state,
              'final_pentachora': reports[-1]['end_pentachora'],
              'standard_sphere_proved': standard,
              'all_segment_joins_checked_by_exact_labelled_state': True,
              'seconds': time.monotonic() - started,
              'verifier_sha256': code_sha, 'segment_verifier_sha256': checker_sha}
    temporary = output / 'CHAIN_REPLAY.json.tmp'
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    temporary.replace(output / 'CHAIN_REPLAY.json')
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--out', type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.manifest, arguments.root, arguments.out)

#!/usr/bin/env python3
"""Extract the exact logged prefix ending at a closed search's best state.

This is certificate curation, not move verification. Run the separate
verify_closed_bursts.py checker on the resulting directory before relying
on it. The original search and its complete history are left untouched.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil

from verify_closed_bursts import digest, is_standard, load, local_path, require, state_digest


def extract(directory, output):
    directory, output = Path(directory).resolve(), Path(output).resolve()
    require(not output.exists(), 'output already exists')
    source_bytes = (directory / 'SOURCE.json').read_bytes()
    check_bytes = (directory / 'CHECKPOINT.json').read_bytes()
    source, check = json.loads(source_bytes), json.loads(check_bytes)
    require(check['status'] != 'RUNNING', 'finish or stop the source run first')
    limit = check['best_moves']
    require(type(limit) is int and 0 <= limit <= check['moves'], 'invalid best index')
    best = local_path(directory, check['best_file'])
    require(digest(best) == check['best_rga_sha256'], 'best RGA checksum differs')
    tri = load(best)
    require(state_digest(tri) == check['best_tight_sha256'] and
            tri.size() == check['best_pentachora'], 'best RGA state differs')
    output.mkdir(parents=True)
    (output / 'moves').mkdir()
    (output / 'states').mkdir()
    (output / 'SOURCE.json').write_bytes(source_bytes)
    shutil.copy2(best, output / 'states/best.rga')
    require(digest(output / 'states/best.rga') == check['best_rga_sha256'],
            'copied best RGA checksum differs')
    n, chunks, counts = 0, [], Counter()
    current = source['start_tight_sha256']
    size = source['source_admission']['pentachora']
    maximum = size
    cycle = 0
    bound_chunks = {}
    for chunk in check['chunks']:
        if n == limit:
            break
        path = local_path(directory, chunk['file'])
        data = path.read_bytes()
        require(hashlib.sha256(data).hexdigest() == chunk['sha256'], 'source chunk differs')
        require(chunk['first_index'] == n, 'nonconsecutive source chunk')
        lines = data.splitlines(keepends=True)
        require(len(lines) == chunk['moves'], 'source chunk length differs')
        selected = lines[:min(len(lines), limit - n)]
        first = n
        for line in selected:
            row = json.loads(line)
            require(row['index'] == n and row['before_sha256'] == current and
                    row['before_pentachora'] == size, 'broken source prefix')
            current, size = row['after_sha256'], row['pentachora']
            maximum = max(maximum, size)
            counts[str(row['operation'][:2])] += 1
            cycle = row['cycle']
            n += 1
        destination = output / 'moves' / ('%06d.jsonl' % len(chunks))
        destination.write_bytes(b''.join(selected))
        chunks.append({'file': str(destination.relative_to(output)),
                       'sha256': digest(destination), 'first_index': first,
                       'moves': len(selected)})
        bound_chunks[chunk['file']] = chunk['sha256']
    require(n == limit and current == check['best_tight_sha256'] and
            size == check['best_pentachora'], 'prefix does not end at the best state')
    result = dict(check)
    result.update({'status': 'EXTRACTED_PREFIX_TO_BEST_REPLAY_REQUIRED',
                   'moves': n, 'current_pentachora': size,
                   'current_tight_sha256': current, 'current_file': 'states/best.rga',
                   'current_rga_sha256': digest(output / 'states/best.rga'),
                   'best_file': 'states/best.rga', 'chunks': chunks,
                   'counts': dict(counts), 'max_pentachora': maximum, 'cycle': cycle,
                   'standard_endpoint': is_standard(tri), 'replay_verified': False,
                   'elapsed_seconds_describes': 'original completed search, not extracted prefix replay'})
    (output / 'CHECKPOINT.json').write_text(json.dumps(result, indent=2) + '\n')
    require(digest(directory / 'SOURCE.json') == hashlib.sha256(source_bytes).hexdigest()
            and digest(directory / 'CHECKPOINT.json') == hashlib.sha256(check_bytes).hexdigest()
            and digest(best) == check['best_rga_sha256']
            and all(digest(local_path(directory, name)) == expected
                    for name, expected in bound_chunks.items()), 'source changed during extraction')
    receipt = {'status': 'PASS_EXACT_PREFIX_EXTRACTION_NOT_A_MOVE_REPLAY',
               'source_metadata_sha256': hashlib.sha256(source_bytes).hexdigest(),
               'original_checkpoint_sha256': hashlib.sha256(check_bytes).hexdigest(),
               'original_run_moves': check['moves'], 'extracted_moves': n,
               'best_pentachora': size, 'best_tight_sha256': current,
               'best_rga_sha256': digest(output / 'states/best.rga'),
               'extracted_checkpoint_sha256': digest(output / 'CHECKPOINT.json'),
               'extractor_sha256': digest(__file__),
               'independent_move_replay_required': True}
    (output / 'EXTRACTION.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    extract(args.directory, args.out)

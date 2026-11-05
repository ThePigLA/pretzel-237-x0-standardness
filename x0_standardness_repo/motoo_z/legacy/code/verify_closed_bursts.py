#!/usr/bin/env python3
"""Independently replay a guarded closed-4D path without importing its runner.

Requires Regina 7.4.1.  Every move is guarded again, every labelled state is
hashed, and both the minimum checkpoint and the final checkpoint are bound.
An optional audited recognition launcher can accelerate source admission;
its presence is explicitly recorded.  No random search or simplify is used.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import time

import regina


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def state_digest(tri):
    return hashlib.sha256(tri.tightEncoding().encode('ascii')).hexdigest()


def load(path):
    packet = regina.open(str(path))
    require(isinstance(packet, regina.Triangulation4), 'not a 4D triangulation')
    return regina.Triangulation4(packet)


def local_path(directory, relative):
    name = Path(relative)
    require(not name.is_absolute() and '..' not in name.parts,
            'certificate file path must be relative and contained')
    return directory / name


def admit(tri):
    require(tri.size() > 0 and tri.isValid() and tri.isConnected() and
            tri.isOrientable() and not tri.isIdeal() and tri.isClosed(),
            'source is not admitted as a finite closed orientable PL manifold')
    return {'pentachora': tri.size(), 'vertices': tri.countVertices(),
            'finite_closed_connected_orientable_valid': True}


def apply_checked(tri, operation):
    require(isinstance(operation, list) and len(operation) == 3,
            'malformed move')
    kind, dimension, index = operation
    require(type(dimension) is int and 0 <= dimension <= 4,
            'invalid face dimension')
    total = tri.size() if dimension == 4 else tri.countFaces(dimension)
    require(type(index) is int and 0 <= index < total, 'invalid face index')
    face = tri.simplex(index) if dimension == 4 else tri.face(dimension, index)
    before, vertices = tri.size(), tri.countVertices()
    if kind == 'P':
        require(tri.hasPachner(face), 'Pachner guard rejects move')
        require(tri.pachner(face), 'Pachner application failed')
        require(tri.size() == before + 2 * dimension - 4, 'wrong Pachner size')
        if dimension in (0, 4):
            require(tri.countVertices() == vertices + (1 if dimension == 4 else -1),
                    'wrong vertex change for extreme Pachner move')
    elif kind == 'Z':
        require(dimension in (0, 1, 2), 'invalid 2-0 face dimension')
        require(tri.has20(face), '2-0 guard rejects move')
        require(tri.move20(face) and tri.size() == before - 2, '2-0 application failed')
    elif kind == 'C':
        require(dimension == 1 and face.vertex(0) != face.vertex(1),
                'collapse does not use an edge with distinct endpoints')
        require(tri.hasCollapseEdge(face), 'edge-collapse guard rejects move')
        require(tri.collapseEdge(face), 'edge-collapse application failed')
        require(tri.size() < before and tri.countVertices() == vertices - 1,
                'wrong collapse size or vertex change')
    elif kind == 'F':
        require(dimension == 1 and tri.has44(face), '4-4 guard rejects move')
        require(tri.move44(face) and tri.size() == before, '4-4 application failed')
    else:
        raise ValueError('unsupported move code')


def is_standard(tri):
    if tri.size() != 2:
        return False
    standard = regina.Triangulation4()
    first, second = standard.newSimplices(2)
    for facet in range(5):
        first.join(facet, second, regina.Perm5())
    return tri.isIsomorphicTo(standard) is not None


def verify(directory, source, output):
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina 7.4.1')
    directory, source, output = Path(directory), Path(source), Path(output)
    require(not output.exists(), 'output receipt already exists')
    metadata_bytes = (directory / 'SOURCE.json').read_bytes()
    checkpoint_bytes = (directory / 'CHECKPOINT.json').read_bytes()
    meta = json.loads(metadata_bytes)
    checkpoint = json.loads(checkpoint_bytes)
    pinned = {str(directory / 'SOURCE.json'): hashlib.sha256(metadata_bytes).hexdigest(),
              str(directory / 'CHECKPOINT.json'): hashlib.sha256(checkpoint_bytes).hexdigest(),
              str(source): digest(source), str(Path(__file__)): digest(__file__)}
    require(pinned[str(source)] == meta['source_sha256'], 'source file checksum differs')
    tri = load(source)
    started = time.monotonic()
    admission = admit(tri)
    current = state_digest(tri)
    require(current == meta['start_tight_sha256'], 'initial labelled state differs')
    minimum, maximum, moves = tri.size(), tri.size(), 0
    counts = Counter()
    minimum_state = current
    minimum_index = 0
    for chunk in checkpoint['chunks']:
        path = local_path(directory, chunk['file'])
        chunk_bytes = path.read_bytes()
        require(hashlib.sha256(chunk_bytes).hexdigest() == chunk['sha256'],
                'move chunk checksum differs')
        pinned[str(path)] = chunk['sha256']
        require(chunk['first_index'] == moves, 'nonconsecutive chunk')
        local_count = 0
        for line in chunk_bytes.decode('utf-8').splitlines():
            row = json.loads(line)
            require(row['index'] == moves and row['before_sha256'] == current,
                    'broken labelled-state chain')
            require(row['before_pentachora'] == tri.size(), 'wrong pre-move size')
            apply_checked(tri, row['operation'])
            current = state_digest(tri)
            require(row['after_sha256'] == current and row['pentachora'] == tri.size(),
                    'post-move state differs')
            moves += 1
            local_count += 1
            counts[str(row['operation'][:2])] += 1
            maximum = max(maximum, tri.size())
            if tri.size() < minimum:
                minimum = tri.size()
                minimum_state = current
                minimum_index = moves
        require(local_count == chunk['moves'], 'chunk move count differs')
    require(moves == checkpoint['moves'] and dict(counts) == checkpoint['counts'],
            'total move count differs')
    require(current == checkpoint['current_tight_sha256'] and
            tri.size() == checkpoint['current_pentachora'], 'wrong final state')
    require(minimum == checkpoint['best_pentachora'] and
            minimum_state == checkpoint['best_tight_sha256'] and
            minimum_index == checkpoint['best_moves'], 'wrong best-state claim')
    require(maximum == checkpoint['max_pentachora'], 'wrong maximum size')
    for prefix, expected_state in [('current', current), ('best', minimum_state)]:
        path = local_path(directory, checkpoint[prefix + '_file'])
        pinned[str(path)] = checkpoint[prefix + '_rga_sha256']
        require(digest(path) == pinned[str(path)],
                prefix + ' RGA checksum differs')
        require(state_digest(load(path)) == expected_state,
                prefix + ' RGA labelled state differs')
    final_admission = admit(tri)
    standard = is_standard(tri)
    require(standard == checkpoint['standard_endpoint'], 'standardness flag differs')
    require(all(digest(path) == expected for path, expected in pinned.items()),
            'a source, certificate, state, or verifier file changed during replay')
    report = {'schema': 'motoo.closed-bursts.independent-replay.v1',
              'status': 'PASS_STANDARDIZATION' if standard else 'PASS_PARTIAL_CLOSED_PATH',
              'source_sha256': pinned[str(source)], 'moves': moves, 'counts': dict(counts),
              'source_admission': admission, 'final_admission': final_admission,
              'best_pentachora': minimum, 'best_at_move': minimum_index,
              'maximum_pentachora': maximum, 'final_tight_sha256': current,
              'standard_sphere_proved': standard,
              'source_metadata_sha256': pinned[str(directory / 'SOURCE.json')],
              'checkpoint_sha256': pinned[str(directory / 'CHECKPOINT.json')],
              'verifier_sha256': pinned[str(Path(__file__))],
              'all_bound_files_unchanged_through_replay': True,
              'seconds_including_source_admission': time.monotonic() - started,
              'recognition_runtime': os.environ.get('MOTOO_RECOGNITION_RUNTIME'),
              'trust': 'Regina 7.4.1 checked local moves; separately implemented replay; exact standard target',
              'runner_imported': False, 'native_4d_simplify_used': False}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(output)
    print(json.dumps(report), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    arguments = parser.parse_args()
    verify(arguments.directory, arguments.source, arguments.out)

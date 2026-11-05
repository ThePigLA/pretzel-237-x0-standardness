#!/usr/bin/env python3
"""Reduce a finite 4-manifold using streamed, guarded Regina moves.

The source may have real boundary or be closed.  Every move, including a
boundary shelling or an edge collapse, is replayable with its guard active.
Search is heuristic; failure has no implication for the PL type.  The only
standard-sphere success criterion is an isomorphism to the explicit double
of a 4-simplex.  No native simplify() is used in the four-dimensional path.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import random
import signal
import time
import regina


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(tri):
    return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()


def save_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def load4(path):
    packet = regina.open(str(path))
    require(isinstance(packet, regina.Triangulation4), 'source is not a 4D triangulation')
    return regina.Triangulation4(packet)


def save4(path, tri):
    path = Path(path)
    temporary = path.with_suffix('.tmp.rga')
    packet = regina.PacketOfTriangulation4(tri)
    packet.setLabel('Guarded finite 4D reduction checkpoint')
    require(packet.save(str(temporary)), 'could not save checkpoint')
    temporary.replace(path)


def standard_sphere():
    tri = regina.Triangulation4()
    a, b = tri.newSimplices(2)
    for facet in range(5):
        a.join(facet, b, regina.Perm5())
    return tri


def is_standard(tri):
    return tri.size() == 2 and tri.isIsomorphicTo(standard_sphere()) is not None


def admission(tri):
    require(tri.size() > 0 and tri.isValid() and tri.isConnected() and
            tri.isOrientable() and not tri.isIdeal(),
            'source must be a finite connected orientable PL 4-manifold')
    return {'pentachora': tri.size(), 'vertices': tri.countVertices(),
            'boundary_facets': tri.countBoundaryFacets(),
            'boundary_components': tri.countBoundaryComponents(),
            'closed': tri.isClosed(),
            'trust': 'Regina 7.4.1 validity and vertex-link recognition'}


def perform(tri, op):
    require(isinstance(op, list) and len(op) == 3, 'malformed move')
    kind, dimension, index = op
    require(type(dimension) is int and dimension in range(5), 'invalid face dimension')
    face_count = tri.size() if dimension == 4 else tri.countFaces(dimension)
    require(type(index) is int and 0 <= index < face_count,
            'invalid face index or dimension')
    face = tri.simplex(index) if dimension == 4 else tri.face(dimension, index)
    before = tri.size()
    if kind == 'P':
        require(tri.hasPachner(face) and tri.pachner(face), 'Pachner guard failed')
        require(tri.size() == before + 2 * dimension - 4, 'wrong Pachner size')
    elif kind == 'Z':
        require(dimension in (0, 1, 2), 'bad 2-0 dimension')
        require(tri.has20(face) and tri.move20(face), '2-0 guard failed')
        require(tri.size() == before - 2, 'wrong 2-0 size')
    elif kind == 'C':
        require(dimension == 1, 'collapse must concern an edge')
        vertices = tri.countVertices()
        require(tri.hasCollapseEdge(face) and tri.collapseEdge(face), 'collapse guard failed')
        require(tri.size() < before and tri.countVertices() == vertices - 1,
                'wrong collapse size/vertex change')
    elif kind == 'S':
        require(dimension == 4, 'shelling must concern a pentachoron')
        require(tri.hasShellBoundary(face) and tri.shellBoundary(face), 'shelling guard failed')
        require(tri.size() == before - 1, 'wrong shelling size')
    elif kind == 'F':
        require(dimension == 1, '4-4 must concern an edge')
        require(tri.has44(face) and tri.move44(face), '4-4 guard failed')
        require(tri.size() == before, 'wrong 4-4 size')
    else:
        raise ValueError('unrecognised move code')


def legal(tri, kind, dimension):
    guard = {'P': tri.hasPachner, 'Z': tri.has20, 'C': tri.hasCollapseEdge,
             'S': tri.hasShellBoundary, 'F': tri.has44}[kind]
    options = []
    for face in (tri.simplices() if dimension == 4 else tri.faces(dimension)):
        if kind == 'Z' and face.degree() != 2:
            continue
        if kind == 'P' and dimension < 4 and face.degree() != 5 - dimension:
            continue
        if kind == 'C' and face.vertex(0) == face.vertex(1):
            continue
        if kind == 'S' and not any(face.adjacentSimplex(i) is None for i in range(5)):
            continue
        if guard(face):
            options.append([kind, dimension, face.index()])
    return options


def restore_prefix(directory, original, endpoint):
    """Recover valid Regina caches by replay, never by changing a cached flag."""
    directory, original, endpoint = Path(directory), Path(original), Path(endpoint)
    metadata = json.loads((directory / 'SOURCE.json').read_text())
    checkpoint = json.loads((directory / 'CHECKPOINT.json').read_text())
    require(sha(original) == metadata['source_sha256'], 'prefix source checksum differs')
    tri = load4(original)
    admission(tri)
    current = digest(tri)
    require(current == metadata['start_tight_sha256'], 'wrong prefix start')
    n = 0
    with (directory / 'moves.jsonl').open() as stream:
        for line in stream:
            if n >= checkpoint['moves']:
                break
            row = json.loads(line)
            require(row['index'] == n and row['before_sha256'] == current,
                    'prefix chain discontinuity')
            perform(tri, row['operation'])
            current = digest(tri)
            require(current == row['after_sha256'] and tri.size() == row['pentachora'],
                    'prefix replay differs')
            n += 1
    require(n == checkpoint['moves'] and current == checkpoint['current_tight_sha256'],
            'prefix endpoint differs')
    require(tri == load4(endpoint), 'prefix does not reach requested source')
    return tri, {'directory': str(directory.resolve()), 'original_file': str(original.resolve()),
                 'original_sha256': sha(original), 'moves': n,
                 'endpoint_sha256': sha(endpoint), 'end_tight_sha256': current,
                 'all_move_guards_checked': True, 'cached_flags_manually_modified': False}


def replay(directory, source, best, output, prefix_directory=None, prefix_source=None):
    directory, source = Path(directory), Path(source)
    meta = json.loads((directory / 'SOURCE.json').read_text())
    end = json.loads((directory / 'CHECKPOINT.json').read_text())
    require(sha(source) == meta['source_sha256'], 'source file checksum differs')
    if prefix_directory is not None:
        require(prefix_source is not None, 'prefix source is required')
        tri, prefix_receipt = restore_prefix(prefix_directory, prefix_source, source)
    else:
        tri, prefix_receipt = load4(source), None
    checked_source = admission(tri)
    current = digest(tri)
    require(current == meta['start_tight_sha256'], 'wrong labelled source')
    limit = end['best_moves'] if best else end['moves']
    expected = end['best_tight_sha256'] if best else end['current_tight_sha256']
    counts, count = Counter(), 0
    started = time.monotonic()
    with (directory / 'moves.jsonl').open() as stream:
        for line in stream:
            if count >= limit:
                break
            row = json.loads(line)
            require(row['index'] == count and row['before_sha256'] == current,
                    'discontinuous move log')
            perform(tri, row['operation'])
            current = digest(tri)
            require(current == row['after_sha256'] and tri.size() == row['pentachora'],
                    'move replay differs')
            counts[str(row['operation'][:2])] += 1
            count += 1
            if count % 100 == 0:
                print(json.dumps({'replayed': count, 'pentachora': tri.size()}), flush=True)
    require(count == limit and current == expected, 'incomplete or wrong endpoint')
    saved = load4(directory / ('best.rga' if best else 'current.rga'))
    require(saved == tri, 'replay differs from saved labelled checkpoint')
    result = {'status': 'PASS_STANDARDIZATION' if is_standard(tri) else 'PASS_PARTIAL_PATH',
              'moves': count, 'counts': dict(counts), 'source_sha256': sha(source),
              'source_admission': checked_source, 'end_admission': admission(tri),
              'end_tight_sha256': current, 'standard_endpoint': is_standard(tri),
              'elapsed_seconds': time.monotonic() - started,
              'restored_by_checked_prefix': prefix_receipt,
              'native_4d_simplify_used': False,
              'trust': 'Regina 7.4.1 guarded moves and source manifold admission'}
    save_json(output, result)
    print(json.dumps(result), flush=True)
    return result


def search(source, out, seconds, seed, growth, neutral_limit, collapse,
           prefix_directory=None, prefix_source=None):
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina 7.4.1')
    source, out = Path(source), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    stop_requested = False
    def request_stop(signum, frame):
        nonlocal stop_requested
        stop_requested = True
    previous_handlers = {sig: signal.signal(sig, request_stop)
                         for sig in (signal.SIGINT, signal.SIGTERM)}
    if prefix_directory is not None:
        require(prefix_source is not None, 'prefix source is required')
        print(json.dumps({'stage': 'replaying_source_prefix'}), flush=True)
        tri, prefix_receipt = restore_prefix(prefix_directory, prefix_source, source)
    else:
        tri, prefix_receipt = load4(source), None
    print(json.dumps({'stage': 'source_admission', 'pentachora': tri.size()}), flush=True)
    source_info = admission(tri)
    rng = random.Random(seed)
    started = time.monotonic()
    current = digest(tri)
    best_hash, best_size, best_move, n = current, tri.size(), 0, 0
    counts, neutral, max_size = Counter(), 0, tri.size()
    last_update = started
    save_json(out / 'SOURCE.json', {
        'schema': 'motoo.guarded-finite-path.v1', 'source_file': str(source.resolve()),
        'source_sha256': sha(source), 'start_tight_sha256': current,
        'source_admission': source_info, 'regina_distribution': '7.4.1',
        'restored_by_checked_prefix': prefix_receipt,
        'seed': seed, 'policy': {'growth': growth, 'neutral_limit': neutral_limit,
                               'collapse': collapse},
        'native_4d_simplify_used': False})
    save4(out / 'best.rga', tri)

    def checkpoint(status):
        save4(out / 'current.rga', tri)
        result = {'status': status, 'moves': n, 'best_moves': best_move,
                  'current_pentachora': tri.size(), 'best_pentachora': best_size,
                  'current_tight_sha256': current, 'best_tight_sha256': best_hash,
                  'current_rga_sha256': sha(out / 'current.rga'),
                  'best_rga_sha256': sha(out / 'best.rga'),
                  'counts': dict(counts), 'max_pentachora': max_size,
                  'elapsed_seconds': time.monotonic() - started,
                  'standard_endpoint': is_standard(tri)}
        save_json(out / 'CHECKPOINT.json', result)
        return result

    checkpoint('RUNNING')
    status = 'WALL_BUDGET'
    with (out / 'moves.jsonl').open('w') as log:
        try:
            while time.monotonic() - started < seconds:
                if stop_requested:
                    status = 'STOP_REQUESTED_AFTER_LAST_COMPLETE_MOVE'
                    break
                if is_standard(tri):
                    status = 'STANDARD_SPHERE_REACHED'
                    break
                options = []
                priorities = ([('C', 1)] if collapse else [])
                priorities += [('S', 4), ('Z', 2), ('Z', 1), ('P', 1), ('P', 0), ('Z', 0)]
                for kind, dimension in priorities:
                    options = legal(tri, kind, dimension)
                    if options:
                        break
                reducing = bool(options)
                if not options and neutral < neutral_limit:
                    options = legal(tri, 'P', 2) + legal(tri, 'F', 1)
                    neutral += 1
                if not options and growth and tri.size() + 2 <= best_size + growth:
                    options = legal(tri, 'P', 3)
                    neutral = 0
                if not options:
                    status = 'LOCAL_PLATEAU'
                    break
                op, before = rng.choice(options), current
                perform(tri, op)
                current = digest(tri)
                log.write(json.dumps({'index': n, 'operation': op,
                                      'before_sha256': before, 'after_sha256': current,
                                      'pentachora': tri.size()}, separators=(',', ':')) + '\n')
                log.flush()
                n += 1
                counts[str(op[:2])] += 1
                max_size = max(max_size, tri.size())
                if reducing:
                    neutral = 0
                if tri.size() < best_size:
                    best_size, best_hash, best_move = tri.size(), current, n
                    save4(out / 'best.rga', tri)
                checkpoint('RUNNING')
                if time.monotonic() - last_update >= 10:
                    print(json.dumps({'moves': n, 'current': tri.size(), 'best': best_size,
                                      'seconds': time.monotonic() - started}), flush=True)
                    last_update = time.monotonic()
        except KeyboardInterrupt:
            # SIGINT normally sets a flag above.  An unexpected exception
            # must never commit a possibly mutated, unlogged triangulation.
            raise RuntimeError('Unexpected interruption; replay the previously committed checkpoint')
    for sig, handler in previous_handlers.items():
        signal.signal(sig, handler)
    result = checkpoint(status)
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    s = sub.add_parser('search')
    s.add_argument('source', type=Path)
    s.add_argument('--out', type=Path, required=True)
    s.add_argument('--seconds', type=float, default=600)
    s.add_argument('--seed', type=int, default=0)
    s.add_argument('--growth', type=int, default=0)
    s.add_argument('--neutral-limit', type=int, default=40)
    s.add_argument('--no-collapse', action='store_true')
    s.add_argument('--prefix-directory', type=Path)
    s.add_argument('--prefix-source', type=Path)
    r = sub.add_parser('replay')
    r.add_argument('directory', type=Path)
    r.add_argument('--source', type=Path, required=True)
    r.add_argument('--best', action='store_true')
    r.add_argument('--out', type=Path, required=True)
    r.add_argument('--prefix-directory', type=Path)
    r.add_argument('--prefix-source', type=Path)
    a = p.parse_args()
    if a.command == 'search':
        search(a.source, a.out, a.seconds, a.seed, a.growth, a.neutral_limit, not a.no_collapse,
               a.prefix_directory, a.prefix_source)
    else:
        replay(a.directory, a.source, a.best, a.out, a.prefix_directory, a.prefix_source)

#!/usr/bin/env python3
"""Closed guarded excursions with bounded random candidate sampling.

This changes candidate selection only. Every accepted operation uses the
unchanged guarded kernel and immutable journal from reduce_closed_bursts.
After random probes fail, a complete cyclic scan still detects any legal
candidate. Selection is deliberately not promised uniform among legal moves.
"""
from __future__ import annotations
import argparse
from collections import Counter
import importlib.metadata
import json
from pathlib import Path
import random
import signal
import time
import regina
import reduce_closed_bursts as kernel


class Sampler:
    def __init__(self, rng, probes):
        self.rng, self.probes = rng, probes
        self.stats = Counter()

    def one(self, tri, kind, dim):
        if kind == 'C' and tri.countVertices() < 2:
            self.stats['single_vertex_collapse_skips'] += 1
            return None
        count = tri.size() if dim == 4 else tri.countFaces(dim)
        if not count:
            return None
        guard = {'P': tri.hasPachner, 'Z': tri.has20,
                 'C': tri.hasCollapseEdge, 'F': tri.has44}[kind]

        def candidate(index):
            self.stats['face_probes'] += 1
            face = tri.simplex(index) if dim == 4 else tri.face(dim, index)
            if kind == 'P' and dim < 4 and face.degree() != 5 - dim:
                return None
            if kind == 'Z' and face.degree() != 2:
                return None
            if kind == 'C' and face.vertex(0) == face.vertex(1):
                return None
            if kind == 'F' and face.degree() != 4:
                return None
            self.stats['guard_calls'] += 1
            if guard(face):
                return [kind, dim, index]
            return None

        seen = set()
        for _ in range(min(self.probes, count)):
            index = self.rng.randrange(count)
            if index in seen:
                continue
            seen.add(index)
            answer = candidate(index)
            if answer is not None:
                self.stats['sample_hits'] += 1
                return answer
        self.stats['fallback_scans'] += 1
        start = self.rng.randrange(count)
        for offset in range(count):
            index = (start + offset) % count
            if index in seen:
                continue
            answer = candidate(index)
            if answer is not None:
                self.stats['fallback_hits'] += 1
                return answer
        self.stats['empty_candidate_sets'] += 1
        return None


def search(a):
    kernel.require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina7.4.1')
    kernel.require(a.seconds > 0 and a.growth >= 4 and a.burst24 > 0 and
                   a.mix >= 0 and a.birth_every > 0 and a.birth_vertices >= 2 and
                   a.chunk_size > 0 and a.probes > 0, 'invalid policy')
    source, out = a.source.resolve(), a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    pins = {str(path): kernel.sha_file(path) for path in
            (source, Path(__file__).resolve(), Path(kernel.__file__).resolve())}
    tri = kernel.load_tri(source)
    print(json.dumps({'stage': 'source_admission', 'pentachora': tri.size()}), flush=True)
    admission = kernel.admit(tri)
    kernel.require(all(kernel.sha_file(p) == h for p, h in pins.items()),
                   'source or code changed during admission')
    started = time.monotonic()
    rng = random.Random(a.seed)
    sampler = Sampler(rng, a.probes)
    two_zero_order = [('Z', 1), ('Z', 2)] if a.edge20_first else [('Z', 2), ('Z', 1)]
    policy = {'growth': a.growth, 'burst24': a.burst24, 'mix': a.mix,
              'birth_every': a.birth_every, 'birth_target_vertices': a.birth_vertices,
              'chunk_size': a.chunk_size, 'initial_burst': a.initial_burst,
              'birth_cycles_preserve_vertices_from_edge_collapse': True,
              'candidate_selection': 'bounded_random_probes_then_complete_cyclic_scan',
              'random_probes': a.probes,
              'edge20_first': a.edge20_first,
              'two_zero_descent_order': two_zero_order,
              'mix_selection': 'random_move_family_first_then_other_family_if_empty',
              'uniform_distribution_among_all_legal_moves_asserted': False}
    kernel.atomic_json(out / 'SOURCE.json', {
        'schema': 'motoo.closed-bursts.source.v1',
        'source_file': str(source), 'source_sha256': pins[str(source)],
        'start_tight_sha256': kernel.state_sha(tri), 'source_admission': admission,
        'regina': '7.4.1', 'seed': a.seed, 'policy': policy,
        'code_sha256': pins[str(Path(__file__).resolve())],
        'guard_and_journal_code_sha256': pins[str(Path(kernel.__file__).resolve())],
        'pinned_files': pins,
        'target': 'explicit double of Delta4, all five gluings identity',
        'native_4d_simplify_used': False})
    journal = kernel.ChunkJournal(out, tri, a.chunk_size, started)
    stopping = False

    def stop(signum, frame):
        nonlocal stopping
        stopping = True

    handlers = {s: signal.signal(s, stop) for s in (signal.SIGINT, signal.SIGTERM)}
    cycle, reason = 0, 'WALL_BUDGET'

    def live():
        return not stopping and time.monotonic() - started < a.seconds and not kernel.standard_target(tri)

    def choose(kind, dim, phase):
        operation = sampler.one(tri, kind, dim)
        if operation is None:
            return False
        journal.record(tri, operation, phase, cycle)
        return True

    def descent(allow_collapse):
        changed = False
        while live():
            order = ([('C', 1)] if allow_collapse else []) + two_zero_order + [
                ('P', 1), ('P', 0), ('Z', 0)]
            if not any(choose(kind, dim, 'descent') for kind, dim in order):
                break
            changed = True
        return changed

    try:
        if not a.initial_burst:
            descent(True)
        while live():
            cycle += 1
            made_move = False
            ceiling = journal.best.size() + a.growth
            birth_cycle = (cycle - 1) % a.birth_every == 0
            if birth_cycle:
                for _ in range(max(1, a.birth_vertices - tri.countVertices())):
                    if not live() or tri.size() + 4 > ceiling:
                        break
                    made_move = choose('P', 4, 'forced_vertex_birth') or made_move
            for _ in range(a.burst24):
                if not live() or tri.size() + 2 > ceiling:
                    break
                if not choose('P', 3, 'growth_burst'):
                    break
                made_move = True
            for _ in range(a.mix):
                if not live():
                    break
                families = [('P', 2), ('F', 1)]
                if rng.randrange(2):
                    families.reverse()
                found = any(choose(kind, dim, 'mix') for kind, dim in families)
                if not found:
                    break
                made_move = True
            made_move = descent(not birth_cycle) or made_move
            if not made_move:
                reason = 'NO_ALLOWED_MOVE_WITHIN_POLICY'
                break
        if kernel.standard_target(tri):
            reason = 'STANDARD_SPHERE_REACHED'
        elif stopping:
            reason = 'STOP_REQUESTED_AFTER_COMPLETE_MOVE'
        result = journal.commit(tri, reason, cycle)
    finally:
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        stable = {p: kernel.sha_file(p) == h for p, h in pins.items()}
        kernel.atomic_json(out / 'SELECTION.json', {
            'schema': 'motoo.closed-sampled.selection.v1',
            'counts': dict(sampler.stats), 'policy': policy,
            'pinned_files_unchanged': stable,
            'seconds_after_admission': time.monotonic() - started})
        kernel.require(all(stable.values()), 'source or code changed during search')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--seconds', type=float, default=600)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--growth', type=int, default=40)
    p.add_argument('--burst24', type=int, default=12)
    p.add_argument('--mix', type=int, default=80)
    p.add_argument('--birth-every', type=int, default=3)
    p.add_argument('--birth-vertices', type=int, default=5)
    p.add_argument('--chunk-size', type=int, default=50)
    p.add_argument('--probes', type=int, default=64)
    p.add_argument('--edge20-first', action='store_true',
                   help='try edge 2-0 before triangle 2-0 during descent')
    p.add_argument('--initial-burst', action='store_true')
    search(p.parse_args())


if __name__ == '__main__':
    main()

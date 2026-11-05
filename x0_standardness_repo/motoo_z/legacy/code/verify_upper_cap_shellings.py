#!/usr/bin/env python3
"""Replay every upper-cap boundary shelling back to the exact nine-cell reference."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import regina


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(tri):
    return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()


def verify(upper, intermediate, reference, out):
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina 7.4.1')
    upper, intermediate, reference, out = map(Path, (upper, intermediate, reference, out))
    files = {path: sha(path) for path in (upper, intermediate, reference)}
    require(out.resolve() not in {path.resolve() for path in files}, 'output would overwrite input')
    started = time.monotonic()
    tri = regina.Triangulation4(regina.open(str(upper)))
    middle = regina.Triangulation4(regina.open(str(intermediate)))
    target = regina.Triangulation4(regina.open(str(reference)))
    require((tri.size(), middle.size(), target.size()) == (483, 474, 9),
            'unexpected upper-cap sizes')
    require(tri.isValid() and not tri.isIdeal() and tri.isConnected() and
            tri.isOrientable() and tri.countBoundaryComponents() == 1 and
            tri.boundaryComponent(0).isReal(), 'upper cap is not an admitted finite PL manifold')
    source_digest = digest(tri)
    events = []
    middle_seen = False
    while tri.size() > target.size():
        index = tri.size() - 1
        before = digest(tri)
        simplex = tri.simplex(index)
        require(tri.hasShellBoundary(simplex), 'upper-cap shelling guard rejected')
        require(tri.shellBoundary(simplex), 'upper-cap checked shelling failed')
        events.append({'index': len(events), 'removed_simplex': index,
                       'before_sha256': before, 'after_sha256': digest(tri)})
        if tri.size() == middle.size():
            require(tri.tightEncoding() == middle.tightEncoding(), '474-cell intermediate differs')
            middle_seen = True
    require(middle_seen and tri.tightEncoding() == target.tightEncoding(),
            'upper shellings do not reach exact reference')
    require(all(sha(path) == value for path, value in files.items()), 'an input changed during replay')
    result = {'status': 'PASS_ALL_474_UPPER_CAP_SHELLINGS', 'moves': len(events),
              'source_sha256': files[upper], 'intermediate_sha256': files[intermediate],
              'reference_sha256': files[reference], 'source_tight_sha256': source_digest,
              'end_tight_sha256': digest(tri), 'intermediate_recovered_exactly': middle_seen,
              'reference_recovered_exactly': True, 'all_native_shelling_guards_checked': True,
              'source_finite_PL_admission_checked': True,
              'native_4d_simplify_used': False, 'input_hashes_stable': True,
              'events': events, 'seconds': time.monotonic() - started,
              'checker_sha256': sha(__file__)}
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = out.with_suffix(out.suffix + '.tmp')
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    temporary.replace(out)
    print(json.dumps({k: v for k, v in result.items() if k != 'events'}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upper', required=True)
    parser.add_argument('--intermediate', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    verify(args.upper, args.intermediate, args.reference, args.out)

#!/usr/bin/env python3
"""Fail-closed version of the supplied marked-link comparison diagnostic.

Requires SnapPy. It preserves the comparison criterion of the supplied
cert.py but returns nonzero on an unestablished comparison and records the
positive cusp maps. A negative result means NOT_CERTIFIED, not non-isotopy.
The PD-to-component-order convention is the same as the supplied checker;
the declared archival order must agree with the exterior's cusp order.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    data = json.loads(path.read_text())
    pd = data['pd_code']
    require(isinstance(pd, list) and all(len(row) == 4 for row in pd), 'invalid PD row lengths')
    decorations = data.get('decorations')
    if decorations is None:
        decorations = [c.get('decoration', c.get('framing')) for c in data['components']]
    decorations = [None if x == 'x' else x for x in decorations]
    require(all(x is None or type(x) is int for x in decorations), 'invalid decorations')
    return data, pd, decorations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('first', type=Path)
    parser.add_argument('second', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'NOT_CERTIFIED',
              'first_path': str(args.first), 'second_path': str(args.second),
              'negative_result_is_nonisotopy_certificate': False,
              'component_order_scope': 'Identity on the cusp order induced by the supplied PD inputs; archive component correspondence is a prerequisite.'}
    code = 1
    try:
        A, pd_a, decorations_a = load(args.first)
        B, pd_b, decorations_b = load(args.second)
        result.update(first_sha256=sha(args.first), second_sha256=sha(args.second),
                      first_decorations=decorations_a, second_decorations=decorations_b)
        if decorations_a != decorations_b:
            result['reason'] = 'Framings or dotted roles differ in the declared component order.'
        elif 'order' in A and 'order' in B and A['order'] != B['order']:
            result['reason'] = 'The declared archive component orders differ.'
        else:
            import snappy
            result['snappy_version'] = getattr(snappy, '__version__', 'unknown')
            MA = snappy.Link(pd_a).exterior()
            MB = snappy.Link(pd_b).exterior()
            require(MA.num_cusps() == MB.num_cusps() == len(decorations_a), 'component/cusp count mismatch')
            sa = MA.isometry_signature(of_link=True)
            sb = MB.isometry_signature(of_link=True)
            result.update(first_decorated_signature=sa, second_decorated_signature=sb,
                          signatures_equal=(sa == sb))
            good = []
            for iso in MA.is_isometric_to(MB, return_isometries=True):
                images = list(iso.cusp_images())
                if images != list(range(MA.num_cusps())):
                    continue
                matrices = [[[int(m[i, j]) for j in range(2)] for i in range(2)] for m in iso.cusp_maps()]
                if all(abs(m[0][0]) == 1 and m[0][1] == 0 and m[1][0] == 0 and m[1][1] == m[0][0] for m in matrices):
                    good.append({'cusp_images': images, 'cusp_maps': matrices})
            result['accepted_marked_isometries'] = good
            if sa == sb and good:
                result['status'] = 'CERTIFIED_SAME_MARKED_LINK'
                result['reason'] = 'Positive exterior isometry with identity cusp correspondence and each peripheral matrix equal to plus or minus I; framing/dotted-role lists agree.'
                code = 0
            else:
                result['reason'] = 'The required positive marked identification was not established.'
    except Exception as exc:
        result['status'] = 'CHECK_ERROR'
        result['reason'] = type(exc).__name__ + ': ' + str(exc)
        code = 2
    result['exit_code'] = code
    payload = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(payload)
    sys.stdout.write(payload)
    return code


if __name__ == '__main__':
    raise SystemExit(main())

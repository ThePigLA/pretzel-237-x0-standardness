#!/usr/bin/env python3
"""Exact-query encoding and fallback controls for the recognition bridge."""
import hashlib
import json
import os
from pathlib import Path
import sys
import regina


def require(ok, message):
    if not ok:
        raise ValueError(message)


runtime_path = Path(os.environ['MOTOO_RECOGNITION_RUNTIME'])
root = runtime_path.parent
cases = [('sphere', regina.Example3.threeSphere(), 'isSphere', True, True),
         ('ball', regina.Example3.ball(), 'isBall', True, True),
         ('rp3', regina.Example3.lens(2, 1), 'isSphere', False, False),
         ('poincare', regina.Example3.poincare(), 'isSphere', False, False),
         ('solid_torus', regina.Example3.handlebody(1), 'isBall', False, False)]
records = []
for name, tri, method, expected_result, expected_shortcut in cases:
    before = len(json.loads(runtime_path.read_text())['calls'])
    encoding = tri.tightEncoding()
    digest = hashlib.sha256(encoding.encode('ascii')).hexdigest()
    result = bool(getattr(tri, method)())
    runtime = json.loads(runtime_path.read_text())
    calls = runtime['calls'][before:]
    require(bool(calls), ('recognition symbol was not intercepted', name))
    receipt = json.loads((root / calls[0]['directory'] / 'CALL.json').read_text())
    require(receipt['input_tight_sha256'] == digest, ('serialized a different native object', name))
    require(receipt['exact_encoding_roundtrip'], name)
    require(result == expected_result, (name, result))
    require(receipt['accepted_positive'] == expected_shortcut, (name, receipt))
    require(tri.tightEncoding() == encoding, ('input mutated', name))
    reloaded = regina.Triangulation3(regina.open(str(root / calls[0]['directory'] / 'source.rga')))
    require(reloaded.tightEncoding() == encoding, ('saved source differs', name))
    records.append({'name': name, 'result': result,
                    'expected_result': expected_result,
                    'shortcut_accepted': receipt['accepted_positive'],
                    'native_input_tight_sha256': digest,
                    'exact_native_query_serialization_verified': True,
                    'calls': len(calls)})

# Exercise the actual internal 4D->3D PLT calls without changing any kernel.
before = len(json.loads(runtime_path.read_text())['calls'])
closed = regina.Triangulation4()
a, b = closed.newSimplices(2)
for facet in range(5):
    a.join(facet, b, regina.Perm5())
require(closed.isValid() and not closed.isIdeal() and closed.isClosed(),
        'explicit doubled-simplex control failed manifold admission')
calls = json.loads(runtime_path.read_text())['calls'][before:]
require(len(calls) == closed.countVertices() == 5, 'expected five internal4D link calls')
expected_links = [hashlib.sha256(v.buildLink().tightEncoding().encode('ascii')).hexdigest()
                  for v in closed.vertices()]
require([row['input_tight_sha256'] for row in calls] == expected_links,
        'internal4D link serialization differs')
require(all(row['accepted_positive'] for row in calls), 'a control link was not proved')
report = {'status': 'PASS_EXACT_SERIALIZATION_POSITIVE_AND_ORIGINAL_FALLBACK_CONTROLS',
          'cases': records, 'native_4d_internal_calls_verified': 5,
          'native_4d_vertex_link_encodings_match': True,
          'four_dimensional_guard_kernel_modified': False,
          'python_optimization_level': sys.flags.optimize,
          'all_tests_use_explicit_runtime_checks': True}
(root / 'CONTROLS.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report), flush=True)

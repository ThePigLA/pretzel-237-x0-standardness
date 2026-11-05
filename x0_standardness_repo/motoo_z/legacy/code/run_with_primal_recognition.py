#!/usr/bin/env python3
"""Run a Python task with an audited positive 3D recognition shortcut.

Linux/x86-64, pinned Regina7.4.1 wheel only. The original4D move kernel is
unchanged. Each accepted sphere/ball result has a frozen primal certificate
and independent replay. Unknown/errors fall back to original recognition.
"""
from __future__ import annotations
import argparse
import contextlib
import ctypes
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import runpy
import shutil
import sys
import tempfile
import threading
import time
import traceback

EXPECTED_ENGINE_SHA256 = '42e0c5a630d71f98c6e5ff3c70c520d1c3d136151aad4ac37fe807004fbaca1d'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def run(bridge, proofs, command):
    bridge, proofs = Path(bridge).resolve(), Path(proofs).resolve()
    require(sys.platform == 'linux', 'this exact ELF bridge requires Linux')
    require(bridge.is_file(), 'compile the bridge before launching')
    if os.environ.get('MOTOO_PRIMAL_REEXEC') != '1':
        environment = os.environ.copy()
        existing = environment.get('LD_PRELOAD', '')
        require(not existing, 'unexpected existing LD_PRELOAD; explicit review is required')
        environment['LD_PRELOAD'] = str(bridge)
        environment['MOTOO_PRIMAL_REEXEC'] = '1'
        optimization = ['-O'] if sys.flags.optimize == 1 else \
                       (['-OO'] if sys.flags.optimize >= 2 else [])
        os.execve(sys.executable, [sys.executable, *optimization, *sys.argv], environment)
    require(os.environ.get('LD_PRELOAD') == str(bridge), 'wrong preloaded helper')
    import regina
    from primal_triangle_group import solve
    from verify_primal_triangle_group import verify
    require(importlib.metadata.version('regina') == '7.4.1', 'requires Regina7.4.1')
    engines = list(Path(regina.__file__).parent.glob('engine*.so'))
    require(len(engines) == 1, 'ambiguous Regina engine shared object')
    engine = engines[0].resolve()
    require(sha(engine) == EXPECTED_ENGINE_SHA256, 'unreviewed engine binary')
    proofs.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    runtime = {'schema': 'motoo.primal-recognition-runtime.v1',
               'recognition_backend': 'positive primal-group shortcut with independent proof replay; original Regina fallback',
               'original_engine': str(engine), 'original_engine_sha256': sha(engine),
               'regina_distribution': '7.4.1', 'bridge_sha256': sha(bridge),
               'python_optimization_level': sys.flags.optimize,
               'bridge_source_sha256': sha(Path(__file__).with_name('primal_recognition_bridge.cpp')),
               'bridge_export_map_sha256': sha(Path(__file__).with_name('primal_recognition_bridge.map')),
               'launcher_sha256': sha(__file__),
               'solver_sha256': sha(Path(__file__).with_name('primal_triangle_group.py')),
               'independent_verifier_sha256': sha(Path(__file__).with_name('verify_primal_triangle_group.py')),
               'intercepted_operations': ['Triangulation3.isSphere', 'Triangulation3.isBall'],
               'four_dimensional_guard_kernel_modified': False,
               'cached_object_fields_written_by_bridge': False,
               'command': command, 'calls': [], 'status': 'CONFIGURING'}
    pinned_files = {str(path): sha(path) for path in (
        engine, bridge, Path(__file__),
        Path(__file__).with_name('primal_recognition_bridge.cpp'),
        Path(__file__).with_name('primal_recognition_bridge.map'),
        Path(__file__).with_name('primal_triangle_group.py'),
        Path(__file__).with_name('verify_primal_triangle_group.py'))}
    runtime['pinned_files'] = pinned_files
    runtime['proof_execution'] = 'in-process decoded copy in private temporary directory; independently replayed and atomically published before acceptance'
    private_proofs = tempfile.TemporaryDirectory(prefix='motoo-primal-recognition-')
    private_root = Path(private_proofs.name)
    runtime['private_proof_directory'] = str(private_root)
    atomic_json(proofs / 'RUNTIME.json', runtime)
    library = ctypes.CDLL(str(bridge))
    callback_type = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int)
    library.motoo_primal_configure.argtypes = [ctypes.c_char_p, callback_type]
    library.motoo_primal_configure.restype = ctypes.c_int
    library.motoo_primal_error.restype = ctypes.c_char_p
    library.motoo_primal_counter.argtypes = [ctypes.c_uint]
    library.motoo_primal_counter.restype = ctypes.c_ulonglong
    successful = {}
    call_number = 0
    class Ineligible(Exception):
        pass
    def check_pins():
        require(all(sha(path) == expected for path, expected in pinned_files.items()),
                'a pinned code or binary file changed during execution')
    def callback(pointer, size, kind):
        nonlocal call_number
        call_number += 1
        call_name = 'call_%05d_%s' % (call_number, 'sphere' if kind == 0 else 'ball')
        directory = private_root / call_name
        published = proofs / call_name
        directory.mkdir()
        then = time.monotonic()
        receipt = {'call': call_number, 'kind': 'sphere' if kind == 0 else 'ball',
                   'accepted_positive': False, 'native_fallback_required': True}
        verdict = 0
        try:
            check_pins()
            require(kind in (0, 1), 'unknown recognition operation')
            raw = ctypes.string_at(pointer, size)
            encoding = raw.decode('ascii')
            tight_hash = hashlib.sha256(raw).hexdigest()
            receipt['input_tight_sha256'] = tight_hash
            tri = regina.Triangulation3.tightDecoding(encoding)
            require(tri.tightEncoding() == encoding, 'public encoding round trip differs')
            receipt['exact_encoding_roundtrip'] = True
            receipt['tetrahedra'] = tri.size()
            source = directory / 'source.rga'
            require(regina.PacketOfTriangulation3(tri).save(str(source)), 'source export failed')
            receipt['source_sha256'] = sha(source)
            if not (tri.isValid() and not tri.isIdeal() and tri.isConnected()
                    and tri.isOrientable()):
                raise Ineligible('link hypotheses fail')
            if kind == 0:
                if not tri.isClosed():
                    raise Ineligible('sphere link is not closed')
            else:
                if not (tri.hasBoundaryFacets() and tri.countBoundaryComponents() == 1):
                    raise Ineligible('ball link does not have one real boundary')
                bc = tri.boundaryComponent(0)
                if not (bc.isOrientable() and bc.eulerChar() == 2):
                    raise Ineligible('ball boundary is not S2')
            receipt['finite_valid_connected_orientable_hypotheses'] = True
            cache_key = (kind, encoding)
            if cache_key in successful:
                prior = successful[cache_key]
                receipt['proof_reference'] = prior
                receipt['previous_exact_encoding_proof_reused'] = True
                verdict = 1
            else:
                transcript = io.StringIO()
                with contextlib.redirect_stdout(transcript):
                    solved = solve(source, directory / 'proof', 120)
                    replayed = verify(source, directory / 'proof', directory / 'INDEPENDENT_REPLAY.json')
                (directory / 'proof_output.txt').write_text(transcript.getvalue())
                receipt['proof_result_sha256'] = sha(directory / 'proof/RESULT.json')
                receipt['independent_replay_sha256'] = sha(directory / 'INDEPENDENT_REPLAY.json')
                require(solved['source_sha256'] == replayed['source_sha256'] ==
                        receipt['source_sha256'], 'proof source differs from exported query')
                expected = 'sphere' if kind == 0 else 'ball'
                if solved['all_original_generators_proved_identity'] and \
                        replayed['status'] == 'PASS_INDEPENDENT_TRIVIAL_GROUP_PROOF' and \
                        replayed['recognised_link_kind_by_poincare_and_3d_schoenflies'] == expected:
                    verdict = 1
            require(sha(source) == receipt['source_sha256'], 'query source file changed')
            reloaded = regina.Triangulation3(regina.open(str(source)))
            require(reloaded.tightEncoding() == encoding and
                    sha(source) == receipt['source_sha256'],
                    'proof source is not the exact serialized native query')
            receipt['proof_source_reloaded_equals_native_query'] = True
            check_pins()
            receipt['accepted_positive'] = verdict == 1
            receipt['native_fallback_required'] = verdict != 1
            receipt['status'] = 'PROVED_POSITIVE' if verdict == 1 else 'INCONCLUSIVE_USE_ORIGINAL'
        except Ineligible as error:
            receipt['status'] = 'HYPOTHESES_NOT_MET_USE_ORIGINAL'
            receipt['reason'] = str(error)
            verdict = 0
        except BaseException as error:
            receipt['status'] = 'SHORTCUT_ERROR_USE_ORIGINAL'
            receipt['error'] = type(error).__name__ + ': ' + str(error)
            receipt['traceback'] = traceback.format_exc()
            verdict = -1
        finally:
            receipt['seconds'] = time.monotonic() - then
            atomic_json(directory / 'CALL.json', receipt)
            manifest = {str(path.relative_to(directory)): sha(path)
                        for path in directory.rglob('*') if path.is_file()}
            staging = proofs / (call_name + '.publishing')
            shutil.copytree(directory, staging)
            require(all(sha(staging / path) == expected
                        for path, expected in manifest.items()),
                    'published proof copy differs from closed private proof')
            staging.replace(published)
            require(all(sha(published / path) == expected
                        for path, expected in manifest.items()),
                    'atomic published proof differs from private proof')
            runtime['calls'].append({'directory': call_name,
                                     'receipt_sha256': sha(published / 'CALL.json'),
                                     'published_files_sha256': manifest,
                                     'accepted_positive': receipt['accepted_positive'],
                                     'input_tight_sha256': receipt.get('input_tight_sha256')})
            runtime['seconds'] = time.monotonic() - started
            atomic_json(proofs / 'RUNTIME.json', runtime)
            print(json.dumps({'recognition_call': call_number, 'kind': receipt['kind'],
                              'tetrahedra': receipt.get('tetrahedra'),
                              'accepted_positive': receipt['accepted_positive'],
                              'seconds': receipt['seconds']}), flush=True)
            if verdict == 1:
                successful.setdefault(cache_key, call_name)
        return verdict
    callback_lock = threading.RLock()
    runtime['callback_serialized_by_python_lock'] = True
    def safe_callback(pointer, size, kind):
        # This outer boundary also covers receipt writes and printing in the
        # inner finally block. CFUNCTYPE must never see an escaping exception.
        try:
            with callback_lock:
                return callback(pointer, size, kind)
        except BaseException:
            return -1
    callback_reference = callback_type(safe_callback)
    require(library.motoo_primal_configure(str(engine).encode(), callback_reference) == 1,
            'bridge configuration failed: ' + str(library.motoo_primal_error()))
    runtime['status'] = 'RUNNING'
    atomic_json(proofs / 'RUNTIME.json', runtime)
    os.environ['MOTOO_RECOGNITION_RUNTIME'] = str(proofs / 'RUNTIME.json')
    require(command, 'missing target Python script')
    target = Path(command[0]).resolve()
    require(target.is_file(), 'target script is missing')
    sys.path.insert(0, str(target.parent))
    sys.argv = [str(target), *command[1:]]
    try:
        runpy.run_path(str(target), run_name='__main__')
        runtime['status'] = 'TARGET_COMPLETED'
    except BaseException as error:
        runtime['status'] = 'TARGET_EXITED_OR_RAISED'
        runtime['target_exception'] = type(error).__name__ + ': ' + str(error)
        raise
    finally:
        runtime['bridge_counters'] = dict(zip(
            ('sphere_calls', 'ball_calls', 'proved_positive', 'original_fallbacks', 'bridge_errors'),
            [library.motoo_primal_counter(i) for i in range(5)]))
        runtime['seconds'] = time.monotonic() - started
        runtime['pinned_files_unchanged_at_exit'] = {
            path: sha(path) == expected for path, expected in pinned_files.items()}
        if not all(runtime['pinned_files_unchanged_at_exit'].values()):
            runtime['status'] = 'CODE_OR_BINARY_BINDING_CHANGED'
        atomic_json(proofs / 'RUNTIME.json', runtime)
        require(all(runtime['pinned_files_unchanged_at_exit'].values()),
                'pinned code or binary changed during run')
        # Retain private originals throughout the run and until every final
        # published file has been rebound against its recorded checksum.
        require(all(sha(proofs / call['directory'] / relative) == expected
                    for call in runtime['calls']
                    for relative, expected in call['published_files_sha256'].items()),
                'a published proof file changed before run completion')
        private_proofs.cleanup()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bridge', type=Path,
                        default=Path(__file__).with_name('primal_recognition_bridge.so'))
    parser.add_argument('--proofs', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    run(args.bridge, args.proofs, command)

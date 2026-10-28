#!/usr/bin/env python3
"""Replay the fixed X0 reduction and its connected paths to the standard sphere.

Full replay uses the exact construction output, the elementary certificates,
the recorded native macro chain, and two independent terminal audits. It does
not rerun a discovery search or prove the geometric construction identification.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from common import (ROOT, CERTIFICATE_SHA256, ENGINE_SHA256, FACETS_SHA256,
                    G6, NATIVE_SOURCE_SHA256, S2, candidate_identity,
                    check_release_manifest, code_integrity, file_record,
                    input_integrity, read_json, require, require_new_output,
                    sha256)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def write_report(path, report):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def child_environment():
    env = dict(os.environ)
    env.pop('PYTHONOPTIMIZE', None)
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    return env


def native_environment(interpreter, full):
    probe = r'''
import hashlib,importlib.metadata,json,platform,sys
from pathlib import Path
import regina
p=Path(regina.engine.__file__)
print(json.dumps({'python':platform.python_version(),'python_major_minor':list(sys.version_info[:2]),
 'platform':platform.platform(),'system':platform.system(),'machine':platform.machine(),
 'regina_distribution_version':importlib.metadata.version('regina'),
 'regina_api_version':regina.versionString(),'engine_filename':p.name,
 'engine_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'optimization':sys.flags.optimize}))
'''
    result = subprocess.run([interpreter, '-E', '-B', '-c', probe],
                            env=child_environment(), capture_output=True, text=True, timeout=30)
    require(result.returncode == 0,
            'native interpreter could not load Regina; install environment/requirements.txt.\n' +
            result.stderr[-3000:])
    env = json.loads(result.stdout)
    require(env['regina_distribution_version'] == '7.4.1' and env['regina_api_version'] == '7.4',
            'this release requires Regina distribution 7.4.1 (API 7.4)')
    require(env['optimization'] == 0, 'native interpreter has assertions disabled')
    if full:
        require(env['system'] == 'Linux' and env['machine'] == 'x86_64'
                and env['python_major_minor'] == [3, 12],
                'full macro replay requires Linux x86_64 and CPython 3.12; '
                'see environment/README.md for the container command')
        require(env['engine_sha256'] == ENGINE_SHA256,
                'native engine differs from the recorded Linux CPython 3.12 wheel')
    return env


def run_stage(label, command, output, report, report_path, wall_limit):
    log = output / 'logs' / (label + '.log')
    require(not log.exists(), 'stage log already exists: ' + str(log))
    start = time.monotonic()
    print('REPLAY ' + label, flush=True)
    entry = {'status': 'RUNNING', 'started_utc': utc_now(),
             'command': list(map(str, command)), 'log': str(log.relative_to(output))}
    report['stages'][label] = entry
    write_report(report_path, report)
    child = None
    try:
        with log.open('x') as stream:
            child = subprocess.Popen(list(map(str, command)), stdout=stream,
                                     stderr=subprocess.STDOUT, start_new_session=True,
                                     env=child_environment(), cwd=ROOT)
            while True:
                try:
                    code = child.wait(timeout=25)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - start
                    print(f'{label}: {elapsed:.0f}s elapsed; log {log.name}', flush=True)
                    require(elapsed <= wall_limit, label + ' exceeded the stage wall limit')
        require(code == 0, f'{label} failed with exit {code}; see {log}')
        entry.update({'status': 'PASS', 'returncode': code})
    except BaseException:
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                child.wait()
        entry['status'] = 'FAIL_OR_INTERRUPTED'
        raise
    finally:
        entry['wall_seconds'] = time.monotonic() - start
        entry['finished_utc'] = utc_now()
        if log.is_file():
            entry['log_identity'] = file_record(log)
        write_report(report_path, report)
    print(f'{label}: PASS ({entry["wall_seconds"]:.1f}s)', flush=True)


def bind_elementary(output, candidate):
    path = output / 'elementary/audit_steps012.json'
    e = read_json(path)
    require(e['status'] == 'PASS' and e['source_files_unchanged'] is True
            and e['replay_equals_macro_source_exactly_including_ordered_facets_and_gluings'] is True,
            'elementary replay is incomplete')
    require(e['master_facets_sha256'] == FACETS_SHA256 and e['prefix_contractions'] == 3308
            and e['step2_moves'] == 629874 and e['link_moves'] == 529460
            and e['endpoint_facets'] == 1332 and e['endpoint_vertices'] == 45
            and e['endpoint_f_vector'] == [45, 795, 2750, 3330, 1332],
            'elementary counts or master interface differ')
    require(e['input_file_sha256']['inputs/marked_control.json.gz'] ==
            candidate['files']['marked_candidate.json.gz']['sha256'], 'selected master was not replayed')
    require(e['input_file_sha256']['inputs/cover_regina.json.gz'] ==
            candidate['files']['cover_regina.json.gz']['sha256'], 'selected cover was not checked')
    require(e['audit_script_sha256'] == sha256(ROOT / 'source/elementary/audit_steps012.py'),
            'elementary checker identity differs')
    for stage, script in [('step0', 'step0_links/check_links.py'),
                          ('step1', 'step1_contractions/check_contractions.py'),
                          ('step2', 'step2_flips/replay_moves.py')]:
        require(e[stage]['returncode'] == 0 and e[stage]['assertions_enabled'] is True
                and e[stage]['script_sha256'] == sha256(ROOT / 'elementary' / script),
                'elementary subchecker did not complete: ' + stage)
    replay_path = output / 'elementary/replayed_1332.json.gz'
    require(sha256(replay_path) == e['replayed_export_sha256'], 'fresh elementary endpoint changed')
    actual = read_json(replay_path)
    source_path = ROOT / 'elementary/step3_macro/x0_C1332_regina.json.gz'
    native_source = read_json(source_path)
    require(all(actual[key] == native_source[key] for key in ('simplex_count', 'facets', 'gluings')),
            'fresh ordered elementary endpoint does not equal the macro source')
    require(e['macro_source_sha256'] == sha256(source_path) == NATIVE_SOURCE_SHA256,
            'wrong macro source bytes')
    return e, {'status': 'PASS', 'receipt': file_record(path),
               'fresh_endpoint': file_record(replay_path),
               'selected_candidate_replayed': True,
               'ordered_facets_and_gluings_equal_native_source': True,
               'macro_source_sha256': NATIVE_SOURCE_SHA256}


def bind_native(output, elementary):
    path = output / 'native/independent_replay.json'
    n = read_json(path)
    require(n['status'] == 'PASS' and n['input_files_unchanged'] is True,
            'native replay is incomplete')
    require(n['source']['file_sha256'] == elementary['macro_source_sha256'] == NATIVE_SOURCE_SHA256,
            'elementary/native interface mismatch')
    require(n['verifier_file_sha256'] == sha256(ROOT / 'source/native/replay_native_and_cap.py'),
            'native checker identity differs')
    require(n['certificate_file_sha256'] == sha256(ROOT / 'elementary/step3_macro/x0_macro_cert.json'),
            'wrong native certificate')
    require(n['macros']['blocks_replayed'] == n['macros']['numbered_hashes_checked'] == 5173
            and n['macros']['guarded_2_4_moves'] == 54563
            and n['macros']['final_size'] == 5, 'native macro count differs')
    require(n['capping']['status'] == 'PASS' and n['capping']['target_isoSig'] == G6
            and n['capping']['no_simplification_performed_after_capping'] is True,
            'direct capped endpoint was not verified')
    for name, record in n['outputs'].items():
        require(Path(name).name == name and name not in ('.', '..'), 'unsafe native output name')
        require(file_record(output / 'native' / name) == record, 'fresh native output changed: ' + name)
    require(set(n['outputs']) == {'replayed_T5.json', 'capped_G6.json',
                                 'macro_block_sizes.json', 'capping_receipt.json'},
            'native output set differs')
    return n, {'status': 'PASS', 'receipt': file_record(path),
               'capped_endpoint': file_record(output / 'native/capped_G6.json'),
               'source_matches_fresh_elementary_endpoint': True,
               'direct_cap_target': G6, 'post_cap_simplification_used': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path,
                        help='construction output candidate/ directory; default: supplied exact reference')
    parser.add_argument('--output', type=Path, help='new directory for this run, logs, and fresh receipts')
    parser.add_argument('--scope', choices=('full', 'elementary', 'terminal'), default='full',
                        help='full chain, elementary part only, or small terminal paths only')
    parser.add_argument('--native-python', default=sys.executable,
                        help='Python interpreter containing the pinned Regina package')
    parser.add_argument('--integrity-only', action='store_true',
                        help='check all released bytes; performs no mathematical replay')
    parser.add_argument('--keep-work', action='store_true', help='retain the isolated input staging tree')
    parser.add_argument('--stage-timeout', type=int, default=1800,
                        help='outer per-stage wall limit in seconds; imported checkers retain their own limits')
    args = parser.parse_args()
    require(sys.flags.optimize == 0, 'run without -O/-OO; elementary checkers require assertions')
    require(os.name == 'posix', 'this driver uses POSIX subprocess and resource controls')
    require(args.stage_timeout > 0, '--stage-timeout must be positive')
    require(not (args.scope == 'terminal' and args.candidate is not None),
            '--scope terminal does not consume a construction candidate')
    inputs_before = input_integrity()
    code_before = code_integrity()
    if args.integrity_only:
        release = check_release_manifest()
        print(json.dumps({'status': 'INTEGRITY_PASS', 'mathematical_replay_performed': False,
                          'imported_file_count': inputs_before['imported_file_count'], **release}, indent=2))
        return 0
    require(args.output is not None, '--output is required for a replay')
    candidate_path = None
    candidate = None
    if args.scope != 'terminal':
        candidate_path = (args.candidate or ROOT / 'reference/original/candidate').absolute()
        candidate = candidate_identity(candidate_path)
        candidate_path = Path(candidate['directory'])
        candidate['selection'] = 'explicit_construction_candidate' if args.candidate else 'packaged_reference'
    output = require_new_output(args.output, candidate_path)
    native_python = args.native_python
    if os.path.dirname(native_python):
        native_python = os.path.abspath(os.path.expanduser(native_python))
    else:
        native_python = shutil.which(native_python)
    require(native_python is not None, 'native Python interpreter not found')
    environment = None
    if args.scope != 'elementary':
        environment = native_environment(native_python, args.scope == 'full')
    output.mkdir(parents=True)
    (output / 'logs').mkdir()
    start = time.monotonic()
    report_path = output / 'REDUCTION_REPLAY.json'
    report = {'schema': 'x0-reduction-fresh-composed-replay-v1', 'status': 'RUNNING',
              'scope': args.scope, 'started_utc': utc_now(), 'candidate': candidate,
              'code': code_before, 'inputs': inputs_before, 'environment': environment,
              'stages': {}, 'interfaces': {},
              'construction_geometric_identification_replayed': False,
              'historical_success_receipts_used_for_acceptance': False,
              'discovery_search_performed': False, 'whole_chain_is_pure_Pachner': False}
    write_report(report_path, report)
    py = [sys.executable, '-E', '-B', '-u']
    native_py = [native_python, '-E', '-B', '-u']
    work = output / 'work'
    try:
        elementary = None
        if args.scope != 'terminal':
            shutil.copytree(ROOT / 'elementary', work / 'elementary')
            shutil.copytree(ROOT / 'reference/original', work / 'original')
            for src_name, review_name in [('marked_candidate.json.gz', 'marked_control.json.gz'),
                                          ('cover_regina.json.gz', 'cover_regina.json.gz')]:
                shutil.copy2(candidate_path / src_name, work / 'elementary/inputs' / review_name)
                shutil.copy2(candidate_path / src_name, work / 'original/candidate' / src_name)
            run_stage('elementary', py + [ROOT / 'source/elementary/audit_steps012.py',
                      '--reviewer-root', work / 'elementary', '--original-root', work / 'original',
                      '--output-dir', output / 'elementary'], output, report, report_path, args.stage_timeout)
            elementary, binding = bind_elementary(output, candidate)
            report['interfaces']['construction_to_elementary_to_native'] = binding
            if args.scope == 'full':
                run_stage('native_macros_and_cap', native_py + [ROOT / 'source/native/replay_native_and_cap.py',
                          '--source', ROOT / 'elementary/step3_macro/x0_C1332_regina.json.gz',
                          '--certificate', ROOT / 'elementary/step3_macro/x0_macro_cert.json',
                          '--output', output / 'native'], output, report, report_path, args.stage_timeout)
                _, binding = bind_native(output, elementary)
                report['interfaces']['native_to_closed_six'] = binding
        if args.scope != 'elementary':
            terminal_dir = output / 'terminal'
            terminal_dir.mkdir()
            cert = ROOT / 'certificates/standardization.json'
            run_stage('terminal_supplied_replay', native_py + [ROOT / 'source/terminal/verify_x0_standard.py',
                      cert], output, report, report_path, args.stage_timeout)
            run_stage('terminal_independent_replay', native_py + [ROOT / 'source/terminal/independent_replay.py',
                      '--certificate', cert, '--output', terminal_dir / 'independent_replay.json'],
                      output, report, report_path, args.stage_timeout)
            run_stage('terminal_independent_two_zero', native_py + [ROOT / 'source/terminal/independent_two_zero.py',
                      cert, '--output', terminal_dir / 'independent_two_zero.json'],
                      output, report, report_path, args.stage_timeout)
            command = native_py + [ROOT / 'terminal.py', '--certificate', cert,
                                   '--replay', terminal_dir / 'independent_replay.json',
                                   '--two-zero', terminal_dir / 'independent_two_zero.json',
                                   '--output', terminal_dir / 'bindings.json']
            if args.scope == 'full':
                command += ['--capped', output / 'native/capped_G6.json']
            run_stage('terminal_composition', command, output, report, report_path, args.stage_timeout)
            terminal_binding = read_json(terminal_dir / 'bindings.json')
            require(terminal_binding['status'] == 'PASS'
                    and terminal_binding['schema'] == 'x0-reduction-terminal-binding-v1'
                    and terminal_binding['checker_sha256'] == sha256(ROOT / 'terminal.py')
                    and terminal_binding['certificate']['sha256'] == CERTIFICATE_SHA256
                    and terminal_binding['standard_anchor']['isoSig'] == S2
                    and terminal_binding['exact_path_intersections_checked'] is True
                    and terminal_binding['all_moves'] == 8436
                    and terminal_binding['two_zero_moves'] == 594
                    and terminal_binding['capped_endpoint_bound'] == (args.scope == 'full'),
                    'terminal composition is incomplete or has the wrong scope')
            if args.scope == 'full':
                require(terminal_binding['capped_endpoint']['sha256'] ==
                        sha256(output / 'native/capped_G6.json'), 'terminal cap binding differs')
            for label, filename in [('independent_replay', 'independent_replay.json'),
                                    ('independent_two_zero', 'independent_two_zero.json')]:
                require(terminal_binding['receipts'][label]['sha256'] == sha256(terminal_dir / filename),
                        'fresh terminal receipt changed: ' + label)
            standard_export = {
                'simplex_count': 2,
                'gluings': [[0, facet, 1, [0, 1, 2, 3, 4]] for facet in range(5)],
                'isoSig': S2,
                'f_vector': terminal_binding['standard_anchor']['properties']['f_vector'],
                'construction': 'Two distinct four-simplices joined along corresponding facets by identity',
            }
            (terminal_dir / 'standard_sphere.json').write_text(json.dumps(standard_export, indent=2) + '\n')
            report['interfaces']['terminal_paths_to_explicit_standard_double'] = {
                'status': 'PASS', 'receipt': file_record(terminal_dir / 'bindings.json'),
                'certificate_sha256': CERTIFICATE_SHA256, 'standard_sphere_isoSig': S2,
                'fresh_capped_endpoint_supplied': args.scope == 'full'}
        require(input_integrity() == inputs_before, 'module proof inputs changed during replay')
        require(code_integrity() == code_before, 'replay code changed during run')
        if candidate is not None:
            after = candidate_identity(candidate_path)
            require(after['files'] == candidate['files'], 'selected construction inputs changed during run')
        report['all_code_and_inputs_unchanged'] = True
        if work.exists() and not args.keep_work:
            shutil.rmtree(work)
        report['worktree_retained'] = work.exists()
        if args.scope == 'full':
            report['status'] = 'PASS_MASTER_PL_HOMEOMORPHIC_TO_STANDARD_S4'
            report['claim'] = ('The exact 103744-pentachoron master is PL homeomorphic to the '
                               'explicit double of a four-simplex, with the documented native move dependencies. '
                               'Its identification with the geometric X0 construction is supplied by construction/.')
        elif args.scope == 'elementary':
            report['status'] = 'PASS_ELEMENTARY_ONLY'
            report['claim'] = 'The exact master is a closed PL manifold and reduces to the exact 1332-facet macro source.'
        else:
            report['status'] = 'PASS_TERMINAL_ONLY'
            report['claim'] = 'The specified closed six-pentachoron triangulation is PL homeomorphic to the explicit standard double.'
        report['master_to_standard_sphere_replayed'] = args.scope == 'full'
        report['result_isoSig'] = S2 if args.scope != 'elementary' else None
        report['fresh_output_files'] = {
            str(p.relative_to(output)): file_record(p)
            for folder in ('elementary', 'native', 'terminal', 'logs')
            for p in sorted((output / folder).rglob('*')) if p.is_file()}
        return_code = 0
    except BaseException as error:
        report['status'] = 'INTERRUPTED' if isinstance(error, KeyboardInterrupt) else 'FAIL'
        report['error_type'] = type(error).__name__
        report['error'] = str(error)
        report['master_to_standard_sphere_replayed'] = False
        report['worktree_retained'] = work.exists()
        print(report['status'] + ': ' + str(error), file=sys.stderr, flush=True)
        return_code = 130 if isinstance(error, KeyboardInterrupt) else 1
    report['finished_utc'] = utc_now()
    report['wall_seconds'] = time.monotonic() - start
    write_report(report_path, report)
    print(report['status'], flush=True)
    print('Receipt: ' + str(report_path), flush=True)
    return return_code


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    try:
        raise SystemExit(main())
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print('PREFLIGHT_FAIL: ' + str(error), file=sys.stderr)
        raise SystemExit(1)

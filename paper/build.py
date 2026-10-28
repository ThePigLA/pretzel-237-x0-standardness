#!/usr/bin/env python3
"""Compile the editable manuscript in a fresh directory without changing its sources.

The LaTeX source is self-contained and embeds its bibliography. This command
runs pdfLaTeX twice with shell escape disabled; it performs no topology replay.
When Poppler's tools are available, it records PDF page counts and compares
normalized extracted text with the supplied reading copy. Visual review is
a separate step.
"""
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import unicodedata


ROOT = Path(__file__).resolve().parent
STEM = 'x0_standardness_proof_note'
SOURCE = ROOT / (STEM + '.tex')
REFERENCE = ROOT / (STEM + '.pdf')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    raw = Path(path).read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def absolute_without_symlinks(path):
    path = Path(os.path.abspath(Path(path).expanduser()))
    for ancestor in (path, *path.parents):
        require(not ancestor.is_symlink(), 'symlinked output path is not accepted: ' + str(ancestor))
    return path


def source_inventory():
    require(SOURCE.is_file(), 'missing LaTeX source: ' + str(SOURCE))
    paths = [SOURCE, REFERENCE, ROOT / (STEM + '.txt'), ROOT / 'bibliography.bib', Path(__file__)]
    return {p.name: identity(p) for p in paths if p.is_file()}


def executable(name):
    expanded = os.path.expanduser(name)
    if os.path.dirname(expanded):
        path = Path(os.path.abspath(expanded))
        return str(path) if path.is_file() and os.access(path, os.X_OK) else None
    return shutil.which(expanded)


def write_report(path, report):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def command_output(command, timeout=30):
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    require(result.returncode == 0,
            'command failed: ' + str(command[0]) + '\n' + (result.stdout + result.stderr)[-4000:])
    return result.stdout


def normalized_text(text):
    # Ignore whitespace/layout differences only. Mathematical symbols, words,
    # identifiers, command flags, and printed numbers are retained.
    return ' '.join(unicodedata.normalize('NFC', text).split())


def pdf_information(path, tool):
    raw = command_output([tool, str(path)])
    match = re.search(r'^Pages:\s+(\d+)\s*$', raw, re.MULTILINE)
    require(match is not None, 'pdfinfo did not report a page count')
    return {'pages': int(match.group(1)), 'pdfinfo': raw}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path,
                        help='new output directory outside paper/')
    parser.add_argument('--engine', default='pdflatex',
                        help='pdfLaTeX-compatible executable or absolute path (default: pdflatex)')
    parser.add_argument('--timeout', type=int, default=180,
                        help='wall-time limit in seconds for each compilation pass')
    parser.add_argument('--verify-reference', '--require-reference-match', dest='verify_reference', action='store_true',
                        help='require matching page counts and normalized PDF text; needs pdfinfo and pdftotext')
    parser.add_argument('--keep-work', action='store_true', help='retain TeX auxiliary files and the copied source')
    args = parser.parse_args()
    require(args.timeout > 0, '--timeout must be positive')
    output = absolute_without_symlinks(args.output)
    require(not output.exists(), '--output must name a new directory: ' + str(output))
    require(output != ROOT and ROOT not in output.parents and output not in ROOT.parents,
            '--output must be outside paper/ and must not contain the manuscript directory')
    engine = executable(args.engine)
    require(engine is not None, 'pdfLaTeX executable not found: ' + args.engine)
    pdfinfo = shutil.which('pdfinfo')
    pdftotext = shutil.which('pdftotext')
    if args.verify_reference:
        require(REFERENCE.is_file(), 'reference PDF is required for comparison')
        require(pdfinfo is not None and pdftotext is not None,
                '--verify-reference needs Poppler pdfinfo and pdftotext')
    before = source_inventory()
    version_result = subprocess.run([engine, '--version'], capture_output=True, text=True, timeout=30)
    require(version_result.returncode == 0, 'cannot query the selected TeX engine')
    version = version_result.stdout.splitlines()[0] if version_result.stdout else ''
    require('pdfTeX' in version, 'the source requires pdfLaTeX/pdfTeX (it uses pdfTeX Unicode primitives)')
    output.mkdir(parents=True)
    work = output / 'work'
    logs = output / 'logs'
    work.mkdir()
    logs.mkdir()
    shutil.copy2(SOURCE, work / SOURCE.name)
    started = time.monotonic()
    report_path = output / 'BUILD_REPORT.json'
    report = {'schema': 'x0-paper-fresh-build-v1', 'status': 'RUNNING',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'engine': {'executable': engine, 'version': version},
              'source_files': before, 'passes': [], 'shell_escape_enabled': False,
              'topology_replay_performed': False, 'visual_review_performed_by_builder': False,
              'reference_match_required': args.verify_reference}
    write_report(report_path, report)
    returncode = 0
    try:
        for number in (1, 2):
            command = [engine, '-no-shell-escape', '-halt-on-error',
                       '-interaction=nonstopmode', '-file-line-error', SOURCE.name]
            log_path = logs / ('pdflatex_pass_' + str(number) + '.log')
            pass_start = time.monotonic()
            print('COMPILE pass ' + str(number), flush=True)
            with log_path.open('x') as log:
                result = subprocess.run(command, cwd=work, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=args.timeout)
            entry = {'pass': number, 'command': command,
                     'returncode': result.returncode, 'wall_seconds': time.monotonic() - pass_start,
                     'log': str(log_path.relative_to(output)), 'log_identity': identity(log_path)}
            report['passes'].append(entry)
            write_report(report_path, report)
            require(result.returncode == 0, 'compilation failed; see ' + str(log_path))
        latex_log = work / (STEM + '.log')
        require(latex_log.is_file(), 'TeX did not produce a log')
        log_text = latex_log.read_text(errors='replace')
        warning_patterns = (
            r'Overfull \\[hv]box', r'Underfull \\[hv]box', r'Missing character:',
            r'Undefined control sequence', r'LaTeX Warning:', r'LaTeX Font Warning:',
            r'Package \S+ Warning:', r'Class \S+ Warning:',
        )
        warnings = [line for line in log_text.splitlines()
                    if any(re.search(pattern, line) for pattern in warning_patterns)]
        report['typesetting_warnings'] = warnings
        shutil.copy2(latex_log, logs / 'final_latex.log')
        require(not warnings, 'typesetting warnings require review; see logs/final_latex.log')
        built_pdf = work / (STEM + '.pdf')
        require(built_pdf.is_file() and built_pdf.stat().st_size > 0,
                'TeX did not produce a nonempty PDF')
        final_pdf = output / built_pdf.name
        shutil.copy2(built_pdf, final_pdf)
        report['output_pdf'] = {'path': final_pdf.name, **identity(final_pdf)}
        report['source_copy_unchanged'] = identity(work / SOURCE.name) == before[SOURCE.name]
        require(report['source_copy_unchanged'], 'copied TeX source changed during compilation')
        comparison = {'reference_present': REFERENCE.is_file(),
                      'byte_equal': identity(REFERENCE) == identity(final_pdf) if REFERENCE.is_file() else None,
                      'pdfinfo_available': pdfinfo is not None,
                      'pdftotext_available': pdftotext is not None,
                      'page_count_equal': None, 'normalized_text_equal': None,
                      'normalization': 'Unicode NFC and whitespace folding only'}
        if pdfinfo is not None:
            built_info = pdf_information(final_pdf, pdfinfo)
            report['output_pages'] = built_info['pages']
            (logs / 'built_pdfinfo.txt').write_text(built_info['pdfinfo'])
            if REFERENCE.is_file():
                ref_info = pdf_information(REFERENCE, pdfinfo)
                (logs / 'reference_pdfinfo.txt').write_text(ref_info['pdfinfo'])
                comparison['reference_pages'] = ref_info['pages']
                comparison['page_count_equal'] = built_info['pages'] == ref_info['pages']
        if pdftotext is not None:
            built_text = output / (STEM + '.txt')
            command_output([pdftotext, '-layout', str(final_pdf), str(built_text)])
            report['output_text'] = {'path': built_text.name, **identity(built_text)}
            if REFERENCE.is_file():
                ref_text = logs / 'reference_extracted_text.txt'
                command_output([pdftotext, '-layout', str(REFERENCE), str(ref_text)])
                comparison['normalized_text_equal'] = (
                    normalized_text(built_text.read_text()) == normalized_text(ref_text.read_text()))
                comparison['reference_extracted_text'] = identity(ref_text)
        report['reference_comparison'] = comparison
        if args.verify_reference:
            require(comparison['page_count_equal'] is True and comparison['normalized_text_equal'] is True,
                    'rebuilt PDF differs from the reference in page count or normalized extracted text')
        require(source_inventory() == before, 'manuscript source/reference files changed during build')
        report['source_files_unchanged'] = True
        if not args.keep_work:
            shutil.rmtree(work)
        report['worktree_retained'] = work.exists()
        report['status'] = 'PASS'
    except BaseException as error:
        report['status'] = 'INTERRUPTED' if isinstance(error, KeyboardInterrupt) else 'FAIL'
        report['error_type'] = type(error).__name__
        report['error'] = str(error)
        report['source_files_unchanged'] = source_inventory() == before
        report['worktree_retained'] = work.exists()
        returncode = 130 if isinstance(error, KeyboardInterrupt) else 1
        print(report['status'] + ': ' + str(error), file=sys.stderr, flush=True)
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    report['wall_seconds'] = time.monotonic() - started
    write_report(report_path, report)
    print(report['status'] + ': ' + str(report_path), flush=True)
    return returncode


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print('PREFLIGHT_FAIL: ' + str(error), file=sys.stderr)
        raise SystemExit(1)

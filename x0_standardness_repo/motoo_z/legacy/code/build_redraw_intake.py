#!/usr/bin/env python3
"""Package the redraw visual intake and tracing work; no new topology claim."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = args.out.resolve()
    require(not target.exists(), 'Choose a new archive output')
    selected = set()
    for rel in [
        'REDRAW_INTAKE.md', 'requirements.txt',
        'inputs/motoo_rollspun.pdf', 'inputs/figure1_link.json', 'inputs/figure2_link.json',
        'inputs/redraw_20261102', 'results/redraw_20261102',
        'scratch/redraw_20261102/centerline_dev',
        'code/extract_redraw.py', 'code/trace_redraw_centerlines.py',
        'code/check_redraw_markings.py', 'code/check_archive_hashes.py',
        'code/build_redraw_intake.py',
    ]:
        item = root / rel
        require(item.exists(), 'Missing selected file: ' + rel)
        files = item.rglob('*') if item.is_dir() else [item]
        selected.update(p for p in files if p.is_file() and p.suffix != '.pyc'
                        and '__pycache__' not in p.parts)
    status = json.loads((root / 'results/redraw_20261102/INTAKE_STATUS.json').read_text())
    require(status['status'] == 'VISUALLY_SUITABLE_CHECKED_PD_PENDING', 'Unexpected intake status')
    require(not any(status[k] for k in ('independent_isotopy_verified', 'new_pd_generated',
                'native_katie_admission_completed', 'new_redraw_group_or_boundary_computed',
                'new_redraw_4d_triangulation_generated')), 'Package scope unexpectedly changed')
    for rel, digest in status['bindings'].items():
        require(root / rel in selected and sha(root / rel) == digest, 'Receipt binding mismatch: ' + rel)
    source = json.loads((root / 'inputs/redraw_20261102/SOURCE_RECEIPT.json').read_text())
    for item in source['files']:
        require(sha(root / item['file']) == item['sha256'], 'Uploaded PDF bytes changed')
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='motoo_redraw_archive_') as temporary:
        parent = Path(temporary)
        stage = parent / 'motoo_redraw_intake_v1'
        manifest = []
        for original in sorted(selected):
            rel = original.relative_to(root)
            destination = stage / rel
            destination.parent.mkdir(parents=True, exist_ok=True)
            digest = sha(original)
            shutil.copy2(original, destination)
            require(sha(destination) == digest == sha(original), 'Input changed during staging')
            manifest.append(f'{digest}  {rel.as_posix()}')
        (stage / 'SHA256SUMS').write_text('\n'.join(manifest) + '\n')
        pending = parent / 'bundle.zip'
        with zipfile.ZipFile(pending, 'w', compression=zipfile.ZIP_DEFLATED,
                             compresslevel=9, allowZip64=True) as z:
            for path in sorted(stage.rglob('*')):
                if path.is_file():
                    z.write(path, path.relative_to(parent).as_posix())
        extracted = parent / 'extracted'
        with zipfile.ZipFile(pending) as z:
            require(z.testzip() is None, 'ZIP CRC failed')
            require(len(z.namelist()) == len(selected) + 1, 'Wrong archive inventory')
            z.extractall(extracted)
        project = extracted / stage.name
        subprocess.run([sys.executable, '-O', str(project / 'code/check_archive_hashes.py')], check=True)
        actual = {p.relative_to(project).as_posix() for p in project.rglob('*') if p.is_file()}
        expected = {p.relative_to(root).as_posix() for p in selected} | {'SHA256SUMS'}
        require(actual == expected, 'Extracted inventory differs')
        shutil.copy2(pending, target)
    receipt = {'file': str(target), 'bytes': target.stat().st_size, 'sha256': sha(target),
               'files_in_manifest': len(selected), 'zip_crc_passed': True,
               'extracted_file_hashes_passed': True, 'extracted_inventory_passed': True,
               'status': status['status'], 'new_4d_triangulation_included': False}
    target.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

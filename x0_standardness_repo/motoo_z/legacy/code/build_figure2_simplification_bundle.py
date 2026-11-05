#!/usr/bin/env python3
"""Create and byte-check a frozen Figure2 handoff from an explicit file list."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def build(root, selection, output):
    root, selection, output = Path(root).resolve(), Path(selection), Path(output)
    data = json.loads(selection.read_text())
    paths = set()
    for item in data['include']:
        path = root / item
        if not path.exists():
            raise FileNotFoundError(path)
        for candidate in (path.rglob('*') if path.is_dir() else [path]):
            if not candidate.is_file():
                continue
            relative = candidate.relative_to(root)
            if '__pycache__' in relative.parts or candidate.suffix in ('.pyc', '.part', '.tmp'):
                continue
            if any(part.startswith('.') for part in relative.parts):
                continue
            paths.add(relative)
    prefix = data.get('archive_directory', 'motoo_figure2_simplification_v1')
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = {'schema': 'motoo.figure2.bundle-manifest.v1', 'files': [],
                'scope': data['scope'], 'archive_directory': prefix}
    renames = data.get('rename', {})
    written_names = set()
    temporary = output.with_suffix(output.suffix + '.tmp')
    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED,
                         compresslevel=9, allowZip64=True) as archive:
        for relative in sorted(paths):
            path = root / relative
            destination = Path(renames.get(relative.as_posix(), relative.as_posix()))
            if destination.is_absolute() or '..' in destination.parts or destination.as_posix() in written_names:
                raise ValueError('Unsafe or duplicate archive path: ' + str(destination))
            written_names.add(destination.as_posix())
            before = path.stat()
            content = path.read_bytes()
            after = path.stat()
            if before.st_mtime_ns != after.st_mtime_ns or before.st_size != after.st_size:
                raise RuntimeError('Selected file changed during packaging: ' + str(relative))
            archive.writestr(prefix + '/' + destination.as_posix(), content)
            manifest['files'].append({'path': destination.as_posix(), 'bytes': len(content),
                                      'sha256': hashlib.sha256(content).hexdigest()})
        archive.writestr(prefix + '/MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    with zipfile.ZipFile(temporary) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise RuntimeError('ZIP CRC failure: ' + bad)
        for item in manifest['files']:
            content = archive.read(prefix + '/' + item['path'])
            if len(content) != item['bytes'] or hashlib.sha256(content).hexdigest() != item['sha256']:
                raise RuntimeError('Archived file differs: ' + item['path'])
    temporary.replace(output)
    receipt = {'status': 'PASS_ZIP_CRC_AND_EVERY_FILE_HASH', 'file': str(output.resolve()),
               'files': len(manifest['files']), 'bytes': output.stat().st_size,
               'sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
               'uncompressed_bytes': sum(x['bytes'] for x in manifest['files'])}
    output.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=Path(__file__).resolve().parents[1], type=Path)
    parser.add_argument('--selection', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    build(args.root, args.selection, args.out)

#!/usr/bin/env python3
"""Integrity and file-interface checks for the X0 reduction module."""
from pathlib import Path
import gzip
import hashlib
import json
import os


ROOT = Path(__file__).resolve().parent
MASTER_SHA256 = '251d7ad38235869fa268e9ee6ec9da2f397cf1ca71db267d0bf0ca9b434c5e4c'
COVER_SHA256 = '798c4c5cc26571d3b3cf05655ce800bd4a62b2e68a9ac82e1a3d8db24479be98'
FACETS_SHA256 = 'f40b43c09373dff94cde542a3076a8c2f482b1632726b63eed2c5610ddd1e92b'
CERTIFICATE_SHA256 = 'e9ba9f14831bb0da32043d7d93eb0148f0683e2a6e4f67d5156e8e9fc94adcdf'
NATIVE_SOURCE_SHA256 = 'bc1f378af75c8cd85a2be9ac8966aa6547abc8e760248af2765b629f60e34512'
ENGINE_SHA256 = '42e0c5a630d71f98c6e5ff3c70c520d1c3d136151aad4ac37fe807004fbaca1d'
G6 = 'gAMLMQaabcdeffff2a+a2avbmaoaxaOb0aUa'
S2 = 'cPkbbbbaaaaaaaa'


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def file_record(path):
    p = Path(path)
    return {'bytes': p.stat().st_size, 'sha256': sha256(p)}


def read_json(path):
    def no_duplicates(pairs):
        out = {}
        for key, value in pairs:
            require(key not in out, 'duplicate JSON key: ' + key)
            out[key] = value
        return out
    path = Path(path)
    raw = path.read_bytes()
    if path.suffix == '.gz':
        raw = gzip.decompress(raw)
    return json.loads(raw, object_pairs_hook=no_duplicates)


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def absolute_no_symlinks(path):
    """Keep lexical ancestors visible while rejecting symlinked input trees."""
    p = Path(os.path.abspath(Path(path).expanduser()))
    for q in (p, *p.parents):
        require(not q.is_symlink(), 'symlinked path is not accepted: ' + str(q))
    return p


def safe_relative(root, name):
    rel = Path(name)
    require(not rel.is_absolute() and rel.parts and '..' not in rel.parts,
            'unsafe inventory path: ' + str(name))
    p = absolute_no_symlinks(root / rel)
    require(p.is_relative_to(root), 'inventory path escapes module: ' + str(name))
    return p


def input_integrity():
    """Check the exact imported proof payload, including its original manifests."""
    imported = {}
    inventories = {}
    for name, schema in (
        ('provenance/master_sources.json', 'x0-reduction-imported-master-files-v1'),
        ('provenance/terminal_sources.json', 'x0-reduction-imported-terminal-files-v1'),
    ):
        p = safe_relative(ROOT, name)
        inventory = read_json(p)
        require(inventory['schema'] == schema, 'wrong imported-source inventory schema')
        inventories[name] = file_record(p)
        for rel, entry in inventory['files'].items():
            require(rel not in imported, 'duplicate imported file: ' + rel)
            q = safe_relative(ROOT, rel)
            require(q.is_file(), 'missing imported file: ' + rel)
            actual = file_record(q)
            require(actual == {'bytes': entry['bytes'], 'sha256': entry['sha256']},
                    'imported file changed: ' + rel)
            imported[rel] = actual
    # These primary proof identities are explicit even if an inventory is edited.
    for rel, expected in (
        ('reference/original/candidate/marked_candidate.json.gz', MASTER_SHA256),
        ('reference/original/candidate/cover_regina.json.gz', COVER_SHA256),
        ('elementary/step3_macro/x0_C1332_regina.json.gz', NATIVE_SOURCE_SHA256),
        ('certificates/standardization.json', CERTIFICATE_SHA256),
    ):
        require(sha256(safe_relative(ROOT, rel)) == expected, 'primary proof input changed: ' + rel)
    return {'inventories': inventories, 'imported_files': imported,
            'imported_file_count': len(imported)}


def code_integrity():
    names = ['verify.py', 'common.py', 'terminal.py']
    names += [str(p.relative_to(ROOT)) for p in sorted((ROOT / 'source').rglob('*.py'))]
    names += [str(p.relative_to(ROOT)) for p in sorted((ROOT / 'elementary').rglob('*.py'))]
    return {name: file_record(safe_relative(ROOT, name)) for name in names}


def candidate_identity(directory):
    """Bind the exact candidate; the fixed certificate is not a generic reducer."""
    directory = absolute_no_symlinks(directory)
    require(directory.is_dir(), 'candidate directory does not exist: ' + str(directory))
    records = {}
    for name, expected in (('marked_candidate.json.gz', MASTER_SHA256),
                           ('cover_regina.json.gz', COVER_SHA256)):
        p = absolute_no_symlinks(directory / name)
        require(p.is_file(), 'candidate input is missing: ' + str(p))
        record = file_record(p)
        require(record['sha256'] == expected,
                name + ' is not the exact released construction output; expected SHA-256 ' + expected)
        records[name] = record
    master = read_json(directory / 'marked_candidate.json.gz')
    require(len(master['facets']) == 103744 and len(master['vertices']) == 4310,
            'candidate has unexpected size')
    require(canonical_sha256(master['facets']) == master['facets_sha256'] == FACETS_SHA256,
            'candidate facet identity mismatch')
    return {'directory': str(directory), 'files': records,
            'facets_sha256': FACETS_SHA256, 'pentachora': 103744, 'vertices': 4310}


def require_new_output(path, candidate=None):
    out = absolute_no_symlinks(path)
    require(not out.exists(), '--output must name a new directory: ' + str(out))
    require(out != ROOT and out not in ROOT.parents, 'output contains the module')
    for name in ('source', 'elementary', 'reference', 'certificates', 'provenance',
                 'verification', 'environment', 'docs'):
        protected = ROOT / name
        require(out != protected and protected not in out.parents,
                'output is inside protected module inputs: ' + str(out))
    if candidate is not None:
        candidate = absolute_no_symlinks(candidate)
        require(out != candidate and candidate not in out.parents and out not in candidate.parents,
                'output and selected candidate directories must be disjoint')
    return out


def check_release_manifest():
    p = ROOT / 'MANIFEST_SHA256.json'
    require(p.is_file(), 'release manifest is missing')
    manifest = read_json(p)
    require(manifest['schema'] == 'x0-reduction-release-manifest-v1', 'release manifest schema')
    for name, record in manifest['files'].items():
        q = safe_relative(ROOT, name)
        require(q.is_file() and file_record(q) == record, 'release payload mismatch: ' + name)
    return {'manifest_sha256': sha256(p), 'payload_files_checked': len(manifest['files'])}


if __name__ == '__main__':
    imported = input_integrity()
    release = check_release_manifest()
    print(json.dumps({'status': 'INTEGRITY_PASS',
                      'imported_files_checked': imported['imported_file_count'],
                      **release}, indent=2))

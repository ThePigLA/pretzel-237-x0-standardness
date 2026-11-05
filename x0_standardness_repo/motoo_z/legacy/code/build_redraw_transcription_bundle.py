#!/usr/bin/env python3
"""Package the checked redraw data, preserving explicit failed/conditional gates."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
PREFIX='motoo_redraw_transcription_v1'
OUT=ROOT/'output/motoo_redraw_transcription_v1.zip'
FILES={}

def add(path, alias=None):
    p=ROOT/path
    if p.is_dir():
        for f in sorted(p.rglob('*')):
            if f.is_file() and '__pycache__' not in f.parts and f.suffix not in ('.pyc','.lock'):
                add(f.relative_to(ROOT))
    elif p.is_file():
        FILES[str(alias or path)]=p
    else:
        raise FileNotFoundError(p)

add('REDRAW_TRANSCRIPTION.md','README.md')
for p in [
    'inputs/motoo_rollspun.pdf','inputs/figure1_link.json','inputs/figure2_link.json',
    'inputs/redraw_20261102','source',
    'results/redraw_20261102/extraction_dev/figure1',
    'results/redraw_20261102/extraction_dev/figure2',
    'results/redraw_20261102/TRANSCRIPTION_REPLAY.json',
    'results/redraw_20261102/TRANSCRIPTION_STATUS.json',
    'results/redraw_20261102/katie_figure1_source_v1',
    'results/redraw_20261102/katie_figure2_dashed_v1',
    'scratch/redraw_20261102/assembled_f1_v1',
    'scratch/redraw_20261102/assembled_f2_v1',
    'scratch/redraw_20261102/assembled_f2_dashed_v1',
    'scratch/redraw_20261102/reconstruction_dev/figure1/auto_v2',
    'scratch/redraw_20261102/reconstruction_dev/figure1/auto_v3',
    'scratch/redraw_20261102/reconstruction_dev/figure1/selfcrossing_review',
    'scratch/redraw_20261102/reconstruction_dev/figure1/cyan_marking_audit',
    'scratch/redraw_20261102/reconstruction_dev/figure1/cyan_join_audit',
    'scratch/redraw_20261102/reconstruction_dev/figure1/SELF_CROSSING_AUDIT.json',
    'scratch/redraw_20261102/reconstruction_dev/figure1/CYAN_SOURCE_AUDIT.json',
    'scratch/redraw_20261102/manual_f2_dev',
    'scratch/redraw_20261102/crossings_f1_v2',
    'scratch/redraw_20261102/crossings_f1_v4',
    'scratch/redraw_20261102/crossings_f2_v1',
    'scratch/redraw_20261102/crossings_f2_v2',
    'scratch/redraw_20261102/f1_v4_review',
    'output/figure2_crossing_question.png',
]: add(p)
for p in sorted((ROOT/'code').glob('*.py')):add(p.relative_to(ROOT))
for pattern in ['manual_f1*.json','manual_f2*.json','f1_overrides*.json','f2_overrides*.json',
                'figure1_cyan_independent_audit.json','trace_f1_dots_red.py','trace_f2_dots.py']:
    for p in sorted((ROOT/'scratch/redraw_20261102').glob(pattern)):add(p.relative_to(ROOT))
completion=ROOT/'results/redraw_20261102/figure2_completion_conditional_v1'
if completion.is_dir():
    add(completion.relative_to(ROOT))
    for directory in ['results/boundary/standard_cap_six_boundary_v1',
                      'results/boundary/upper_cap_selftest_v1']:
        add(directory)

manifest={'schema':'motoo.distribution-manifest.v1','archive_root':PREFIX,
          'scope':'Checked redraw transcriptions, explicit identification discrepancy, and conditional Figure2 conversion; no new standardness claim.',
          'files':[{'path':name,'size_bytes':path.stat().st_size,
                    'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
                   for name,path in sorted(FILES.items())]}
manifest_bytes=(json.dumps(manifest,indent=2)+'\n').encode()
OUT.parent.mkdir(exist_ok=True);tmp=OUT.with_suffix('.tmp.zip')
with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for name,path in sorted(FILES.items()):z.writestr(PREFIX+'/'+name,path.read_bytes())
    z.writestr(PREFIX+'/MANIFEST.json',manifest_bytes)
with zipfile.ZipFile(tmp) as z:
    if z.testzip() is not None:raise ValueError('ZIP CRC verification failed')
    for row in manifest['files']:
        b=z.read(PREFIX+'/'+row['path'])
        if len(b)!=row['size_bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:
            raise ValueError('Archived byte verification failed: '+row['path'])
    if len(z.namelist())!=len(manifest['files'])+1:raise ValueError('Unexpected archive entries')
tmp.replace(OUT)
receipt={'archive':str(OUT),'size_bytes':OUT.stat().st_size,
         'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest(),
         'files':len(manifest['files']),'zip_crc_passed':True,'every_archived_sha256_checked':True}
(ROOT/'output/motoo_redraw_transcription_v1_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))

#!/usr/bin/env python3
"""Freeze a self-contained source/certificate archive, omitting the local venv."""
from pathlib import Path
import hashlib,json,os,zipfile
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'pretzel_237_even_control_v1_20261026.zip'
payload={}
for folder,dirs,files in os.walk(ROOT):
    dirs[:]=sorted(d for d in dirs if d not in ('venv','.venv','__pycache__','.git'))
    for name in sorted(files):
        p=Path(folder)/name;rel=p.relative_to(ROOT).as_posix()
        if name=='MANIFEST_SHA256.json' or name.endswith(('.pyc','.tmp')):continue
        payload[rel]=p.read_bytes()
manifest={'schema':'exact-release-file-manifest-v1','files':{p:{'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()} for p,b in sorted(payload.items())}}
m=(json.dumps(manifest,sort_keys=True,indent=2)+'\n').encode();(ROOT/'MANIFEST_SHA256.json').write_bytes(m)
payload['MANIFEST_SHA256.json']=m
with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p,b in sorted(payload.items()):z.writestr('pretzel_237_even_control_v1/'+p,b)
print(json.dumps({'path':str(OUT),'files':len(payload),'bytes':OUT.stat().st_size,'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}))

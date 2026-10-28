#!/usr/bin/env python3
"""Native finite source admission; simplification uses only an unmarked copy."""
from pathlib import Path
import json,gzip,time,hashlib,regina

ROOT=Path(__file__).resolve().parent

def load(p):return json.loads(gzip.decompress(p.read_bytes()))

def build(data):
    t=regina.Triangulation3();t.newTetrahedra(data['simplex_count'])
    for a,f,b,p in data['gluings']:
        t.tetrahedron(a).join(f,t.tetrahedron(b),regina.Perm4(*p))
    return t

def main():
    out={'scope':'Native admission of explicit facet models; original marked source unchanged','regina':regina.versionString(),'cases':{}}
    for name in ('quotient','cover'):
        t0=time.monotonic();p=ROOT/'data'/f'{name}_regina.json.gz';raw=load(p);t=build(raw)
        row={'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'tetrahedra':t.size(),
             'valid':t.isValid(),'closed':t.isClosed(),'orientable':t.isOrientable(),'connected':t.isConnected(),
             'f_vector':list(t.fVector())}
        if not all(row[k] for k in ('valid','closed','orientable','connected')):raise ValueError(row)
        print(name,'admitted',json.dumps(row),flush=True)
        c=regina.Triangulation3(t);c.simplify();row['unmarked_simplified_size']=c.size();row['unmarked_simplified_isoSig']=c.isoSig()
        print(name,'simplified',c.size(),flush=True)
        row['H1_on_unmarked_simplified_copy']=str(c.homology())
        if name=='quotient':row['S3_recognition']=c.isSphere()
        row['seconds']=time.monotonic()-t0;out['cases'][name]=row
        print(name,json.dumps(row),flush=True)
        (ROOT/'data'/'native_admission.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()

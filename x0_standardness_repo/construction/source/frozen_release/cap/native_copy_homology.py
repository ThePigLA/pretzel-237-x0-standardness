#!/usr/bin/env python3
"""Bounded invariant cross-check on an unmarked simplification copy.

This is not a replayable standardization certificate, and does not modify the
source cap.  Each native child has a wall cap and writes stage progress.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, subprocess, sys, time


def child(path,output):
    import regina
    raw=path.read_bytes()
    with gzip.open(path,'rt') as f:c=json.load(f)
    r={'schema':'cap-unmarked-copy-invariants-v1','twist':c['twist'],'source_file_sha256':hashlib.sha256(raw).hexdigest(),'scope':'invariant cross-check; simplify is trusted, no portable move certificate claimed'}
    def save(stage):r['stage']=stage;output.write_text(json.dumps(r,indent=2)+'\n')
    save('constructing')
    tr=regina.Triangulation4();tr.newPentachora(len(c['facets']))
    for i,k,j,p in c['gluings']:tr.pentachoron(i).join(k,tr.pentachoron(j),regina.Perm5(*p))
    r['source_size']=tr.size();r['source_valid']=tr.isValid();r['source_euler']=tr.eulerCharTri();save('simplifying_copy')
    cp=regina.Triangulation4(tr);changed=cp.simplify();r['copy_size']=cp.size();r['copy_simplified']=changed;save('copy_invariants')
    r.update({'H1':str(cp.homology()),'H2':str(cp.homology(2)),'euler':cp.eulerCharTri(),'copy_valid':cp.isValid(),'copy_boundary_components':cp.countBoundaryComponents()})
    if r['H1']!='0' or r['H2']!='Z' or r['euler']!=2 or not r['copy_valid'] or r['copy_boundary_components']!=1:raise ValueError(r)
    r['result']='PASS';save('completed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--child',type=Path);p.add_argument('--output',type=Path);p.add_argument('--seconds',type=int,default=60);a=p.parse_args()
    if a.child:child(a.child,a.output)
    else:
        root=Path(__file__).resolve().parent;results=[]
        for m in (0,2):
            op=root/f'native_copy_m{m}.json';start=time.monotonic()
            try:
                r=subprocess.run([sys.executable,__file__,'--child',str(root/'data'/f'cap_m{m}.json.gz'),'--output',str(op)],timeout=a.seconds,check=False,capture_output=True,text=True)
                result={'twist':m,'returncode':r.returncode,'seconds':time.monotonic()-start,'receipt':op.name}
                if r.returncode:result['stderr']=r.stderr[-2000:]
            except subprocess.TimeoutExpired:result={'twist':m,'result':'WALL_CAP','seconds':time.monotonic()-start,'receipt':op.name}
            results.append(result);print(json.dumps(result),flush=True)
        (root/'native_copy_summary.json').write_text(json.dumps(results,indent=2)+'\n')

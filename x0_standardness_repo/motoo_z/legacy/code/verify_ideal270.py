#!/usr/bin/env python3
"""Portable binding check plus independent replay of the ideal-270 certificate."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import regina
from audit_figure2_snaps import audit


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(source,certificate,out):
    source,certificate,out=map(Path,(source,certificate,out))
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina 7.4.1')
    md_path=certificate/'SOURCE.json';final_path=certificate/'FINAL.json'
    md=json.loads(md_path.read_text());final=json.loads(final_path.read_text())
    require(final['status']=='PASS_EXACT_DETERMINISTIC_TARGET','certificate did not reach its target')
    log=certificate/'moves.jsonl';endpoint=certificate/final['endpoint_file']
    bindings={p:sha(p) for p in (source,md_path,final_path,log,endpoint)}
    require(bindings[source]==md['source_sha256']==final['source_sha256'],'wrong source RGA')
    require(bindings[log]==final['log_sha256'],'move log checksum differs')
    require(bindings[endpoint]==final['endpoint_rga_sha256'],'endpoint RGA checksum differs')
    require(out.resolve() not in {p.resolve() for p in bindings},'output would overwrite an input')
    out.parent.mkdir(parents=True,exist_ok=True)
    temporary=out.with_name(out.name+'.independent.tmp')
    audit(source,log,temporary)
    result=json.loads(temporary.read_text())
    require(result['moves']==final['moves'],'replay move count differs')
    require(result['final_tight_sha256']==final['current_tight_sha256'],'replay endpoint differs')
    require(result['final_pentachora']==final['current_pentachora']==270,'unexpected final size')
    packet=regina.open(str(endpoint))
    require(isinstance(packet,regina.Triangulation4),'endpoint is not a 4D triangulation')
    require(hashlib.sha256(packet.tightEncoding().encode()).hexdigest()==result['final_tight_sha256'],
            'saved endpoint does not match replay')
    require(all(sha(p)==value for p,value in bindings.items()),'a certificate input changed during replay')
    result['certificate_binding_passed']=True
    result['endpoint_file_bound_to_replay']=True
    result['regina_distribution']='7.4.1'
    result['verifier_sha256']=sha(__file__)
    result['source_path_supplied_explicitly']=True
    temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(out)
    print(json.dumps({'status':result['status'],'certificate_binding_passed':True,
        'moves':result['moves'],'pentachora':result['final_pentachora'],
        'final_tight_sha256':result['final_tight_sha256'],'report':str(out)}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--certificate',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args();verify(a.source,a.certificate,a.out)

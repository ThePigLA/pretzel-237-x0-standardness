#!/usr/bin/env python3
"""Independently replay immutable ideal-search chunks while production runs.

Retains the live triangulation between sealed manifests and reconstructs
every edge-snap gadget independently. Only closed immutable chunks are read.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import regina
from audit_figure2_snaps import explicit_snap_gadget,require,digest


def sha_file(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save_json(path,data):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def audit(directory,source,out,seconds):
    directory,source,out=Path(directory),Path(source),Path(out)
    out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    meta=json.loads((directory/'SOURCE.json').read_text())
    policy=json.loads((directory/'POLICY.json').read_text())
    require(sha_file(source)==meta['source_sha256'],'source checksum mismatch')
    tri=regina.Triangulation4(regina.open(str(source)))
    require(tri.isValid() and tri.isConnected() and tri.isOrientable(),'raw ideal source admission failed')
    require(digest(tri)==meta['start_tight_sha256'],'source labelled-state mismatch')
    seen=[];n=0;counts=Counter();snap_checks=[];last_print=started
    target=policy['target_moves'];expected=policy['target_tight_sha256']
    while time.monotonic()-started<seconds:
        manifest=json.loads((directory/'SEALED.json').read_text())
        chunks=manifest['chunks']
        require(chunks[:len(seen)]==seen,'previously sealed chunk list changed')
        for chunk in chunks[len(seen):]:
            require(chunk['start']==n,'chunk sequence gap')
            path=directory/chunk['file']
            require(sha_file(path)==chunk['sha256'],'sealed chunk checksum mismatch')
            rows=path.read_text().splitlines()
            require(len(rows)==chunk['end']-chunk['start'],'chunk line count mismatch')
            for line in rows:
                row=json.loads(line);operation=row['operation'];kind=operation[0]
                require(row['index']==n and row['before_sha256']==digest(tri),'broken state chain')
                before=tri.size()
                if kind=='S':
                    index=operation[-1];edge=tri.edge(index)
                    count=tri.countVertices()
                    endpoints=[{'index':edge.vertex(i).index(),'valid':edge.vertex(i).isValid(),
                                'boundary':edge.vertex(i).isBoundary(),'ideal':edge.vertex(i).isIdeal()} for i in (0,1)]
                    gadget=explicit_snap_gadget(tri,index)
                    require(tri.hasSnapEdge(edge) and tri.snapEdge(edge),'snap guard rejected')
                    require(tri.tightEncoding()==gadget.tightEncoding(),'independent snap gadget differs')
                    require(tri.size()==before+4 and tri.countVertices()==count-1,'snap size/vertex counts differ')
                    snap_checks.append({'index':n,'edge':index,'endpoints':endpoints,'explicit_gadget_matches':True})
                else:
                    require(len(operation)==3,'malformed move')
                    _,dim,index=operation
                    face=tri.simplex(index) if dim==4 else tri.face(dim,index)
                    if kind=='P':
                        require(tri.hasPachner(face) and tri.pachner(face),'Pachner guard rejected')
                        require(tri.size()==before+2*dim-4,'Pachner size differs')
                    elif kind=='Z':
                        require(dim in (0,1,2) and tri.has20(face) and tri.move20(face),'2-0 guard rejected')
                        require(tri.size()==before-2,'2-0 size differs')
                    else:raise ValueError('unsupported move')
                require(digest(tri)==row['after_sha256'] and tri.size()==row['pentachora'],'replayed state differs')
                n+=1;counts['S' if kind=='S' else str(operation[:2])]+=1
            require(n==chunk['end'],'chunk endpoint count mismatch')
            seen.append(chunk)
            if time.monotonic()-last_print>10:
                progress={'status':'VERIFIED_SEALED_PREFIX','moves':n,'pentachora':tri.size(),'tight_sha256':digest(tri),'seconds':time.monotonic()-started}
                save_json(out/'PROGRESS.json',progress);print(json.dumps(progress),flush=True);last_print=time.monotonic()
        require(n==manifest['moves'] and digest(tri)==manifest['current_tight_sha256'],'sealed manifest state differs')
        if n==target:
            require(digest(tri)==expected,'target state differs')
            final_path=directory/'FINAL.json'
            if final_path.exists():
                final=json.loads(final_path.read_text())
                require(final['moves']==n and final['current_tight_sha256']==expected,'FINAL state differs')
                producer_endpoint=directory/final['endpoint_file']
                require(sha_file(producer_endpoint)==final['endpoint_rga_sha256'],'producer RGA checksum differs')
                require(regina.Triangulation4(regina.open(str(producer_endpoint))).tightEncoding()==tri.tightEncoding(),'producer RGA differs from independent replay')
                require(sha_file(directory/'moves.jsonl')==final['log_sha256'],'assembled producer log checksum differs')
                endpoint=out/'replayed_endpoint.rga'
                require(regina.PacketOfTriangulation4(tri).save(str(endpoint)),'independent endpoint save failed')
                report={'status':'PASS_SEALED_IDEAL_CHAIN_AND_EXPLICIT_SNAP_GADGETS','source_sha256':sha_file(source),'moves':n,'counts':dict(counts),'snap_checks':snap_checks,'chunks':seen,'end_pentachora':tri.size(),'end_vertices':tri.countVertices(),'end_valid':tri.isValid(),'end_ideal':tri.isIdeal(),'end_boundary_components':tri.countBoundaryComponents(),'end_real_boundary_facets':tri.countBoundaryFacets(),'end_tight_sha256':digest(tri),'independent_endpoint_sha256':sha_file(endpoint),'producer_FINAL_sha256':sha_file(final_path),'seconds':time.monotonic()-started,'standardness_asserted':False,'closed_completion_asserted':False,'trust_boundary':'Regina checked move kernel and source recognition; independent explicit snap gluing and driver'}
                save_json(out/'AUDIT.json',report)
                print(json.dumps({k:v for k,v in report.items() if k not in ('snap_checks','chunks')}),flush=True)
                return report
        require(n<=target,'producer passed expected target')
        time.sleep(.25)
    raise TimeoutError('streaming audit wall budget expired; verified prefix report remains')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seconds',type=float,default=900)
    a=p.parse_args();audit(a.directory,a.source,a.out,a.seconds)

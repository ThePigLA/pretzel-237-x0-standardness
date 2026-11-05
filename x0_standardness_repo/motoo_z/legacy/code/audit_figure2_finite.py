#!/usr/bin/env python3
"""Independent-driver replay of a frozen finite-move checkpoint.

Does not import the search/replay producer. Uses Regina's guarded move kernel
and recognition, so it is not an independent implementation of PL moves.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import regina


def require(ok,message):
    if not ok:raise ValueError(message)


def sha_file(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def sha_state(tri):return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()


def audit(directory,out,source_override=None):
    start=time.monotonic();directory=Path(directory)
    meta=json.loads((directory/'SOURCE.json').read_text())
    final=json.loads((directory/'CHECKPOINT.json').read_text())
    source=Path(source_override) if source_override is not None else Path(meta['source_file'])
    require(sha_file(source)==meta['source_sha256'],'source checksum mismatch')
    tri=regina.Triangulation4(regina.open(str(source)))
    require(tri.isValid() and tri.isConnected() and tri.isOrientable() and not tri.isIdeal(),'source finite-manifold admission failed')
    source_properties={'pentachora':tri.size(),'vertices':tri.countVertices(),'real_boundary_facets':tri.countBoundaryFacets(),'boundary_components':tri.countBoundaryComponents()}
    require(sha_state(tri)==meta['start_tight_sha256'],'labelled source mismatch')
    print(json.dumps({'admitted_source':source_properties,'seconds':time.monotonic()-start}),flush=True)
    limit=final['moves'];counts=Counter();checked=0
    for line in (directory/'moves.jsonl').read_text().splitlines():
        if checked==limit:break
        row=json.loads(line);kind,dim,index=row['operation']
        require(row['index']==checked,'move numbering mismatch')
        require(sha_state(tri)==row['before_sha256'],'before hash mismatch')
        require(type(dim) is int and dim in range(5),'invalid face dimension')
        require(type(index) is int and 0<=index<(tri.size() if dim==4 else tri.countFaces(dim)),'invalid face index')
        face=tri.simplex(index) if dim==4 else tri.face(dim,index)
        size=tri.size();vertices=tri.countVertices()
        if kind=='C':
            require(dim==1,'collapse dimension is not one')
            require(face.vertex(0)!=face.vertex(1),'collapse endpoints identical')
            require(tri.hasCollapseEdge(face),'edge-collapse guard rejected')
            require(tri.collapseEdge(face),'checked edge collapse failed')
            require(tri.size()<size and tri.countVertices()==vertices-1,'collapse count changes wrong')
        elif kind=='S':
            require(dim==4,'shelling dimension is not four')
            require(tri.hasShellBoundary(face),'boundary-shelling guard rejected')
            require(tri.shellBoundary(face),'checked boundary shelling failed')
            require(tri.size()==size-1,'shelling count change wrong')
        elif kind=='P':
            require(tri.hasPachner(face),'Pachner guard rejected')
            require(tri.pachner(face),'checked Pachner move failed')
            require(tri.size()==size+2*dim-4,'Pachner count change wrong')
        elif kind=='Z':
            require(dim in (0,1,2),'bad 2-0 face dimension')
            require(tri.has20(face),'2-0 guard rejected')
            require(tri.move20(face),'checked 2-0 failed')
            require(tri.size()==size-2,'2-0 count change wrong')
        elif kind=='F':
            require(dim==1 and tri.has44(face),'4-4 guard rejected')
            require(tri.move44(face) and tri.size()==size,'checked 4-4 failed')
        else:raise ValueError('unsupported move code')
        require(sha_state(tri)==row['after_sha256'] and tri.size()==row['pentachora'],'result state hash or size mismatch')
        checked+=1;counts[kind]+=1
        if checked%10==0:print(json.dumps({'replayed':checked,'pentachora':tri.size(),'seconds':time.monotonic()-start}),flush=True)
    require(checked==limit,'incomplete log')
    require(sha_state(tri)==final['current_tight_sha256'],'checkpoint endpoint hash mismatch')
    saved=regina.Triangulation4(regina.open(str(directory/'current.rga')))
    require(saved.tightEncoding()==tri.tightEncoding(),'saved endpoint differs')
    require(tri.isValid() and tri.isConnected() and tri.isOrientable() and not tri.isIdeal(),'endpoint admission failed')
    result={'status':'PASS_INDEPENDENT_DRIVER_FINITE_REPLAY','directory':str(directory.resolve()),'source_sha256':sha_file(source),'move_log_sha256':sha_file(directory/'moves.jsonl'),'source_properties':source_properties,'moves':checked,'counts':dict(counts),'endpoint_pentachora':tri.size(),'endpoint_vertices':tri.countVertices(),'endpoint_real_boundary_facets':tri.countBoundaryFacets(),'endpoint_boundary_components':tri.countBoundaryComponents(),'endpoint_tight_sha256':sha_state(tri),'seconds':time.monotonic()-start,'standardness_asserted':False,'trust_boundary':'Regina guarded move kernel and finite-manifold admission; independent driver and exact state replay'}
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--source',type=Path,help='relocated original source; certificate SHA-256 remains enforced')
    a=p.parse_args();audit(a.directory,a.out,a.source)

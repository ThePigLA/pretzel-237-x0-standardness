#!/usr/bin/env python3
"""Freeze a committed finite prefix by independent replay, then export boundary.

Copies only atomic JSON metadata and committed move-log prefixes from a live
run. Regenerates the lower-body RGA by replay; never copies a live endpoint.
Keeps the admitted triangulation alive through boundary extraction so no
cached topological fact is asserted or manually manufactured.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import time
import regina
from complete_upper_handles import boundary_data


def require(ok,message):
    if not ok:raise ValueError(message)


def file_sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def state_sha(tri):return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()
def save_json(path,data):Path(path).write_text(json.dumps(data,indent=2)+'\n')


def guarded(tri,op):
    require(len(op)==3,'bad move length')
    kind,dim,index=op
    require(type(dim) is int and dim in range(5),'bad move dimension')
    require(type(index) is int and 0<=index<(tri.size() if dim==4 else tri.countFaces(dim)),'bad move index')
    face=tri.simplex(index) if dim==4 else tri.face(dim,index)
    before=tri.size()
    if kind=='C':
        require(dim==1 and face.vertex(0)!=face.vertex(1),'bad collapse endpoints')
        vertices=tri.countVertices()
        require(tri.hasCollapseEdge(face),'collapse guard rejected')
        require(tri.collapseEdge(face),'checked collapse failed')
        require(tri.size()<before and tri.countVertices()==vertices-1,'wrong collapse counts')
    elif kind=='S':
        require(dim==4 and tri.hasShellBoundary(face),'shell guard rejected')
        require(tri.shellBoundary(face) and tri.size()==before-1,'checked shell failed')
    elif kind=='P':
        require(tri.hasPachner(face),'Pachner guard rejected')
        require(tri.pachner(face) and tri.size()==before+2*dim-4,'checked Pachner failed')
    elif kind=='Z':
        require(dim in (0,1,2) and tri.has20(face),'2-0 guard rejected')
        require(tri.move20(face) and tri.size()==before-2,'checked 2-0 failed')
    elif kind=='F':
        require(dim==1 and tri.has44(face),'4-4 guard rejected')
        require(tri.move44(face) and tri.size()==before,'checked 4-4 failed')
    else:raise ValueError('unknown move code')


def committed_snapshot(directory,out,prefix):
    """Read a single atomic checkpoint, then copy exactly its committed log."""
    directory=Path(directory)
    metadata_bytes=(directory/'SOURCE.json').read_bytes()
    checkpoint_bytes=(directory/'CHECKPOINT.json').read_bytes()
    meta=json.loads(metadata_bytes);check=json.loads(checkpoint_bytes)
    limit=check['moves'];lines=[]
    with (directory/'moves.jsonl').open() as stream:
        for line in stream:
            if len(lines)==limit:break
            lines.append(line)
    require(len(lines)==limit,'committed log prefix incomplete')
    (out/(prefix+'SOURCE.json')).write_bytes(metadata_bytes)
    (out/(prefix+'LIVE_CHECKPOINT_SNAPSHOT.json')).write_bytes(checkpoint_bytes)
    (out/(prefix+'moves.jsonl')).write_text(''.join(lines))
    return meta,check,lines


def replay_segment(tri,meta,check,lines,label,start):
    require(state_sha(tri)==meta['start_tight_sha256'],'wrong '+label+' starting state')
    counts=Counter()
    for i,line in enumerate(lines):
        row=json.loads(line)
        require(row['index']==i and row['before_sha256']==state_sha(tri),'broken '+label+' log chain')
        guarded(tri,row['operation'])
        require(state_sha(tri)==row['after_sha256'] and tri.size()==row['pentachora'],'wrong '+label+' move result')
        counts[row['operation'][0]]+=1
        if (i+1)%100==0:
            print(json.dumps({'segment':label,'replayed':i+1,'pentachora':tri.size(),'seconds':time.monotonic()-start}),flush=True)
    require(len(lines)==check['moves'] and state_sha(tri)==check['current_tight_sha256'],'wrong '+label+' endpoint')
    return dict(counts)


def freeze(live,out):
    started=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    meta,check,lines=committed_snapshot(live,out,'')
    restore=meta['restored_by_checked_prefix']
    pm,pc,pl=committed_snapshot(restore['directory'],out,'prefix_')
    require(pc['moves']==restore['moves'],'restored prefix changed')
    original=Path(restore['original_file'])
    require(file_sha(original)==restore['original_sha256']==pm['source_sha256'],'original source checksum mismatch')
    require(file_sha(meta['source_file'])==meta['source_sha256'],'immediate source checksum mismatch')
    save_json(out/'STATUS.json',{'status':'SNAPSHOT_FROZEN_REPLAY_PENDING','moves':len(lines),'expected_pentachora':check['current_pentachora'],'snapshot_utc':datetime.now(timezone.utc).isoformat()})
    print(json.dumps({'frozen_phase_moves':len(lines),'expected_pentachora':check['current_pentachora']}),flush=True)
    tri=regina.Triangulation4(regina.open(str(original)))
    require(tri.isValid() and tri.isConnected() and tri.isOrientable() and not tri.isIdeal(),'original finite admission failed')
    print(json.dumps({'original_admitted':tri.size(),'seconds':time.monotonic()-started}),flush=True)
    prefix_counts=replay_segment(tri,pm,pc,pl,'source_prefix',started)
    immediate=regina.Triangulation4(regina.open(str(meta['source_file'])))
    require(tri.tightEncoding()==immediate.tightEncoding(),'prefix does not equal immediate source')
    counts=replay_segment(tri,meta,check,lines,'phase3',started)
    require(tri.isValid() and tri.isConnected() and tri.isOrientable() and not tri.isIdeal(),'endpoint finite admission failed')
    require(tri.countBoundaryComponents()==1 and tri.boundaryComponent(0).isReal(),'expected one real boundary')
    lower=out/'lower_reference.rga'
    require(regina.PacketOfTriangulation4(tri).save(str(lower)),'lower save failed')
    boundary,locations=boundary_data(tri)
    boundary_file=out/'boundary.rga'
    require(regina.PacketOfTriangulation3(boundary).save(str(boundary_file)),'boundary save failed')
    save_json(out/'boundary_locations.json',locations)
    boundary_sig=boundary.isoSig()
    (out/'boundary.isosig').write_text(boundary_sig+'\n')
    require(len(locations)==boundary.size()==tri.countBoundaryFacets(),'boundary face-location count mismatch')
    report={'status':'PASS_FROZEN_PREFIX_AND_BOUNDARY_EXPORT','live_directory':str(Path(live).resolve()),'original_source':str(original.resolve()),'original_source_sha256':file_sha(original),'prefix_moves':len(pl),'prefix_counts':prefix_counts,'phase_moves':len(lines),'phase_counts':counts,'phase_log_sha256':file_sha(out/'moves.jsonl'),'lower_pentachora':tri.size(),'lower_vertices':tri.countVertices(),'lower_tight_sha256':state_sha(tri),'lower_rga_sha256':file_sha(lower),'lower_file':str(lower.resolve()),'boundary_tetrahedra':boundary.size(),'boundary_vertices':boundary.countVertices(),'boundary_valid':boundary.isValid(),'boundary_closed':boundary.isClosed(),'boundary_connected':boundary.isConnected(),'boundary_orientable':boundary.isOrientable(),'boundary_isosig':boundary_sig,'boundary_rga_sha256':file_sha(boundary_file),'boundary_file':str(boundary_file.resolve()),'boundary_locations_sha256':file_sha(out/'boundary_locations.json'),'seconds':time.monotonic()-started,'regina':'7.4.1','native_4d_simplify_used':False,'cached_flags_manually_modified':False,'live_endpoint_copied':False,'closed_completion_asserted':False,'standardness_asserted':False,'scope':'Figure2 red over dashed green working input; finite lower body, awaiting actual upper handle gluing','trust_boundary':'Regina manifold admission and guarded moves; independent replay driver; direct boundary face maps'}
    save_json(out/'STATUS.json',report)
    frozen=dict(check);frozen.update({'status':'FROZEN_REGENERATED_BY_GUARDED_REPLAY','current_rga_sha256':file_sha(lower),'current_file':'lower_reference.rga'})
    save_json(out/'CHECKPOINT.json',frozen)
    # Also provide one directly replayable chain from the easy original
    # source, so consumers do not cold-load the difficult intermediate state.
    combined=out/'from_original';combined.mkdir()
    combined_rows=[]
    for index,line in enumerate(pl+lines):
        row=json.loads(line);row['index']=index
        combined_rows.append(json.dumps(row,separators=(',',':'))+'\n')
    (combined/'moves.jsonl').write_text(''.join(combined_rows))
    combined_meta=dict(pm)
    combined_meta.update({'schema':'motoo.guarded-finite-path.v1','independently_frozen':True,'segments':[len(pl),len(lines)]})
    save_json(combined/'SOURCE.json',combined_meta)
    require(regina.PacketOfTriangulation4(tri).save(str(combined/'current.rga')),'combined endpoint save failed')
    combined_counts=Counter(prefix_counts);combined_counts.update(counts)
    save_json(combined/'CHECKPOINT.json',{'status':'FROZEN_REGENERATED_BY_GUARDED_REPLAY','moves':len(combined_rows),'current_pentachora':tri.size(),'current_tight_sha256':state_sha(tri),'current_rga_sha256':file_sha(combined/'current.rga'),'counts':dict(combined_counts),'standard_endpoint':False})
    print(json.dumps({k:v for k,v in report.items() if k!='boundary_isosig'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('live',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();freeze(a.live,a.out)

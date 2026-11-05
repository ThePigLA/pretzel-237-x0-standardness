#!/usr/bin/env python3
"""Bounded isolated Regina microprobes; never modifies the input triangulation.

Native simplification endpoints are exploratory, not explicit move certificates.
Every potentially slow operation runs in a child, with JSON-lines phase timings
flushed before entry. A hard timeout leaves the exact source and receipts intact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def emit(event, **kw):
    print(json.dumps({"event": event, "monotonic": time.monotonic(), "cpu_seconds":time.process_time(), **kw}), flush=True)


def worker(args):
    import regina
    started = time.monotonic()
    tri = regina.Triangulation4(regina.open(str(args.source)))
    emit("loaded", size=tri.size(), seconds=time.monotonic()-started)
    t = time.monotonic()
    emit("enter_skeleton")
    vertices = tri.countVertices()
    emit("skeleton", seconds=time.monotonic()-t, size=tri.size(), vertices=vertices,
         valid=tri.isValid(), orientable=tri.isOrientable(),
         real_boundary_facets=tri.countBoundaryFacets())
    t = time.monotonic()
    emit("enter_operation", method=args.method)
    if args.method == "baseline":
        value = True
    elif args.method == "local-query":
        value = tri.simplifyToLocalMinimum(False)
    elif args.method == "local-minimum":
        value = tri.simplifyToLocalMinimum()
    elif args.method == "updown-small":
        value = tri.simplifyUpDown(2, 4, False)
    elif args.method == "collapse-scan":
        candidates = sorted((e.degree(),e.index()) for e in tri.edges()
                            if e.vertex(0)!=e.vertex(1))
        emit("collapse_candidates", count=len(candidates))
        value = 0
        for degree, index in candidates:
            s = time.monotonic()
            edge = tri.edge(index)
            endpoints = [edge.vertex(i).index() for i in range(2)]
            emit("enter_collapse_guard", edge=index, degree=degree, endpoints=endpoints)
            legal = tri.hasCollapseEdge(edge)
            emit("collapse_guard", edge=index, legal=legal, seconds=time.monotonic()-s)
            value += int(legal)
    elif args.method == "classify-links":
        value = []
        for vertex in tri.vertices():
            s = time.monotonic()
            emit("enter_vertex_type", vertex=vertex.index(), degree=vertex.degree())
            ideal = vertex.isIdeal()
            emit("vertex_type", vertex=vertex.index(), ideal=ideal, seconds=time.monotonic()-s)
            value.append(ideal)
    elif args.method == "snap-single":
        value=0
        (args.out/'SNAP_API.txt').write_text(regina.Triangulation4.snapEdge.__doc__+"\n\n"+regina.Triangulation4.hasSnapEdge.__doc__)
        def state_hash():
            return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()
        log=(args.out/'guarded_moves.jsonl').open('w')
        operation_index=0
        while tri.countVertices()>1:
            edge=next((e for e in tri.edges() if e.vertex(0)!=e.vertex(1) and tri.hasSnapEdge(e)),None)
            if edge is None:raise RuntimeError("no guarded snap candidate")
            s=time.monotonic();index=edge.index();before=tri.size()
            before_hash=state_hash()
            endpoints=[{"index":edge.vertex(v).index(),"ideal":edge.vertex(v).isIdeal(),"boundary":edge.vertex(v).isBoundary(),"valid":edge.vertex(v).isValid()} for v in (0,1)]
            count_before=tri.countVertices()
            guarded=tri.hasSnapEdge(edge)
            if not tri.snapEdge(edge):raise RuntimeError("guarded snap failed")
            value+=1
            after_hash=state_hash()
            row={"index":operation_index,"operation":["S",1,index],"before_sha256":before_hash,"after_sha256":after_hash,"before_pentachora":before,"pentachora":tri.size(),"guard_accepted":guarded,"endpoints_before":endpoints,"vertices_before":count_before}
            log.write(json.dumps(row)+"\n");log.flush();operation_index+=1
            regina.PacketOfTriangulation4(tri).save(str(args.out/'endpoint.rga'))
            emit("snap",edge=index,before=before,after=tri.size(),seconds=time.monotonic()-s)
            s=time.monotonic()
            count=tri.countVertices()
            if count!=count_before-1 or tri.size()!=before+4:raise RuntimeError("snap size or vertex count mismatch")
            emit("snap_rebuild",vertices=count,seconds=time.monotonic()-s)
        single_info={"vertices":tri.countVertices(),"ideal":tri.isIdeal(),"valid":tri.isValid(),"pentachora":tri.size(),"boundary_components":tri.countBoundaryComponents(),"real_boundary_facets":tri.countBoundaryFacets(),"tight_sha256":state_hash(),"snaps":value}
        regina.PacketOfTriangulation4(tri).save(str(args.out/'single_ideal.rga'))
        (args.out/'SINGLE_IDEAL.json').write_text(json.dumps(single_info,indent=2)+"\n")
        emit("all_ideal",**single_info)
        for step in range(30):
            s=time.monotonic();face=None
            for dimension in (2,1):
                face=next((f for f in tri.faces(dimension) if f.degree()==2 and tri.has20(f)),None)
                if face is not None:break
            if face is None:break
            index=face.index();before=tri.size()
            before_hash=state_hash()
            if not tri.move20(face):raise RuntimeError("guarded 2-0 failed")
            row={"index":operation_index,"operation":["Z",dimension,index],"before_sha256":before_hash,"after_sha256":state_hash(),"before_pentachora":before,"pentachora":tri.size(),"guard_accepted":True}
            log.write(json.dumps(row)+"\n");log.flush();operation_index+=1
            emit("post_snap_move",index=step,dimension=dimension,face=index,before=before,after=tri.size(),seconds=time.monotonic()-s)
            regina.PacketOfTriangulation4(tri).save(str(args.out/'endpoint.rga'))
            s=time.monotonic();count=tri.countVertices()
            emit("post_snap_rebuild",index=step,vertices=count,seconds=time.monotonic()-s)
        log.close()
    elif args.method == "ideal-first":
        from ideal_first_relabelling import ideal_first, encode_isomorphism
        original=regina.Triangulation4(tri)
        emit("ideal_vertex_details",vertices=[{"index":v.index(),"degree":v.degree(),"ideal":v.isIdeal()} for v in tri.vertices()])
        s=time.monotonic()
        iso,attempts=ideal_first(tri)
        (args.out/'relabelling.json').write_text(json.dumps(encode_isomorphism(iso))+"\n")
        if iso(original).tightEncoding()!=tri.tightEncoding():
            raise RuntimeError("relabelling composite failed exact-state verification")
        emit("relabeled",attempts=attempts,seconds=time.monotonic()-s)
        s=time.monotonic()
        count=tri.countVertices()
        emit("ideal_first_skeleton",seconds=time.monotonic()-s,vertices=count,first_ideal=tri.simplex(0).vertex(0).isIdeal())
        value=0
        for step in range(30):
            s=time.monotonic()
            face=None
            for dimension in (2,1):
                face=next((f for f in tri.faces(dimension)
                           if f.degree()==2 and
                           all(e.simplex().index()!=0 for e in f.embeddings()) and
                           tri.has20(f)),None)
                if face is not None:break
            if face is None:break
            index=face.index();before=tri.size()
            if not tri.move20(face):raise RuntimeError("guarded 2-0 failed")
            emit("ideal_first_move",index=step,dimension=dimension,face=index,before=before,after=tri.size(),seconds=time.monotonic()-s)
            regina.PacketOfTriangulation4(tri).save(str(args.out/'endpoint.rga'))
            s=time.monotonic()
            count=tri.countVertices()
            emit("ideal_first_rebuild",index=step,seconds=time.monotonic()-s,vertices=count,first_ideal=tri.simplex(0).vertex(0).isIdeal())
            value+=1
    elif args.method == "shell-greedy":
        value = 0
        while True:
            candidate = next((p.index() for p in tri.simplices()
                              if p.hasBoundary() and tri.hasShellBoundary(p)), None)
            if candidate is None:
                break
            before = tri.size()
            if not tri.shellBoundary(tri.simplex(candidate)):
                raise RuntimeError("guarded shell failed")
            value += 1
            emit("shell", before=before, after=tri.size())
            if value % 10 == 0:
                regina.PacketOfTriangulation4(tri).save(str(args.out / "endpoint.rga"))
    else:
        raise ValueError(args.method)
    emit("operation", method=args.method, value=value, seconds=time.monotonic()-t,
         size=tri.size())
    if not regina.PacketOfTriangulation4(tri).save(str(args.out / "endpoint.rga")):
        raise RuntimeError("cannot save endpoint")
    (args.out / "endpoint.tight").write_text(tri.tightEncoding()+"\n")
    emit("complete", total_seconds=time.monotonic()-started, size=tri.size())


def main():
    p=argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--method", choices=["baseline","local-query","local-minimum","updown-small","collapse-scan","classify-links","shell-greedy","ideal-first","snap-single"], required=True)
    p.add_argument("--seconds", type=float, default=45)
    p.add_argument("--worker", action="store_true")
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.worker:
        worker(a)
        return
    command=[sys.executable,str(Path(__file__).resolve()),str(a.source.resolve()),
             "--out",str(a.out.resolve()),"--method",a.method,"--worker"]
    start=time.monotonic()
    timed_out=False
    with (a.out / "probe.jsonl").open("w") as log:
        child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            rc=child.wait(timeout=a.seconds)
        except (subprocess.TimeoutExpired,KeyboardInterrupt):
            timed_out=True
            os.killpg(child.pid,signal.SIGTERM)
            try: rc=child.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGKILL)
                rc=child.wait()
    receipt={"source":str(a.source.resolve()),
             "source_sha256":hashlib.sha256(a.source.read_bytes()).hexdigest(),
             "method":a.method,"timeout_seconds":a.seconds,"timed_out":timed_out,
             "elapsed_seconds":time.monotonic()-start,"returncode":rc,
             "endpoint_exists":(a.out/"endpoint.rga").exists(),
             "scope":"performance diagnostic on a copy; native endpoints are not move certificates"}
    (a.out/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt),flush=True)
    print((a.out/"probe.jsonl").read_text(),end="")


if __name__=="__main__": main()

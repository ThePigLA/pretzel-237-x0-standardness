#!/usr/bin/env python3
"""Replay a diagnostic S/Z/P log; independently reconstruct each snap gadget.

This audits the explicit four-pentachoron gluing and endpoint hypotheses. It
does not independently prove the topological theorem underlying edge snaps;
that theorem and Regina's standard Pachner/2-0 move semantics remain premises.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import regina


def require(condition,message):
    if not condition:raise ValueError(message)


def digest(tri):return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()


def explicit_snap_gadget(tri,edge_index):
    """Return a copied triangulation with the four cells glued by public APIs."""
    answer=regina.Triangulation4(tri)
    edge=answer.edge(edge_index)
    require(edge.vertex(0)!=edge.vertex(1),"snap endpoints are identical")
    require(any(not edge.vertex(v).isBoundary() for v in (0,1)),"both endpoints are boundary/ideal")
    require(all(edge.vertex(v).isValid() for v in (0,1)),"invalid snap endpoint")
    embedding=edge.embedding(0)
    opened=embedding.simplex()
    perm=embedding.vertices()
    adjacent=opened.adjacentSimplex(perm[2])
    require(adjacent is not None,"selected gadget tetrahedron has real boundary")
    glue=opened.adjacentGluing(perm[2])
    require(not opened.isFacetLocked(perm[2]),"audit needs an unlocked tetrahedron")
    cells=[answer.newSimplex() for _ in range(4)]
    swap34=regina.Perm5(3,4)
    cells[0].join(0,cells[1],swap34)
    cells[0].join(2,cells[1],regina.Perm5(0,2,4,1,3))
    cells[0].join(3,cells[2],swap34)
    cells[0].join(4,cells[2],swap34)
    cells[1].join(1,cells[2],regina.Perm5(1,2))
    cells[1].join(2,cells[3],swap34)
    cells[1].join(3,cells[3],swap34)
    cells[2].join(0,cells[3],swap34)
    cells[2].join(1,cells[3],swap34)
    opened.unjoin(perm[2])
    cells[0].join(1,opened,perm*regina.Perm5(3,2,0,1,4))
    cells[3].join(3,adjacent,glue*perm*regina.Perm5(3,1,0,2,4))
    return answer


def audit(source,log,out):
    started=time.monotonic()
    tri=regina.Triangulation4(regina.open(str(source)))
    require(tri.isValid(),"invalid starting triangulation")
    source_info={"pentachora":tri.size(),"vertices":tri.countVertices(),"ideal":tri.isIdeal(),"real_boundary_facets":tri.countBoundaryFacets(),"boundary_components":tri.countBoundaryComponents()}
    snap_checks=[];moves=0
    for line in Path(log).read_text().splitlines():
        row=json.loads(line);operation=row['operation']
        if operation[0]=='S' and len(operation)==2:
            kind,dimension,index='S',1,operation[1]
        else:
            kind,dimension,index=operation
        require(row['index']==moves,"nonsequential log")
        require(digest(tri)==row['before_sha256'],"before state hash mismatch")
        before=tri.size()
        if kind=='S':
            require(dimension==1,"snap must be an edge operation")
            count=tri.countVertices();edge=tri.edge(index)
            endpoints=[{"vertex":edge.vertex(v).index(),"boundary":edge.vertex(v).isBoundary(),"ideal":edge.vertex(v).isIdeal(),"valid":edge.vertex(v).isValid()} for v in (0,1)]
            expected=explicit_snap_gadget(tri,index)
            require(tri.hasSnapEdge(edge),"Regina snap guard rejected")
            require(tri.snapEdge(edge),"Regina checked snap rejected")
            require(tri.tightEncoding()==expected.tightEncoding(),"explicit snap gadget differs from Regina output")
            require(tri.size()==before+4 and tri.countVertices()==count-1,"wrong snap cell/vertex count")
            snap_checks.append({"move_index":moves,"edge":index,"endpoints":endpoints,"before_vertices":count,"after_vertices":tri.countVertices(),"explicit_gadget_matches":True,"after_sha256":digest(tri)})
            print(json.dumps({"audited_snap":len(snap_checks),"vertices":tri.countVertices(),"pentachora":tri.size(),"seconds":time.monotonic()-started}),flush=True)
        elif kind=='Z':
            face=tri.face(dimension,index)
            require(tri.has20(face) and tri.move20(face),"2-0 guard rejected")
            require(tri.size()==before-2,"wrong 2-0 size")
        elif kind=='P':
            face=tri.simplex(index) if dimension==4 else tri.face(dimension,index)
            require(tri.hasPachner(face) and tri.pachner(face),"Pachner guard rejected")
            require(tri.size()==before+2*dimension-4,"wrong Pachner size")
        else:raise ValueError("unsupported operation")
        require(digest(tri)==row['after_sha256'],"after state hash mismatch")
        moves+=1
    result={"status":"PASS_GUARDS_HASH_CHAIN_AND_EXPLICIT_SNAP_GADGETS","source":str(Path(source).resolve()),"source_sha256":hashlib.sha256(Path(source).read_bytes()).hexdigest(),"move_log_sha256":hashlib.sha256(Path(log).read_bytes()).hexdigest(),"source_properties":source_info,"moves":moves,"snap_checks":snap_checks,"final_pentachora":tri.size(),"final_vertices":tri.countVertices(),"final_valid":tri.isValid(),"final_ideal":tri.isIdeal(),"final_real_boundary_facets":tri.countBoundaryFacets(),"final_boundary_components":tri.countBoundaryComponents(),"final_tight_sha256":digest(tri),"seconds":time.monotonic()-started,"standardness_asserted":False,"trust_boundary":"snap topology theorem and Regina guarded Pachner/2-0 move semantics; snap gluing independently reconstructed using public simplex joins"}
    Path(out).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('log',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();audit(args.source,args.log,args.out)

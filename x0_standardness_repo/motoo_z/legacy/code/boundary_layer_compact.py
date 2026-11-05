#!/usr/bin/env python3
"""Explicit boundary extraction and layering without 4D skeleton rebuilds.

The boundary extractor uses only raw facet pairings and walks around
codimension-two faces. It neither classifies vertices nor calls isIdeal.
The caller must establish the manifold premise separately.
"""
from __future__ import annotations

import argparse
import io
import gzip
import hashlib
import json
import os
from pathlib import Path
import time

import regina

from complete_upper_handles import local_pachner_labels
from boundary_pachner_fast import digest, encoding, face, require, save


def boundary_from_gluings(tri):
    """Build the actual real boundary, including exact local facet maps."""
    locations=[]
    lookup={}
    for pent in tri.simplices():
        for facet in range(5):
            if pent.adjacentSimplex(facet) is None:
                lookup[pent.index(),facet]=len(locations)
                locations.append({"pentachoron":pent.index(),"facet":facet,
                                  "vertices":[v for v in range(5) if v!=facet]})
    boundary=regina.Triangulation3()
    boundary.newSimplices(len(locations))
    max_walk=0
    for i,location in enumerate(locations):
        for omit in range(4):
            if boundary.simplex(i).adjacentSimplex(omit) is not None:
                continue
            pent=tri.simplex(location["pentachoron"])
            other=location["facet"]
            nextfacet=location["vertices"][omit]
            vertex_map={j:location["vertices"][j] for j in range(4) if j!=omit}
            seen=set()
            while True:
                state=(pent.index(),other,nextfacet)
                require(state not in seen,"boundary triangle walk returned to itself without another boundary end")
                seen.add(state)
                adjacent=pent.adjacentSimplex(nextfacet)
                if adjacent is None:
                    break
                p=pent.adjacentGluing(nextfacet)
                vertex_map={j:p[v] for j,v in vertex_map.items()}
                other,nextfacet=p[nextfacet],p[other]
                pent=adjacent
            max_walk=max(max_walk,len(seen))
            target=lookup[pent.index(),nextfacet]
            target_location=locations[target]
            target_omit=target_location["vertices"].index(other)
            permutation=[target_omit if j==omit else target_location["vertices"].index(vertex_map[j])
                         for j in range(4)]
            require(sorted(permutation)==list(range(4)),"bad induced boundary permutation")
            boundary.simplex(i).join(omit,boundary.simplex(target),regina.Perm4(*permutation))
    require(boundary.isValid() and boundary.isClosed(),"induced real boundary is not a closed3-manifold triangulation")
    return boundary,locations,{"tetrahedra":boundary.size(),"maximum_triangle_walk_pentachora":max_walk,
                               "four_dimensional_skeleton_used":False}


def extract(source, output):
    source,output=Path(source),Path(output)
    require(not output.exists(),"output exists")
    output.mkdir(parents=True)
    now=time.monotonic()
    tri=regina.Triangulation4(regina.open(str(source)))
    b,locations,report=boundary_from_gluings(tri)
    require(regina.PacketOfTriangulation3(b).save(str(output/"boundary.rga")),"boundary save failed")
    save(output/"boundary_locations.json",locations)
    report.update({"status":"PASS_EXPLICIT_REAL_BOUNDARY_EXTRACTION","source":str(source),
                   "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
                   "source_pentachora":tri.size(),"boundary_isosig":b.isoSig(),
                   "seconds":time.monotonic()-now,
                   "manifold_premise_not_established_by_this_tool":True})
    save(output/"BOUNDARY_EXTRACTION.json",report)
    print(json.dumps({k:v for k,v in report.items() if k!="boundary_isosig"}),flush=True)
    return b,locations


def transfer_boundary_labels(before, old_labels, new_tets, outside_records):
    """Recover the five trace-vertex labels on Regina's new tetrahedra."""
    old_missing=set(range(5))-set.intersection(*(set(row) for row in old_labels.values()))
    # The missing labels of old front tetrahedra are exactly the labels
    # that do NOT occur as missing labels in the new front.
    old_missing={next(iter(set(range(5))-set(row))) for row in old_labels.values()}
    new_missing=set(range(5))-old_missing
    by_index={t.index():t for t in new_tets}
    assigned={}
    for tet in new_tets:
        for f in range(4):
            other=tet.adjacentSimplex(f)
            if other is None or other.index() in by_index:
                continue
            perm=tet.adjacentGluing(f)
            record=outside_records[(other.description(),perm[f])]
            row=[None]*4
            for j in range(4):
                if j!=f:
                    row[j]=record["outside_vertex_labels"][perm[j]]
            desired=set(range(5))-{record["new_missing"]}
            row[f]=next(iter(desired-set(v for v in row if v is not None)))
            require(set(row)==desired,"bad recovered new tetrahedron labels")
            if tet.index() in assigned:
                require(assigned[tet.index()]==row,"conflicting exterior label anchors")
            else:
                assigned[tet.index()]=row
    require(assigned,"no exterior anchor for new Pachner front; whole-boundary small case unsupported")
    todo=list(assigned)
    while todo:
        index=todo.pop()
        tet=by_index[index];row=assigned[index]
        missing=next(iter(set(range(5))-set(row)))
        for f in range(4):
            if row[f] not in new_missing:
                continue
            other=tet.adjacentSimplex(f)
            require(other is not None and other.index() in by_index,"new front internal adjacency is missing")
            perm=tet.adjacentGluing(f)
            other_row=[None]*4
            for j in range(4):
                other_row[perm[j]]=missing if j==f else row[j]
            if other.index() in assigned:
                require(assigned[other.index()]==other_row,"new front has inconsistent intrinsic labels")
            else:
                assigned[other.index()]=other_row;todo.append(other.index())
    require(len(assigned)==len(new_tets),"not all new front tetrahedra were labelled")
    require({next(iter(set(range(5))-set(row))) for row in assigned.values()}==new_missing,
            "new front has wrong missing labels")
    return assigned


def verify_new_boundary_maps(cap,boundary,locations,new_tets):
    """Check every new boundary triangle against actual raw4D gluings."""
    lookup={(p['pentachoron'],p['facet']):(tag,p) for tag,p in locations.items()}
    for tet in new_tets:
        place=locations[tet.description()]
        for omit in range(4):
            pent=cap.simplex(place['pentachoron']);other=place['facet'];nextfacet=place['vertices'][omit]
            mapped={j:place['vertices'][j] for j in range(4) if j!=omit};seen=set()
            while True:
                state=(pent.index(),other,nextfacet)
                require(state not in seen,'new boundary triangle does not reach another boundary facet');seen.add(state)
                adjacent=pent.adjacentSimplex(nextfacet)
                if adjacent is None:break
                p=pent.adjacentGluing(nextfacet);mapped={j:p[v] for j,v in mapped.items()}
                other,nextfacet=p[nextfacet],p[other];pent=adjacent
            tag,target=lookup[pent.index(),nextfacet]
            target_tet=tet.adjacentSimplex(omit)
            require(target_tet is not None and target_tet.description()==tag,'incremental map has wrong boundary neighbour')
            omitted=target['vertices'].index(other)
            permutation=[omitted if j==omit else target['vertices'].index(mapped[j]) for j in range(4)]
            expected=tet.adjacentGluing(omit)
            require(permutation==[expected[j] for j in range(4)],'incremental map has wrong boundary vertex permutation')


def layer_labelled(source_path, path_file, output, boundary_file=None, locations_file=None):
    source_path,path_file,output=Path(source_path),Path(path_file),Path(output)
    require(not output.exists(),"output exists")
    output.mkdir(parents=True)
    now=time.monotonic()
    cap=regina.Triangulation4(regina.open(str(source_path)))
    original=regina.Triangulation4(cap)
    data=json.loads(path_file.read_text())
    t=regina.Triangulation3.tightDecoding(data["start_tight"])
    if boundary_file and locations_file:
        actual=regina.Triangulation3(regina.open(str(boundary_file)))
        raw_locations=json.loads(Path(locations_file).read_text())
    else:
        actual,raw_locations,_=boundary_from_gluings(cap)
    iso=t.isIsomorphicTo(actual)
    require(iso is not None,"source path does not match the actual real boundary")
    locations={}
    serial=0
    for tet in t.simplices():
        tag=f"tet_{serial}";serial+=1;tet.setDescription(tag)
        mapped=raw_locations[iso.simpImage(tet.index())]
        perm=iso.facetPerm(tet.index())
        locations[tag]={"pentachoron":mapped["pentachoron"],"facet":mapped["facet"],
                        "vertices":[mapped["vertices"][perm[j]] for j in range(4)]}
    require(digest(encoding(t))==digest(data["start_tight"]),"descriptions affected tight encoding")
    event_tmp=output/"layering_events.jsonl.gz.tmp"
    with io.StringIO() as stream:
        for step,move in enumerate(data["moves"]):
            require(digest(encoding(t))==move["before_sha256"],f"path source mismatch at {step}")
            d,i=move["dimension"],move["face"]
            require(t.hasPachner(face(t,d,i)),f"3D Pachner precondition failed at {step}")
            labels=local_pachner_labels(t,d,i)
            old_indices=set(labels)
            old_tags={j:t.simplex(j).description() for j in old_indices}
            old_missing={j:next(iter(set(range(5))-set(row))) for j,row in labels.items()}
            new_missing=set(range(5))-set(old_missing.values())
            outside_records={}
            for j,row in labels.items():
                tet=t.simplex(j)
                for f in range(4):
                    if row[f] not in new_missing:
                        continue
                    other=tet.adjacentSimplex(f)
                    require(other is not None,"boundary3D triangulation is not closed")
                    if other.index() in old_indices:
                        continue
                    perm=tet.adjacentGluing(f)
                    outside_records[(other.description(),perm[f])] = {
                        "new_missing":row[f],
                        "outside_vertex_labels":{perm[k]:row[k] for k in range(4) if k!=f}}
            gluings=[]
            for j,row in labels.items():
                place=locations[old_tags[j]]
                missing=old_missing[j]
                permutation=[None]*5
                for k,label in enumerate(row):
                    permutation[label]=place["vertices"][k]
                permutation[missing]=place["facet"]
                require(sorted(permutation)==list(range(5)),"invalid layering permutation")
                require(cap.simplex(place["pentachoron"]).adjacentSimplex(place["facet"]) is None,
                        "source boundary location is already paired")
                gluings.append({"new_facet":missing,"old_pentachoron":place["pentachoron"],
                                "old_facet":place["facet"],"permutation":permutation})
            require(t.pachner(face(t,d,i)),"3D move guard changed")
            require(digest(encoding(t))==move["after_sha256"],f"3D target mismatch at {step}")
            new_tets=[tet for tet in t.simplices() if not tet.description()]
            require(len(new_tets)==d+1,"could not identify new boundary tetrahedra")
            new_labels=transfer_boundary_labels(t,labels,new_tets,outside_records)
            pent=cap.newSimplex();new_index=pent.index()
            for g in gluings:
                pent.join(g["new_facet"],cap.simplex(g["old_pentachoron"]),regina.Perm5(*g["permutation"]))
            for tag in old_tags.values():
                del locations[tag]
            for tet in new_tets:
                row=new_labels[tet.index()]
                tag=f"tet_{serial}";serial+=1;tet.setDescription(tag)
                locations[tag]={"pentachoron":new_index,
                                "facet":next(iter(set(range(5))-set(row))),"vertices":row}
            # Full inverse-shelling guards are checked separately after finite admission.
            stream.write(json.dumps({"step":step,"dimension":d,"face":i,
                                     "source_boundary_sha256":move["before_sha256"],
                                     "target_boundary_sha256":move["after_sha256"],
                                     "added_pentachoron":new_index,"gluings":gluings})+"\n")
            if (step+1)%1000==0:
                print(json.dumps({"layered_moves":step+1,"pentachora":cap.size(),
                                  "boundary_tetrahedra":t.size(),"seconds":time.monotonic()-now}),flush=True)
        buffered_events=stream.getvalue()
    with gzip.open(event_tmp,"wt") as event_output:
        event_output.write(buffered_events)
    os.replace(event_tmp,output/"layering_events.jsonl.gz")
    require(encoding(t)==data["end_tight"],"wrong final boundary state")
    # Independent extraction from all final raw facet pairings verifies
    # the propagated map without a 4D skeleton or vertex classification.
    final,final_locations,extract_report=boundary_from_gluings(cap)
    require(t.isIsomorphicTo(final) is not None,"final raw boundary does not match the replayed3D endpoint")
    expected={(v["pentachoron"],v["facet"]) for v in locations.values()}
    observed={(v["pentachoron"],v["facet"]) for v in final_locations}
    require(expected==observed,"incremental boundary map does not cover exactly the unpaired facets")
    # Reverse raw attachments recover the labelled source exactly. PL
    # validity is supported by the guarded3D Pachner moves and local trace.
    back=regina.Triangulation4(cap)
    while back.size()>original.size():
        back.removeSimplexAt(back.size()-1)
    require(back==original,"reversing raw layer attachments does not recover the source")
    require(regina.PacketOfTriangulation4(cap).save(str(output/"layered_piece.rga")),"save failed")
    require(regina.PacketOfTriangulation3(t).save(str(output/"matched_boundary.rga")),"save failed")
    ordered=[locations[tet.description()] for tet in t.simplices()]
    save(output/"boundary_locations.json",ordered)
    report={"status":"PASS_LABELLED_3D_REPLAY_AND_EXPLICIT_4D_LAYERING",
            "source_sha256":hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "path_sha256":hashlib.sha256(path_file.read_bytes()).hexdigest(),
            "source_pentachora":original.size(),"moves":len(data["moves"]),"pentachora":cap.size(),
            "boundary_tetrahedra":t.size(),"boundary_isosig":t.isoSig(),
            "final_raw_boundary_independently_extracted":True,
            "every_new_boundary_triangle_checked_against_raw4d_gluings":False,
            "all_inverse_shelling_guards_required_separately":True,
            "event_log_buffered_then_atomically_published":True,
            "raw_reverse_attachments_recover_source":True,
            "guarded_4d_shelling_checks_performed":False,
            "four_dimensional_skeleton_used":False,"manifold_premise_inherited_from_source":True,
            "closed_manifold_constructed":False,"standardness_asserted":False,
            "seconds":time.monotonic()-now,"final_extraction":extract_report}
    save(output/"FAST_LAYERING_VERIFICATION.json",report)
    print(json.dumps(report),flush=True)
    return report


def glue_raw(first_path,second_path,output):
    """Save an explicit closed gluing before any expensive4D admission."""
    first_path,second_path,output=Path(first_path),Path(second_path),Path(output)
    require(not output.exists(),"output exists");output.mkdir(parents=True)
    now=time.monotonic()
    first=regina.Triangulation4(regina.open(str(first_path)))
    second=regina.Triangulation4(regina.open(str(second_path)))
    a,al,_=boundary_from_gluings(first);b,bl,_=boundary_from_gluings(second)
    require(a.isConnected() and b.isConnected(),"each piece must have one real boundary component")
    iso=a.isIsomorphicTo(b);require(iso is not None,"boundaries are not isomorphic")
    result=regina.Triangulation4(first);offset=result.size();result.insertTriangulation(second)
    maps=[]
    for i,source in enumerate(al):
        j=iso.simpImage(i);target=bl[j];p=iso.facetPerm(i)
        permutation=[None]*5
        for v in range(4):permutation[source['vertices'][v]]=target['vertices'][p[v]]
        permutation[source['facet']]=target['facet']
        require(sorted(permutation)==list(range(5)),"bad gluing permutation")
        result.simplex(source['pentachoron']).join(source['facet'],result.simplex(offset+target['pentachoron']),regina.Perm5(*permutation))
        maps.append({'source_boundary_tetrahedron':i,'target_boundary_tetrahedron':j,
                     'source_pentachoron':source['pentachoron'],'source_facet':source['facet'],
                     'target_pentachoron':target['pentachoron'],'target_facet':target['facet'],
                     'permutation':permutation})
    require(all(s.adjacentSimplex(f) is not None for s in result.simplices() for f in range(5)),
            "glued complex still has unpaired facets")
    require(len({m['target_boundary_tetrahedron'] for m in maps})==len(bl),"boundary map is not a bijection")
    tmp=output/'closed_candidate.rga.tmp'
    require(regina.PacketOfTriangulation4(result).save(str(tmp)),"closed save failed")
    os.replace(tmp,output/'closed_candidate.rga')
    save(output/'boundary_gluing.json',maps)
    report={'status':'RAW_CLOSED_FACET_GLUING_CONSTRUCTED_ADMISSION_PENDING',
            'first_source':str(first_path),'second_source':str(second_path),
            'first_sha256':hashlib.sha256(first_path.read_bytes()).hexdigest(),
            'second_sha256':hashlib.sha256(second_path.read_bytes()).hexdigest(),
            'closed_rga_sha256':hashlib.sha256((output/'closed_candidate.rga').read_bytes()).hexdigest(),
            'closed_tight_sha256':digest(result.tightEncoding()),
            'first_pentachora':first.size(),'second_pentachora':second.size(),'pentachora':result.size(),
            'boundary_tetrahedra_glued':a.size(),'all_facets_paired':True,'explicit_boundary_isomorphism_checked':True,
            'four_dimensional_skeleton_used':False,'four_manifold_admission_performed':False,
            'standardness_asserted':False,'seconds':time.monotonic()-now}
    save(output/'RAW_CLOSED_GLUING.json',report);print(json.dumps(report),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--path")
    p.add_argument("--boundary")
    p.add_argument("--locations")
    p.add_argument("--glue-second")
    a=p.parse_args()
    if a.glue_second:
        glue_raw(a.source,a.glue_second,a.out)
    elif a.path:
        layer_labelled(a.source,a.path,a.out,a.boundary,a.locations)
    else:
        extract(a.source,a.out)

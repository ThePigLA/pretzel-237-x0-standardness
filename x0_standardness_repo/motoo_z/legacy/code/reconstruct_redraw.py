#!/usr/bin/env python3
"""Reconstruct candidate closed coloured curves with an explicit editing audit.

Uses principal ink paths and global endpoint matching. The output is an actual
closed polygonal trace, but its relationship to the drawing must be audited.
No over/under information, PD admission, or isotopy is asserted here.
"""
import argparse
import heapq
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

try:
    import networkx as nx
except ImportError:
    sys.path.append(str(Path(__file__).resolve().parents[1] / 'venv/lib/python3.12/site-packages'))
    import networkx as nx

from trace_redraw_centerlines import thin, make_graph


LABEL_ROIS = {
    1: {'yellow': [(383, 30, 411, 64), (794, 724, 810, 760), (741, 783, 754, 795)],
        'cyan': [(774, 162, 800, 194)], 'blue': [(696, 682, 704, 693)]},
    2: {'yellow': [(846, 0, 891, 53)], 'cyan': [(840, 196, 879, 241)],
        'purple': [(257, 267, 301, 319)], 'blue': [(1060, 116, 1122, 158)]},
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def connected(graph):
    unseen = set(graph); result = []
    while unseen:
        root = min(unseen); unseen.remove(root); todo = [root]; part = {root}
        for p in todo:
            for q in graph[p]:
                if q in unseen:
                    unseen.remove(q); part.add(q); todo.append(q)
        result.append(part)
    return result


def distances(graph, start):
    best = {start: 0}; previous = {}; queue = [(0, start)]
    while queue:
        dist, p = heapq.heappop(queue)
        if dist != best[p]: continue
        for q in graph[p]:
            nd = dist + math.dist(p, q)
            if nd < best.get(q, math.inf):
                best[q] = nd; previous[q] = p; heapq.heappush(queue, (nd, q))
    return best, previous


def longest_path(graph, part):
    leaves = [p for p in part if len(graph[p]) == 1]
    first = min(leaves or list(part))
    d, _ = distances(graph, first)
    first = max(d, key=d.get)
    d, previous = distances(graph, first)
    last = max(d, key=d.get)
    result = [last]
    while result[-1] != first: result.append(previous[result[-1]])
    result.reverse()
    return result


def simplify(points, tolerance=0.65):
    if len(points) <= 2: return points
    # Iterative Ramer--Douglas--Peucker, retaining first/last source points.
    keep = {0, len(points)-1}; pending = [(0, len(points)-1)]
    p = np.asarray(points)
    while pending:
        i, j = pending.pop()
        if j-i < 2: continue
        delta = p[j]-p[i]; length = np.linalg.norm(delta)
        if length:
            d = np.abs(delta[0]*(p[i+1:j,1]-p[i,1])-delta[1]*(p[i+1:j,0]-p[i,0]))/length
        else: d = np.linalg.norm(p[i+1:j]-p[i],axis=1)
        k = int(np.argmax(d))+i+1
        if d[k-i-1] > tolerance:
            keep.add(k); pending.extend([(i,k),(k,j)])
    return [points[i] for i in sorted(keep)]


def outward(points):
    p = np.asarray(points)
    cum = np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))
    k = min(len(p)-1, int(np.searchsorted(cum, 5.0))+1)
    # Average a short terminal part to reduce single-pixel hook sensitivity.
    v = p[0]-p[k]
    length = np.linalg.norm(v)
    return (v/length).tolist() if length else [0.0,0.0]


def gap_cost(a, b):
    pa, pb = np.asarray(a['position']),np.asarray(b['position'])
    delta=pb-pa; d=float(np.linalg.norm(delta))
    if not d:return 0.0
    u=delta/d
    ca=float(np.dot(a['outward'],u));cb=-float(np.dot(b['outward'],u))
    return d + 10*(2-ca-cb) + (max(0,d-20)**2)


def get_cycles(matching, n):
    unseen=set(range(n));cycles=[]
    while unseen:
        start=2*min(unseen);at=start;cycle=[]
        while True:
            fragment=at//2;unseen.discard(fragment);cycle.append(at)
            at=matching[at^1]
            if at==start:break
            if len(cycle)>n:raise ValueError('Invalid alternating cycle')
        cycles.append(cycle)
    return cycles


def reviewed_blue_cleanup(report, output, source, scale):
    """Trim three independently visually reviewed endpoint pen hooks in F1."""
    from assemble_redraw_pd import intersections
    expected='cfc1eb4f163e7c504089dfb96bceaec1ec5ed3074ab97c466db78a751c2eede2'
    inventory=json.loads((source/'GEOMETRY_INVENTORY.json').read_text())
    if inventory['source_sha256']!=expected:raise ValueError('Reviewed hook edits belong to another PDF')
    regions=[(842,185,851,194),(1127,193,1135,202),(1120,551,1130,560)]
    original=report['points_pdf'];points=list(original);kinds=list(report['edge_kinds']);edits=[]
    for region in regions:
        hits,_=intersections({'blue':{'points_pdf':points}})
        matches=[h for h in hits if region[0]<h['position_pdf'][0]<region[2] and region[1]<h['position_pdf'][1]<region[3]]
        if len(matches)!=1:raise ValueError('Reviewed pen-hook region no longer has one trace overlap')
        h=matches[0];i,j=sorted(v['segment'] for v in h['visits']);middle=points[i+1:j+1]
        span=[max(p[k] for p in middle)-min(p[k] for p in middle) for k in (0,1)]
        if max(span)>4:raise ValueError('Refusing to trim a nonlocal loop')
        edits.append({'kind':'reviewed_endpoint_retrace_trim','review_region_pdf':list(region),
                      'intersection_pdf':h['position_pdf'],'removed_points':middle,
                      'before_segment_indices':[i,j],'removed_span_pdf':span,
                      'reason':'Small pen/endpoint retrace reviewed in original composite; not an intended crossing'})
        points=points[:i+1]+[h['position_pdf']]+points[j+1:]
        kinds=kinds[:i]+['trimmed_ink_endpoint','trimmed_gap_endpoint']+kinds[j+1:]
    report=dict(report,points_pdf=points,edge_kinds=kinds,
                reviewed_local_edits=edits,
                original_candidate_points_sha256=hashlib.sha256(json.dumps(original,separators=(',',':')).encode()).hexdigest(),
                fragment_ranges_refer_to_pretrim_candidate=True,
                status='CLOSED_TRACE_WITH_RECORDED_VISUAL_CORRECTIONS_NATIVE_CHECKS_PENDING')
    (output/'blue_closed_reviewed.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    canvas=Image.open(source/'original_page.png').convert('RGB');draw=ImageDraw.Draw(canvas)
    for a,b in zip(points,points[1:]):draw.line(tuple(v*scale for v in a+b),fill=(0,0,0),width=1)
    canvas.save(output/'blue_reviewed_overlay.png')
    return report


def reconstruct(source, output, colour, figure, scale, threshold, minimum):
    started=time.monotonic();image_path=source/f'layer_{colour}.png'
    image=Image.open(image_path).convert('RGB');rgb=np.asarray(image)
    mask=np.max(255-rgb,axis=-1)>=threshold
    exclusion_records=[]
    labels, _ = ndimage.label(mask, structure=np.ones((3,3)))
    objects=ndimage.find_objects(labels)
    for roi in LABEL_ROIS.get(figure,{}).get(colour,[]):
        x0,y0,x1,y1=[int(v*scale) for v in roi]
        selected=[];excluded=0
        # Never crop a curve that merely passes through an annotation region.
        for identity, box in enumerate(objects,1):
            if box and box[0].start>=y0 and box[0].stop<=y1 and box[1].start>=x0 and box[1].stop<=x1:
                local=labels[box]==identity;excluded+=int(local.sum())
                mask[box][local]=False;selected.append(identity)
        exclusion_records.append({'reason':'visually_identified_annotation_or_stray_region',
                                  'whole_components_only':True,'bbox_pdf':roi,
                                  'removed_ink_component_ids':selected,'ink_pixels_removed':excluded})
    skel,iterations=thin(mask);graph=make_graph(skel);parts=connected(graph)
    fragments=[];exclusions=[];source_parts=[]
    for identity, part in enumerate(parts):
        path=longest_path(graph,part)
        xy=[[x/scale,y/scale] for y,x in path]
        length=sum(math.dist(a,b) for a,b in zip(xy,xy[1:]))
        bbox=[min(x for y,x in part)/scale,min(y for y,x in part)/scale,
              max(x for y,x in part)/scale,max(y for y,x in part)/scale]
        record={'ink_component':identity,'pixels':len(part),'retained_path_pixels':len(path),
                'bbox_pdf':bbox,'principal_length':length,
                'removed_branch_pixels':[[x/scale,y/scale] for y,x in sorted(part-set(path))]}
        if length<minimum:
            exclusions.append(record|{'reason':'principal_path_below_minimum_length'})
            continue
        source_parts.append(record)
        fragments.append({'id':len(fragments),'ink_component':identity,
                          'points':simplify(xy),'length':length,'bbox_pdf':bbox})
    if not fragments:raise ValueError('No curve fragments remain')
    endpoints=[]
    for fragment in fragments:
        for end in (0,1):
            p=fragment['points'] if not end else fragment['points'][::-1]
            endpoints.append({'id':2*fragment['id']+end,'fragment':fragment['id'],
                              'end':end,'position':p[0],'outward':outward(p)})
    edgegraph=nx.Graph();edgegraph.add_nodes_from(range(len(endpoints)))
    costs={}
    for i,a in enumerate(endpoints):
        for j,b in enumerate(endpoints[i+1:],i+1):
            if a['fragment']==b['fragment'] and len(fragments)>1:continue
            cost=gap_cost(a,b);costs[i,j]=cost
            edgegraph.add_edge(i,j,weight=cost)
    matching={}
    for a,b in nx.min_weight_matching(edgegraph):matching[a]=b;matching[b]=a
    if len(matching)!=len(endpoints):raise ValueError('No perfect endpoint matching')
    merges=[]
    while len(get_cycles(matching,len(fragments)))>1:
        cycles=get_cycles(matching,len(fragments));choice=None
        for ci,c1 in enumerate(cycles):
            for c2 in cycles[ci+1:]:
                for p in c1:
                    a,b=p^1,matching[p^1]
                    for q in c2:
                        c,d=q^1,matching[q^1]
                        for first,second in [((a,c),(b,d)),((a,d),(b,c))]:
                            if first[0]//2==first[1]//2 or second[0]//2==second[1]//2:continue
                            delta=(gap_cost(endpoints[first[0]],endpoints[first[1]])+
                                   gap_cost(endpoints[second[0]],endpoints[second[1]])-
                                   gap_cost(endpoints[a],endpoints[b])-gap_cost(endpoints[c],endpoints[d]))
                            candidate=(delta,a,b,c,d,first,second)
                            if choice is None or candidate<choice:choice=candidate
        if choice is None:raise ValueError('Cannot merge endpoint cycles')
        delta,a,b,c,d,first,second=choice
        for u,v in (first,second):matching[u]=v;matching[v]=u
        merges.append({'removed_pairs':[[a,b],[c,d]],'added_pairs':[list(first),list(second)],
                       'score_increase':delta})
    ordered=get_cycles(matching,len(fragments))[0];points=[];edge_kinds=[];joins=[];ranges=[]
    for entry in ordered:
        f=fragments[entry//2];p=f['points'] if entry%2==0 else f['points'][::-1]
        if points:
            points.append(p[0]);edge_kinds.append('gap')
        begin=len(points)
        if not points:points.append(p[0])
        for q in p[1:]:points.append(q);edge_kinds.append('ink')
        ranges.append({'fragment':f['id'],'orientation':1 if entry%2==0 else -1,
                       'first_point_index':begin-1 if begin else 0,'last_point_index':len(points)-1})
        exit_id=entry^1;next_id=matching[exit_id]
        a,b=endpoints[exit_id],endpoints[next_id]
        alternatives=sorted((gap_cost(a,z),z['id']) for z in endpoints if z['id']!=exit_id and z['fragment']!=a['fragment'])
        joins.append({'id':len(joins),'endpoint_ids':[exit_id,next_id],
                      'from_fragment':a['fragment'],'to_fragment':b['fragment'],
                      'from':a['position'],'to':b['position'],
                      'distance':math.dist(a['position'],b['position']),
                      'score':gap_cost(a,b),'source_is_chosen_best':bool(alternatives and alternatives[0][1]==next_id),
                      'closest_alternatives':alternatives[:3],'manual_audit_status':'PENDING'})
    points.append(points[0]);edge_kinds.append('gap')
    if len(edge_kinds)!=len(points)-1:raise ValueError('Edge kind count differs')
    report={'colour':colour,'status':'CLOSED_CANDIDATE_TRACE_REQUIRES_DRAWING_AUDIT','closed':True,
            'points_pdf':points,'edge_kinds':edge_kinds,'fragment_ranges':ranges,'joins':joins,
            'retained_fragments':fragments,'source_component_pruning':source_parts,
            'excluded_tiny_components':exclusions,'label_exclusions':exclusion_records,
            'cycle_merge_edits':merges,'layer_source_sha256':digest(image_path),
            'scale_pixels_per_point':scale,'threshold':threshold,'minimum_principal_length':minimum,
            'no_crossing_information_assigned':True,'drawing_equivalence_certified':False,
            'wall_seconds':time.monotonic()-started}
    (output/f'{colour}_closed.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    canvas=Image.open(source/'original_page.png').convert('RGB');draw=ImageDraw.Draw(canvas)
    for a,b,kind in zip(points,points[1:],edge_kinds):
        draw.line(tuple(v*scale for v in a+b),fill=(0,0,0) if kind=='ink' else (255,0,255),width=2 if kind=='gap' else 1)
    for join in joins:
        x,y=[(a+b)*scale/2 for a,b in zip(join['from'],join['to'])]
        draw.text((x+2,y+2),str(join['id']),fill=(0,0,0))
    canvas.save(output/f'{colour}_closed_overlay.png')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--figure',type=int,choices=[1,2],required=True)
    parser.add_argument('--colours',nargs='+',default=['purple','yellow','cyan','blue'])
    parser.add_argument('--scale',type=float,default=2);parser.add_argument('--threshold',type=int,default=20)
    parser.add_argument('--minimum',type=float,default=5.0)
    parser.add_argument('--reviewed-figure1-blue-hooks',action='store_true')
    args=parser.parse_args()
    if args.out.exists():raise ValueError('Use a new output directory')
    args.out.mkdir(parents=True);components={}
    for colour in args.colours:
        result=reconstruct(args.source,args.out,colour,args.figure,args.scale,args.threshold,args.minimum)
        if args.reviewed_figure1_blue_hooks and args.figure==1 and colour=='blue':
            result=reviewed_blue_cleanup(result,args.out,args.source,args.scale)
        components[colour]=result
        bundle={'schema':'motoo.closed-coloured-traces.v1','figure':args.figure,'components':components,
                'status':'CLOSED_TRACES_REQUIRE_DRAWING_AUDIT','pd_generated':False}
        (args.out/'traces.json').write_text(json.dumps(bundle,separators=(',',':'))+'\n')
        print(json.dumps({'colour':colour,'closed':True,'fragments':len(result['retained_fragments']),
                          'points':len(result['points_pdf']),'joins':len(result['joins']),
                          'cycle_merges':len(result['cycle_merge_edits']),
                          'largest_join':max(j['distance'] for j in result['joins']),
                          'pruned_pixels':sum(len(c['removed_branch_pixels']) for c in result['source_component_pruning']),
                          'seconds':result['wall_seconds']}),flush=True)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Targeted all-over/all-under arc isotopies, allowing crossing growth.

Uses the dual graph, strand deletion and crossing insertion primitives in
Spherogram 2.4.1. The replacement is restricted to a strict uniform-height
arc with distinct crossings and surviving endpoint crossings; no recursive
self-crossing extension and no hidden Reidemeister simplification is used.
Each accepted macro records both diagrams and the complete dual route.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import random
import time
import regina

import networkx as nx
from spherogram import Link
from spherogram.links import simplify as simp


def pd(link):return [list(t) for t in link.PD_code()]


def owner_map(link, source_order):
    return {(ep.crossing,ep.strand_index%2):source_order[i]
            for i,comp in enumerate(link.link_components) for ep in comp}


def metric(link,source_order):
    owners=owner_map(link,source_order);dots={1,2,4}
    dotdot=selfcross=dotincident=0
    for c in link.crossings:
        a,b=owners[c,0],owners[c,1]
        if a in dots or b in dots:dotincident+=1
        if a in dots and b in dots:
            if a==b:selfcross+=1
            else:dotdot+=1
    trans={}
    for i,comp in enumerate(link.link_components):
        if source_order[i] not in dots:continue
        kinds=[ep.strand_index%2 for ep in comp]
        trans[source_order[i]]=sum(a!=b for a,b in zip(kinds,kinds[1:]+kinds[:1]))
    excess=sum(abs(x-2) for x in trans.values())
    score=(selfcross+dotdot,dotincident,excess,len(link.crossings))
    return score,{'dotted_dotted':dotdot,'dotted_self':selfcross,'dot_incident':dotincident,'transitions':trans,'crossings':len(link.crossings)}


def candidates(link):
    ids={c:i for i,c in enumerate(link.crossings)}
    result=[]
    for kind in ['over','under']:
        for strand in simp.over_or_under_strands(link,kind):
            # Also allow proper subarcs: a maximal arc can have an endpoint
            # crossing within itself, whereas a shortened arc may be strict.
            variants=[strand]
            if len(strand)>2:variants += [strand[1:],strand[:-1]]
            for part in variants:
                crossset={ep.crossing for ep in part}
                if len(crossset)!=len(part):continue
                if part[0].previous().crossing in crossset or part[-1].next().crossing in crossset:continue
                result.append((kind,[(ids[ep.crossing],ep.strand_index) for ep in part]))
    return result


def reroute(link,source_order,spec,focus=None,prescribed_dual_path=None):
    before_pd=pd(link)
    beforeowners=owner_map(link,source_order)
    beforecross={c:i for i,c in enumerate(link.crossings)}
    beforelinking=link.linking_matrix()
    kind,description=spec
    from spherogram.links.links_base import CrossingEntryPoint
    strand=[CrossingEntryPoint(link.crossings[i],j) for i,j in description]
    if not strand or len({ep.crossing for ep in strand})!=len(strand):raise ValueError('Non-strict pickup arc')
    if any(not getattr(ep,'is_'+kind+'_crossing')() for ep in strand):raise ValueError('Mixed-height arc')
    if any(a.next()!=b for a,b in zip(strand,strand[1:])):raise ValueError('Disconnected pickup arc')
    start,end=strand[0].previous(),strand[-1].next()
    removed_set={ep.crossing for ep in strand}
    if start.crossing in removed_set or end.crossing in removed_set:raise ValueError('Endpoint crossing removed')
    picked_source=beforeowners[strand[0].crossing,strand[0].strand_index%2]
    G=simp.dual_graph_as_nx(link)
    crossed=simp.dual_edges(strand,G)
    start_face,dest_face=crossed[0][0],crossed[-1][0]
    nx.set_edge_attributes(G,1,'weight')
    for a,b in crossed:G[a][b]['weight']=0
    dots={1,2,4}
    def cost(a,b,e):
        if e['weight']==0:return 0
        ep=e['interface'][a]
        source=beforeowners[ep.crossing,ep.strand_index%2]
        return 10000 if source==focus else (1000 if source in dots else 1)
    path=nx.shortest_path(G,start_face,dest_face,weight=cost) if prescribed_dual_path is None else prescribed_dual_path
    if path[0]!=start_face or path[-1]!=dest_face or len(set(path))!=len(path):
        raise ValueError('Replacement dual route has wrong endpoints or repeats a face')
    if any(not G.has_edge(x,y) for x,y in zip(path,path[1:])):
        raise ValueError('Replacement dual route is not a graph path')
    interfaces=[]
    for a,b in zip(path,path[1:]):
        e=G[a][b]
        if e['weight']:
            ep=e['interface'][a]
            interfaces.append({'from_face':a,'to_face':b,'endpoint':[beforecross[ep.crossing],ep.strand_index]})
    # Give the independent checker the exact original endpoints and deletion
    # support before any objects are mutated.
    record={'kind':kind,'old_strand':[list(x) for x in description],
            'start_endpoint':[beforecross[start.crossing],start.strand_index],
            'end_endpoint':[beforecross[end.crossing],end.strand_index],
            'dual_path':list(path),'crossed_interfaces':interfaces,
            'picked_component_source':picked_source,
            'before_pd_code':before_pd,'before_component_to_source':list(source_order),
            'focus':focus,'strict_uniform_height_arc':True,'implicit_basic_simplification':False}
    simp.remove_strand(link,strand)
    loose=start.rotate(2);new=[]
    for a,b in zip(path,path[1:]):
        if G[a][b]['weight']:
            ep=G[a][b]['interface'][a]
            c,loose=simp.cross(link,ep,kind,loose,'pickup_'+str(len(new)))
            new.append(c)
    end.crossing[end.strand_index]=loose
    link.crossings.extend(new)
    if len(link.link_components)!=len(source_order):raise ValueError('Component disappeared during pickup')
    from spherogram.links.links_base import Crossing
    old_rotate=Crossing.rotate
    shifts={c:0 for c in link.crossings if c in beforecross}
    def tracked_rotate(c,s):
        if c in shifts:
            if s%2:raise ValueError('A surviving crossing rotated across its over/under axes')
            shifts[c]=(shifts[c]-s)%4
        return old_rotate(c,s)
    Crossing.rotate=tracked_rotate
    try:link._rebuild(same_components_and_orientations=True)
    finally:Crossing.rotate=old_rotate
    after_order=[]
    for comp in link.link_components:
        candidates={beforeowners[ep.crossing,ep.strand_index%2] for ep in comp
                    if (ep.crossing,ep.strand_index%2) in beforeowners}
        if len(candidates)!=1:raise ValueError('Cannot transport source component unambiguously')
        after_order.append(candidates.pop())
    if sorted(after_order)!=sorted(source_order):raise ValueError('Source component bijection lost')
    afterlinking=link.linking_matrix()
    for i,a in enumerate(after_order):
        for j,b in enumerate(after_order):
            if afterlinking[i][j]!=beforelinking[source_order.index(a)][source_order.index(b)]:
                raise ValueError('Oriented linking number changed')
    # Validate the resulting abstract embedding is planar, not a virtual move.
    Link(pd(link),check_planarity=True)
    record.update({'after_pd_code':pd(link),'after_component_to_source':after_order,
                   'new_crossings':len(new),'removed_crossings':sorted(beforecross[x] for x in removed_set),
                   'surviving_crossing_map':[[beforecross[c],i] for i,c in enumerate(link.crossings) if c in beforecross],
                   'surviving_crossing_port_shifts':[[beforecross[c],shifts[c]] for c in link.crossings if c in beforecross],
                   'new_crossing_indices':[i for i,c in enumerate(link.crossings) if c in new]})
    return after_order,record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seconds',type=float,default=120);p.add_argument('--seed',type=int,default=0)
    p.add_argument('--ceiling',type=int,default=500)
    p.add_argument('--walk',action='store_true',help='Allow bounded neutral/uphill isotopies, retaining the best witnessed prefix')
    p.add_argument('--objective',choices=['unlink','position'],default='unlink')
    a=p.parse_args();raw=a.source.read_bytes();d=json.loads(raw)
    link=Link(d['pd_code']);order=d.get('pd_component_to_source',list(range(8)))
    def scoring(L,O):
        default,mm=metric(L,O)
        if a.objective=='position':
            return (mm['dotted_dotted'],mm['dotted_self'],sum(abs(v-2) for v in mm['transitions'].values()),mm['dot_incident'],mm['crossings']),mm
        return default,mm
    rng=random.Random(a.seed);start=time.monotonic();log=[];score,metrics=scoring(link,order);attempts=0
    best=copy.deepcopy(link);best_order=list(order);best_log=[];best_score=score;best_metrics=metrics
    def search_key(L):
        # This is a heuristic visited-state key, not an exhaustive marked-link
        # classification. Source component labels remain in every proof record.
        R=regina.Link.fromPD([[x+1 for x in row] for row in pd(L)])
        return R.sig(False,False,False)
    seen={search_key(link)}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    event_path=a.out.with_suffix(a.out.suffix+'.events.jsonl')
    event_stream=event_path.open('w')
    from check_pickup_macros import check_macro
    def bundle(L,O,records,mm):
        return {'schema':'motoo.strict_pickup_macros.v1','source_json_sha256':hashlib.sha256(raw).hexdigest(),
                'source_file':a.source.name,'macros':records,'pd_code':pd(L),'pd_component_to_source':O,
                'decorations':[0 if i in [0,5,6,7] else (-1 if i==3 else 'x') for i in O],
                'metrics':mm,'elapsed_seconds':time.monotonic()-start,'candidate_attempts':attempts,
                'walk_mode':a.walk,'objective':a.objective,'heuristic_states_seen':len(seen),
                'independent_local_tangle_checker_applied_to_every_accepted_macro':True,
                'status':'PRODUCER_RESULT_PENDING_FRESH_CHAIN_REPLAY'}
    def save_checkpoint(L,O,records,mm):
        temporary=a.out.with_suffix(a.out.suffix+'.tmp')
        temporary.write_text(json.dumps(bundle(L,O,records,mm),indent=2)+'\n')
        temporary.replace(a.out)
    save_checkpoint(best,best_order,best_log,best_metrics)
    print(json.dumps({'event':'start','figure':d.get('figure'),'score':score,'metrics':metrics}),flush=True)
    while time.monotonic()-start<a.seconds:
        options=candidates(link);rng.shuffle(options)
        accepted=None
        for spec in options:
            for focus in [None,1,2,4]:
                attempts+=1;trial=copy.deepcopy(link)
                try:new_order,record=reroute(trial,order,spec,focus)
                except ValueError as error:
                    # A rejected geometric candidate is diagnostic, never an
                    # accepted move. The original link is untouched.
                    continue
                if len(trial.crossings)>a.ceiling:continue
                newscore,newmetrics=scoring(trial,new_order)
                if a.walk:
                    if newmetrics['dotted_dotted']+newmetrics['dotted_self']>best_metrics['dotted_dotted']+best_metrics['dotted_self']+2:continue
                    if newmetrics['dot_incident']>best_metrics['dot_incident']+12 or newmetrics['crossings']>best_metrics['crossings']+80:continue
                    key=search_key(trial)
                    if key in seen:continue
                    try:check_macro(record)
                    except ValueError:continue
                    seen.add(key)
                    accepted=(trial,new_order,record,newscore,newmetrics);break
                elif newscore<score:
                    try:check_macro(record)
                    except ValueError:continue
                    accepted=(trial,new_order,record,newscore,newmetrics);break
                if time.monotonic()-start>=a.seconds:break
            if accepted or time.monotonic()-start>=a.seconds:break
        if accepted is None:break
        link,order,record,score,metrics=accepted;log.append(record)
        event_stream.write(json.dumps({'sequence':len(log),'macro':record},separators=(',',':'))+'\n')
        event_stream.flush()
        if score<best_score:
            best=copy.deepcopy(link);best_order=list(order);best_log=list(log);best_score=score;best_metrics=metrics
            save_checkpoint(best,best_order,best_log,best_metrics)
            print(json.dumps({'event':'improved','elapsed':time.monotonic()-start,'score':score,'metrics':metrics,'macros':len(log)}),flush=True)
        if metrics['dotted_dotted']==metrics['dotted_self']==0 and all(v==2 for v in metrics['transitions'].values()):break
    link,order,log,score,metrics=best,best_order,best_log,best_score,best_metrics
    save_checkpoint(link,order,log,metrics)
    event_stream.close()
    print(json.dumps({'event':'finish','metrics':metrics,'macros':len(log),'file':str(a.out)}),flush=True)


if __name__=='__main__':main()

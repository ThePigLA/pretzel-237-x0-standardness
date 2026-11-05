#!/usr/bin/env python3
"""Bounded, logged Reidemeister simplification of the marked full link.

This lane never adds crossings. It does not certify counterclockwise planar
orientation for Katie; that is a separate obligation after normalization.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import time

import regina


def state(link):
    components=[]
    for ci in range(link.countComponents()):
        first=link.component(ci)
        visits=[]
        if first:
            ref=first
            while True:
                visits.append((ref.crossing().index()+1)*(1 if ref.strand() else -1))
                ref=ref.next()
                if ref==first: break
        if visits:
            k=min(range(len(visits)),key=lambda i:(abs(visits[i]),visits[i]))
            visits=visits[k:]+visits[:k]
        components.append(visits)
    return {'signs':[link.crossing(i).sign() for i in range(link.size())], 'components':components}


def digest(s):
    return hashlib.sha256(json.dumps(s,separators=(',',':')).encode()).hexdigest()


def metrics(link,decorations):
    s=state(link)
    owners={}
    for c,seq in enumerate(s['components']):
        for x in seq: owners[(abs(x)-1, x>0)]=c
    dots={i for i,x in enumerate(decorations) if x=='x'}
    selfcross=0; dotdot=0
    for i in range(link.size()):
        a,b=owners[i,False],owners[i,True]
        if a in dots and b in dots:
            if a==b:selfcross+=1
            else:dotdot+=1
    transitions=[]
    for i in sorted(dots):
        seq=s['components'][i]
        transitions.append(sum((seq[j]>0)!=(seq[(j+1)%len(seq)]>0) for j in range(len(seq))) if seq else 0)
    # Strong, sufficient diagram condition. Entirely overpassing or empty dots
    # do not qualify as having two marked transition points.
    excess=sum(abs(t-2) for t in transitions)
    return (dotdot,selfcross,excess,link.size()),{'dotted_dotted_crossings':dotdot,'dotted_self_crossings':selfcross,'dotted_transitions':transitions,'crossings':link.size()}


def perform(link,operation):
    before=state(link)
    kind,i,*rest=operation
    if kind=='R1-':
        removed={i}
        success=link.r1(link.crossing(i))
    elif kind=='R2-':
        removed={i,link.crossing(i).upper().next().crossing().index()}
        success=link.r2(link.crossing(i))
    elif kind=='R3':
        removed=set()
        success=link.r3(link.crossing(i),rest[0])
    else:raise ValueError(operation)
    if not success:raise ValueError('Guarded Reidemeister move rejected')
    after=state(link)
    if len(before['components'])!=len(after['components']):raise ValueError('Component count changed')
    survivor_map={old:new for new,old in enumerate(j for j in range(len(before['signs'])) if j not in removed)}
    if after['signs'] != [v for j,v in enumerate(before['signs']) if j not in removed]:raise ValueError('Unexpected crossing reindexing/sign changes')
    for a,b in zip(before['components'],after['components']):
        expected=[(survivor_map[abs(x)-1]+1)*(1 if x>0 else -1) for x in a if abs(x)-1 not in removed]
        if Counter(expected)!=Counter(b):raise ValueError('A marked component changed its surviving strands')
        if kind!='R3' and expected:
            k=min(range(len(expected)),key=lambda i:(abs(expected[i]),expected[i]))
            expected=expected[k:]+expected[:k]
            if expected!=b:raise ValueError('A deletion changed component cyclic order')
    return {'op':operation,'before_sha256':digest(before),'after_sha256':digest(after),'removed_crossings':sorted(removed)}


def descend(link,log):
    while True:
        chosen=None
        for i in range(link.size()):
            if link.hasR1(link.crossing(i)):
                chosen=['R1-',i];break
            if link.hasR2(link.crossing(i)):
                chosen=['R2-',i];break
        if chosen is None:return
        log.append(perform(link,chosen))


def replay(initial,moves,endpoint):
    link=regina.Link.fromData(initial['signs'],initial['components'])
    for i,record in enumerate(moves):
        observed=perform(link,record['op'])
        if observed!=record:raise ValueError(f'Replay mismatch at move{i}')
    if state(link)!=endpoint:raise ValueError('Replay endpoint mismatch')
    return link


def search(path,out,seconds,seed):
    raw=path.read_bytes();source=json.loads(raw)
    link=regina.Link.fromPD(source.get('regina_pd_code',[[x+1 for x in row] for row in source['pd_code']]))
    initial=state(link);decorations=source['decorations']
    log=[];rng=random.Random(seed)
    initial_score,initial_metrics=metrics(link,decorations)
    descend(link,log)
    best_score,best_metrics=metrics(link,decorations)
    best=regina.Link(link);bestlog=list(log)
    seen={digest(state(link))};start=time.monotonic();attempts=0;restarts=0
    figure=source.get('figure',source.get('source_figure'))
    print(json.dumps({'event':'start','figure':figure,'score':best_score,'metrics':best_metrics}),flush=True)
    while time.monotonic()-start<seconds and best_score[:3]!=(0,0,0):
        moves=[['R3',i,side] for i in range(link.size()) for side in [0,1] if link.hasR3(link.crossing(i),side)]
        rng.shuffle(moves)
        chosen=None
        for move in moves:
            trial=regina.Link(link)
            record=perform(trial,move)
            if record['after_sha256'] in seen:continue
            chosen=(trial,record);break
        if chosen is None:
            restarts+=1;link=regina.Link(best);log=list(bestlog);seen={digest(state(link))}
            if not moves:break
            if restarts>200:break
            continue
        link,record=chosen;log.append(record);seen.add(record['after_sha256']);attempts+=1
        n=link.size();descend(link,log)
        if link.size()<n:seen={digest(state(link))}
        score,mm=metrics(link,decorations)
        if score<best_score:
            best=regina.Link(link);bestlog=list(log);best_score,best_metrics=score,mm
            print(json.dumps({'event':'best','figure':figure,'elapsed':time.monotonic()-start,'score':score,'moves':len(log),'metrics':mm}),flush=True)
        if len(seen)>5000:
            restarts+=1;link=regina.Link(best);log=list(bestlog);seen={digest(state(link))}
    endpoint=state(best)
    replay(initial,bestlog,endpoint)
    result={'schema':'motoo.marked_reidemeister.v1','source_figure':figure,
            'source_json_sha256':hashlib.sha256(raw).hexdigest(),'source_path':path.name,
            'decorations':decorations,'source_component_identity':source.get('pd_component_to_source',source.get('source_component_identity',list(range(8)))),
            'initial_state':initial,'endpoint_state':endpoint,'moves':bestlog,
            'initial_metrics':initial_metrics,'endpoint_metrics':best_metrics,
            'elapsed_seconds':time.monotonic()-start,'seed':seed,'r3_attempts':attempts,'restarts':restarts,
            'all_guards_replayed':True,'marked_component_membership_checked_every_move':True,
            'status':'DOTTED_PROJECTION_NORMALIZED_CCW_NOT_YET_CHECKED' if best_score[:3]==(0,0,0) else 'NOT_STANDARD_POSITION',
            'integer_framings_preserved':True}
    if all(endpoint['components']):
        result['regina_pd_code']=best.pdData()
        result['pd_code']=[[a-1 for a in row] for row in result['regina_pd_code']]
    else:result['pd_export_omitted_for_zero_crossing_components']=True
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'event':'finish','status':result['status'],'score':best_score,'moves':len(bestlog),'file':str(out)}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seconds',type=float,default=60);p.add_argument('--seed',type=int,default=0)
    a=p.parse_args();search(a.source,a.out,a.seconds,a.seed)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Erase exact labelled-state loops from a finished ideal move certificate."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
from reduce_ideal_checked import require,save_json,sha_file,move_key,load_tri,state_hash


def shorten(source_directory,out):
    p=Path(source_directory);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    final=json.loads((p/'FINAL.json').read_text());md=json.loads((p/'SOURCE.json').read_text())
    require(final['status']=='PASS_EXACT_DETERMINISTIC_TARGET','input is not a finished retrace')
    require(sha_file(p/'moves.jsonl')==final['log_sha256'],'input log binding differs')
    stack=[];states={md['start_tight_sha256']:0};n=0;current=md['start_tight_sha256'];removed=0
    for line in (p/'moves.jsonl').read_text().splitlines():
        row=json.loads(line)
        require(row['index']==n and row['before_sha256']==current,'input chain discontinuity')
        current=row['after_sha256'];n+=1
        if current in states:
            pos=states[current]
            for old in stack[pos:]:states.pop(old['after_sha256'],None)
            removed+=len(stack)-pos+1;stack=stack[:pos]
        else:stack.append(row);states[current]=len(stack)
    require(n==final['moves'] and current==final['current_tight_sha256'],'incomplete input log')
    require(len(stack)+removed==n,'loop accounting differs')
    counts=Counter();current=md['start_tight_sha256'];original_indices=[]
    temporary=out/'moves.tmp'
    with temporary.open('w') as output:
        for i,row in enumerate(stack):
            require(row['before_sha256']==current,'shortened chain discontinuity')
            original_indices.append(row['index'])
            copy=dict(row);copy['original_event_index']=copy['index'];copy['index']=i
            output.write(json.dumps(copy,separators=(',',':'))+'\n')
            current=row['after_sha256'];counts[move_key(row['operation'])]+=1
    require(current==final['current_tight_sha256'],'shortened endpoint differs')
    temporary.replace(out/'moves.jsonl')
    endpoint=p/final['endpoint_file'];require(sha_file(endpoint)==final['endpoint_rga_sha256'],'input RGA binding differs')
    require(state_hash(load_tri(endpoint))==current,'input RGA state differs')
    shutil.copy2(endpoint,out/'endpoint.rga');save_json(out/'SOURCE.json',md)
    report={'status':'PASS_EXACT_LABELLED_LOOP_ERASURE','input_events':n,'output_events':len(stack),
        'removed_events':removed,'source_sha256':md['source_sha256'],
        'input_log_sha256':final['log_sha256'],'input_manifest_sha256':sha_file(p/'FINAL.json'),
        'output_log_sha256':sha_file(out/'moves.jsonl'),'endpoint_rga_sha256':sha_file(out/'endpoint.rga'),
        'endpoint_tight_sha256':current,'pentachora':final['current_pentachora'],
        'counts':dict(counts),'guarded_replay_of_shortened_path_pending':True,
        'standardness_asserted':False,'closed_completion_asserted':False}
    save_json(out/'SHORTENING.json',report)
    save_json(out/'ORIGINAL_EVENT_INDICES.json',original_indices)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();shorten(a.directory,a.out)

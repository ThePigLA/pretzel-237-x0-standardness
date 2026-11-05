#!/usr/bin/env python3
"""Rebuild both ideal checkpoint files from the committed guarded move prefix.

Preserve any mismatched original RGA files for diagnosis.  This never changes
the move log, its endpoint hashes, or the source binding.
"""
import argparse
import json
from pathlib import Path
import time
from collections import Counter

from reduce_ideal_checked import (basic,load_tri,perform,require,save_json,
                                  save_tri,sha_file,state_hash,move_key)


def export(directory,source_override=None):
    directory=Path(directory)
    md=json.loads((directory/'SOURCE.json').read_text())
    cp=json.loads((directory/'CHECKPOINT.json').read_text())
    source=Path(source_override) if source_override is not None else Path(md['source_file'])
    require(sha_file(source)==md['source_sha256'],'source file hash differs')
    started=time.monotonic();tri=load_tri(source);initial=basic(tri)
    current=state_hash(tri);counts=Counter();n=0
    before={}
    for name in ('current.rga','best.rga'):
        p=directory/name
        if p.exists():
            t=load_tri(p)
            before[name]={'file_sha256':sha_file(p),'tight_sha256':state_hash(t),'pentachora':t.size()}
    if cp['best_moves']==0:
        require(current==cp['best_tight_sha256'],'initial best differs')
        save_tri(directory/'best.replayed.rga',tri)
    with (directory/'moves.jsonl').open() as stream:
        for line in stream:
            if n==cp['moves']:break
            row=json.loads(line)
            require(row['index']==n and row['before_sha256']==current,'prefix chain discontinuity')
            perform(tri,row['operation']);current=state_hash(tri);n+=1
            require(row['after_sha256']==current and row['pentachora']==tri.size(),'guarded replay differs')
            counts[move_key(row['operation'])]+=1
            if n==cp['best_moves']:
                require(current==cp['best_tight_sha256'],'best prefix hash differs')
                save_tri(directory/'best.replayed.rga',tri)
            if n%1000==0:print(json.dumps({'moves':n,'pentachora':tri.size(),'seconds':time.monotonic()-started}),flush=True)
    require(n==cp['moves'] and current==cp['current_tight_sha256'],'final prefix hash differs')
    require(dict(counts)==cp['counts'],'move counts differ')
    save_tri(directory/'current.replayed.rga',tri)
    after={}
    for stem,key in [('best','best_tight_sha256'),('current','current_tight_sha256')]:
        p=directory/(stem+'.replayed.rga');t=load_tri(p)
        require(state_hash(t)==cp[key],'saved RGA round-trip state differs')
        after[stem+'.rga']={'file_sha256':sha_file(p),'tight_sha256':state_hash(t),'pentachora':t.size()}
    # Only after both exact round-trips pass do we replace public checkpoints.
    for stem in ('best','current'):
        original=directory/(stem+'.rga')
        if original.exists():
            saved=directory/('original_before_replay_'+stem+'.rga')
            require(not saved.exists(),'a previous repair exists; do not overwrite evidence')
            original.replace(saved)
        (directory/(stem+'.replayed.rga')).replace(original)
        cp[stem+'_rga_sha256']=sha_file(original)
    cp['status']='FROZEN_VERIFIED_GUARDED_PREFIX'
    cp['guarded_replay_complete']=True
    save_json(directory/'CHECKPOINT.json',cp)
    receipt={'status':'PASS_REBUILT_CURRENT_AND_BEST_FROM_GUARDED_PREFIX',
             'source_sha256':sha_file(source),'move_log_sha256':sha_file(directory/'moves.jsonl'),
             'moves':n,'best_moves':cp['best_moves'],'before':before,'after':after,
             'source_basic':initial,'seconds':time.monotonic()-started,
             'standardness_asserted':False,'closed_completion_asserted':False}
    save_json(directory/'REPAIR_RECEIPT.json',receipt)
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    p.add_argument('--source',type=Path)
    a=p.parse_args();export(a.directory,a.source)

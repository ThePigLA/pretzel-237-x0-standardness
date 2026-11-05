#!/usr/bin/env python3
"""Deterministically retrace the ideal burst search into immutable log chunks.

This avoids keeping append-file descriptors open across workspace updates.
The initial prefix and RNG state are taken from a repaired frozen checkpoint.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import signal
import time

from reduce_ideal_checked import (basic,eligible,load_tri,move_key,perform,
    replay,require,save_json,save_tri,sha_file,state_hash,tuple_tree)


def run(frozen,out,target,expected,chunk_size,source_override=None):
    frozen=Path(frozen);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    (out/'chunks').mkdir()
    source_data=json.loads((frozen/'SOURCE.json').read_text())
    cp=json.loads((frozen/'CHECKPOINT.json').read_text())
    source=Path(source_override) if source_override is not None else Path(source_data['source_file'])
    require(sha_file(source)==source_data['source_sha256'],'source binding differs')
    report,tri=replay(source,frozen/'moves.jsonl',cp['moves'],cp['current_tight_sha256'])
    save_json(out/'START_REPLAY.json',report);save_json(out/'SOURCE.json',source_data)
    n=cp['moves'];current=cp['current_tight_sha256'];best_size=cp['best_pentachora']
    best_hash=cp['best_tight_sha256'];best_move=cp['best_moves']
    counts=Counter(cp['counts']);neutral=cp.get('neutral_since_reduction',0)
    growth_remaining=cp.get('growth_remaining',0);mix_remaining=cp.get('mix_remaining',0)
    growth_cycles=cp.get('growth_cycles',0);rng=random.Random()
    rng.setstate(tuple_tree(cp['rng_state']))
    growth=80;neutral_limit=100;mix_steps=500
    started=time.monotonic();last_print=started;stop=False
    chunks=[];buffer=[];buffer_start=n
    prefix=b''.join((frozen/'moves.jsonl').read_bytes().splitlines(keepends=True)[:n])
    initial=out/'chunks'/f'{0:09d}_{n:09d}.jsonl'
    temporary=initial.with_suffix('.tmp');temporary.write_bytes(prefix);temporary.replace(initial)
    chunks.append({'file':str(initial.relative_to(out)),'start':0,'end':n,'sha256':sha_file(initial)})
    save_json(out/'POLICY.json',{'source_checkpoint':str(frozen.resolve()),
        'source_checkpoint_sha256':sha_file(frozen/'CHECKPOINT.json'),
        'target_moves':target,'target_tight_sha256':expected,'growth':growth,
        'neutral_limit':neutral_limit,'mix_steps':mix_steps,'chunk_size':chunk_size,
        'code_sha256':sha_file(__file__),'helper_code_sha256':sha_file(Path(__file__).with_name('reduce_ideal_checked.py')),
        'claim':'deterministic retrace; no standardness or closed completion assertion'})

    def seal():
        nonlocal buffer,buffer_start
        if buffer:
            p=out/'chunks'/f'{buffer_start:09d}_{n:09d}.jsonl'
            require(not p.exists(),'immutable chunk already exists')
            q=p.with_suffix('.tmp');q.write_text(''.join(buffer));q.replace(p)
            chunks.append({'file':str(p.relative_to(out)),'start':buffer_start,'end':n,'sha256':sha_file(p)})
            buffer=[];buffer_start=n
        checkpoint={'moves':n,'current_tight_sha256':current,'current_pentachora':tri.size(),
            'best_moves':best_move,'best_pentachora':best_size,'best_tight_sha256':best_hash,
            'counts':dict(counts),'rng_state':rng.getstate(),'neutral_since_reduction':neutral,
            'growth_remaining':growth_remaining,'mix_remaining':mix_remaining,
            'growth_cycles':growth_cycles,'chunks':chunks,'seconds':time.monotonic()-started}
        # The manifest is replaced only after its immutable chunks are closed.
        save_json(out/'SEALED.json',checkpoint)
        return checkpoint

    def commit(op):
        nonlocal n,current,best_size,best_hash,best_move,last_print
        before=current;perform(tri,op);current=state_hash(tri)
        row={'index':n,'before_sha256':before,'operation':op,'after_sha256':current,
             'pentachora':tri.size(),'elapsed_seconds':time.monotonic()-started}
        buffer.append(json.dumps(row,separators=(',',':'))+'\n')
        n+=1;counts[move_key(op)]+=1
        if tri.size()<best_size:best_size=tri.size();best_hash=current;best_move=n
        if len(buffer)>=chunk_size:seal()
        if time.monotonic()-last_print>=10:
            print(json.dumps({'moves':n,'current':tri.size(),'best':best_size,'seconds':time.monotonic()-started}),flush=True)
            last_print=time.monotonic()

    def signal_stop(signum,frame):
        nonlocal stop
        stop=True
    signal.signal(signal.SIGINT,signal_stop);signal.signal(signal.SIGTERM,signal_stop)
    seal()
    while n<target and not stop:
        if growth_remaining:
            options=eligible(tri,'P',3)
            if options:
                growth_remaining-=1;commit(rng.choice(options));continue
            growth_remaining=0
        if mix_remaining:
            options=eligible(tri,'P',2)
            if options:
                mix_remaining-=1;commit(rng.choice(options));continue
            mix_remaining=0
        options=[]
        for kind,k in [('Z',2),('Z',1),('P',1),('P',0),('Z',0)]:
            options=eligible(tri,kind,k)
            if options:break
        reducing=bool(options)
        if not options:
            if neutral>=neutral_limit:
                if growth and tri.size()+2<=best_size+growth:
                    growth_remaining=max(1,(best_size+growth-tri.size())//2)
                    mix_remaining=mix_steps;growth_cycles+=1;neutral=0;continue
                raise ValueError('unexpected local plateau before replay target')
            options=eligible(tri,'P',2);neutral+=1
        require(options,'no allowed move before replay target')
        op=rng.choice(options)
        if reducing:neutral=0
        commit(op)
    final=seal()
    require(n==target,'stopped before target; immutable sealed prefix retained')
    require(current==expected,'deterministic target state hash differs')
    final['status']='PASS_EXACT_DETERMINISTIC_TARGET'
    endpoint=out/f'endpoint_{tri.size()}_{current[:12]}.rga'
    save_tri(endpoint,tri)
    require(state_hash(load_tri(endpoint))==current,'saved endpoint round trip differs')
    full=out/'moves.jsonl';temporary=full.with_suffix('.tmp')
    count=0
    with temporary.open('wb') as output:
        for chunk in chunks:
            p=out/chunk['file'];require(sha_file(p)==chunk['sha256'],'immutable chunk changed')
            with p.open('rb') as inp:
                for line in inp:output.write(line);count+=1
    require(count==target,'assembled log length differs');temporary.replace(full)
    final.update({'endpoint_file':endpoint.name,'endpoint_rga_sha256':sha_file(endpoint),
        'source_sha256':sha_file(source),'log_sha256':sha_file(full),
        'standardness_asserted':False,'closed_completion_asserted':False})
    save_json(out/'FINAL.json',final)
    print(json.dumps({k:v for k,v in final.items() if k not in ('rng_state','chunks')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('frozen',type=Path)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--target-moves',type=int,required=True)
    p.add_argument('--expected-hash',required=True);p.add_argument('--chunk-size',type=int,default=100)
    p.add_argument('--source',type=Path)
    a=p.parse_args();run(a.frozen,a.out,a.target_moves,a.expected_hash,a.chunk_size,a.source)

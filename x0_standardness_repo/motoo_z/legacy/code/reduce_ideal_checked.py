#!/usr/bin/env python3
"""Stream guarded Pachner/2-0/edge-snap reductions of a bound ideal 4D input.

No native simplify or edge collapse is used.  Regina's skeleton builder does
implicitly recognize vertex links; this is not a separate link certificate.
This records a path of local Regina moves; it does not assert source admission
as a closed manifold, a cap, or standardness.  Tight-state hashes preserve
simplex labels for deterministic replay.  All guards remain active with -O.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import signal
import time

import regina


def require(value,message):
    if not value: raise ValueError(message)


def sha_file(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def state_hash(tri): return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()


def save_json(path,data):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(data,indent=2)+'\n');q.replace(p)


def save_tri(path,tri):
    p=Path(path);q=p.with_suffix('.tmp.rga')
    packet=regina.PacketOfTriangulation4(tri);packet.setLabel('Guarded ideal lower-body checkpoint; no cap or standardness assertion')
    require(packet.save(str(q)),'could not save checkpoint');q.replace(p)


def load_tri(path):
    packet=regina.open(str(path));require(isinstance(packet,regina.Triangulation4),'not a 4D triangulation')
    return regina.Triangulation4(packet)


def basic(tri):
    require(tri.size()>0 and tri.isValid() and tri.isConnected() and tri.isOrientable(),
            'input/output fails basic combinatorial validity')
    return {'pentachora':tri.size(),'vertices':tri.countVertices(),'valid':True,
            'connected':True,'orientable':True,'real_boundary_facets':tri.countBoundaryFacets(),
            'vertex_link_recognition_explicitly_requested':False,
            'implicit_regina_skeleton_recognition':True,'closed_manifold_asserted':False}


def eligible(tri,kind,dimension,protect_first=False):
    if dimension==4:
        require(kind=='P','only Pachner moves can use a top-dimensional simplex')
        return [['P',4,s.index()] for s in tri.simplices()
                if (not protect_first or s.index()!=0) and tri.hasPachner(s)]
    degree=2 if kind=='Z' else 5-dimension
    guard=tri.has20 if kind=='Z' else tri.hasPachner
    # Test the degree before guards that might otherwise inspect a large star.
    return [[kind,dimension,f.index()] for f in tri.faces(dimension)
            if f.degree()==degree and
            (not protect_first or all(e.simplex().index()!=0 for e in f.embeddings()))
            and guard(f)]


def perform(tri,op):
    if op[0]=='R':
        require(len(op)==2 and len(op[1])==tri.size(),'bad relabelling length')
        require(sorted(row[0] for row in op[1])==list(range(tri.size())),
                'simplex images are not bijective')
        iso=regina.Isomorphism4(tri.size())
        for i,(dest,perm) in enumerate(op[1]):
            require(sorted(perm)==list(range(5)),'vertex permutation is not bijective')
            iso.setSimpImage(i,dest);iso.setFacetPerm(i,regina.Perm5(perm))
        iso.applyInPlace(tri)
        return
    if op[0]=='S':
        require(len(op)==2 and 0<=op[1]<tri.countEdges(),'bad snap edge index')
        e=tri.edge(op[1]);old=tri.size();old_vertices=tri.countVertices()
        require(tri.hasSnapEdge(e),'edge snap guard rejected')
        require(tri.snapEdge(e),'guarded edge snap rejected')
        require(tri.size()==old+4 and tri.countVertices()==old_vertices-1,
                'unexpected snap size or vertex change')
        return
    require(len(op)==3 and op[0] in ('P','Z'),'unsupported operation')
    kind,k,index=op
    require(k in range(5) and (kind!='Z' or k<=2),'bad move dimension')
    require(0<=index<(tri.size() if k==4 else tri.countFaces(k)),'face index out of range')
    face=tri.simplex(index) if k==4 else tri.face(k,index)
    old=tri.size()
    guard=tri.has20 if kind=='Z' else tri.hasPachner
    move=tri.move20 if kind=='Z' else tri.pachner
    require(guard(face),'move guard rejected')
    require(move(face),'guarded move rejected')
    require(tri.size()==old+(-2 if kind=='Z' else 2*k-4),'wrong size change')


def replay(source,log,limit=None,expected=None):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina 7.4.1')
    tri=load_tri(source);initial=basic(tri);current=state_hash(tri);n=0;counts=Counter()
    started=time.monotonic()
    with Path(log).open() as stream:
        for line in stream:
            if limit is not None and n>=limit:break
            row=json.loads(line)
            require(row['index']==n,'log event numbering mismatch')
            require(row['before_sha256']==current,'state chain discontinuity')
            perform(tri,row['operation']);current=state_hash(tri)
            require(current==row['after_sha256'],'replayed state differs')
            require(tri.size()==row['pentachora'],'replayed size differs')
            counts[move_key(row['operation'])]+=1;n+=1
            if n%50==0: print(json.dumps({'replay_moves':n,'pentachora':tri.size(),'seconds':time.monotonic()-started}),flush=True)
    if limit is not None:require(n==limit,'log has too few moves')
    if expected is not None:require(current==expected,'wrong expected endpoint')
    return {'status':'PASS_GUARDED_IDEAL_PATH','moves':n,'counts':dict(counts),
            'source_sha256':sha_file(source),'source_basic':initial,'end_basic':basic(tri),
            'end_tight_sha256':current,'elapsed_seconds':time.monotonic()-started,
            'standardness_asserted':False,'source_interpretation':'Figure2 red over dashed green back arc; working assumption'},tri


def move_key(op):
    return op[0] if op[0] in ('R','S') else str(op[:2])


def tuple_tree(value):
    return tuple(tuple_tree(x) for x in value) if isinstance(value,list) else value


def publish_closed_log(path,expected_count,expected_hash):
    """Replace a closed append log atomically for cross-process readers."""
    path=Path(path);temporary=path.with_name(path.name+'.published.tmp')
    count=0;last=None
    with path.open('rb') as source,temporary.open('wb') as target:
        for line in source:
            target.write(line);count+=1;last=line
        target.flush();os.fsync(target.fileno())
    require(count==expected_count,'closed log does not contain the committed prefix')
    if expected_count:require(json.loads(last)['after_sha256']==expected_hash,'closed log endpoint differs')
    temporary.replace(path)


def search(source,out,seconds,seed,neutral_limit,growth,ideal_first_color=None,resume=False,snap_to_one=False,mix_steps=200,reseed=None):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina7.4.1')
    out=Path(out);started=time.monotonic();previous_elapsed=0
    rng=random.Random(seed);n=0;counts=Counter();neutral=0;mode='reduce'
    growth_remaining=0;mix_remaining=0;growth_cycles=0
    if resume:
        metadata=json.loads((out/'SOURCE.json').read_text())
        old=json.loads((out/'CHECKPOINT.json').read_text())
        require(sha_file(source)==metadata['source_sha256'],'resume source binding differs')
        lines=(out/'moves.jsonl').read_text().splitlines(keepends=True)
        require(len(lines)>=old['moves'],'checkpoint log is incomplete')
        if len(lines)>old['moves']:
            # Preserve any tail that was flushed before the atomic checkpoint.
            tail=out/('uncheckpointed_tail_'+str(time.time_ns())+'.jsonl')
            tail.write_text(''.join(lines[old['moves']:]))
            (out/'moves.jsonl').write_text(''.join(lines[:old['moves']]))
        # Replaying reconstructs Regina's valid, topology-preserving caches.
        # A fresh RGA reload cannot retain its known ideal-vertex count.
        verified,tri=replay(source,out/'moves.jsonl',old['moves'],old['current_tight_sha256'])
        save_json(out/'RESUME_REPLAY.json',verified)
        source_info=metadata['source_basic'];source_hash=metadata['source_sha256']
        current=old['current_tight_sha256'];best_size=old['best_pentachora']
        best_hash=old['best_tight_sha256'];best_move=old['best_moves']
        n=old['moves'];counts=Counter(old['counts']);max_size=old['max_pentachora']
        neutral=old.get('neutral_since_reduction',0)
        growth_remaining=old.get('growth_remaining',0);mix_remaining=old.get('mix_remaining',0)
        growth_cycles=old.get('growth_cycles',0)
        ideal_first_color=metadata['policy'].get('ideal_first_color')
        snap_to_one=metadata['policy'].get('snap_to_one_ideal',False)
        previous_elapsed=old['elapsed_seconds']
        if 'rng_state' in old:rng.setstate(tuple_tree(old['rng_state']))
        else:
            # Legacy checkpoints predate RNG-state recording; the path still
            # resumes exactly, but the subsequent heuristic choices are new.
            rng.seed(metadata['seed']+n*4294967297)
        if reseed is not None:rng.seed(reseed)
        started=time.monotonic()
    else:
        out.mkdir(parents=True,exist_ok=False)
        tri=load_tri(source);source_info=basic(tri)
        source_hash=sha_file(source);current=state_hash(tri);best_size=tri.size();best_hash=current;best_move=0
        max_size=tri.size()
        metadata={'schema':'motoo.guarded-ideal-path.v2','source_file':str(Path(source).resolve()),
              'source_sha256':source_hash,'source_basic':source_info,'start_tight_sha256':current,
              'regina_distribution':'7.4.1','seed':seed,
              'policy':{'neutral_limit':neutral_limit,'growth_ceiling':growth,'ideal_first_color':ideal_first_color,
                        'snap_to_one_ideal':snap_to_one,'mix_steps':mix_steps},
              'source_interpretation':'Figure2 red over dashed green back arc; adopted working assumption',
              'source_author_confirmation_received':False,'native_simplify_used':False,
              'edge_collapse_used':False,'closed_completion_asserted':False,'standardness_asserted':False}
        save_json(out/'SOURCE.json',metadata);save_tri(out/'best.rga',tri);save_tri(out/'current.rga',tri)
    last_update=started;stopping=False
    stop_file=out/'STOP'
    if resume and stop_file.exists():stop_file.unlink()
    with (out/'RUN_SESSIONS.jsonl').open('a') as sessions:
        sessions.write(json.dumps({'start_move':n,'additional_seconds':seconds,
            'code_sha256':sha_file(__file__),'reseed':reseed,'neutral_limit':neutral_limit,
            'growth_ceiling':growth,'mix_steps':mix_steps,'policy':'grow_in_burst_then_mix_then_reduce'})+'\n')
    def stop_after_move(signum,frame):
        nonlocal stopping
        stopping=True
    old_sigint=signal.signal(signal.SIGINT,stop_after_move)
    old_sigterm=signal.signal(signal.SIGTERM,stop_after_move)
    def checkpoint(reason):
        save_tri(out/'current.rga',tri)
        report={'status':reason,'moves':n,'best_moves':best_move,'current_pentachora':tri.size(),
                'best_pentachora':best_size,'current_tight_sha256':current,'best_tight_sha256':best_hash,
                'current_rga_sha256':sha_file(out/'current.rga'),'best_rga_sha256':sha_file(out/'best.rga'),
                'source_sha256':source_hash,'counts':dict(counts),'max_pentachora':max_size,
                'elapsed_seconds':previous_elapsed+time.monotonic()-started,
                'session_seconds':time.monotonic()-started,
                'rng_state':rng.getstate(),'neutral_since_reduction':neutral,
                'growth_remaining':growth_remaining,'mix_remaining':mix_remaining,
                'growth_cycles':growth_cycles,
                'effective_policy':{'neutral_limit':neutral_limit,'growth_ceiling':growth,'mix_steps':mix_steps},
                'standardness_asserted':False,
                'closed_completion_asserted':False,'native_simplify_used':False}
        save_json(out/'CHECKPOINT.json',report);return report
    checkpoint('RUNNING')
    def commit(log,op):
        nonlocal n,current,max_size,best_size,best_hash,best_move
        before=current;perform(tri,op);current=state_hash(tri)
        row={'index':n,'before_sha256':before,'operation':op,'after_sha256':current,
             'pentachora':tri.size(),'elapsed_seconds':previous_elapsed+time.monotonic()-started}
        log.write(json.dumps(row,separators=(',',':'))+'\n');log.flush()
        n+=1;counts[move_key(op)]+=1;max_size=max(max_size,tri.size())
        if tri.size()<best_size:
            best_size=tri.size();best_hash=current;best_move=n;save_tri(out/'best.rga',tri)
        checkpoint('RUNNING')
    with (out/'moves.jsonl').open('a' if resume else 'w') as log:
        try:
            if not resume and ideal_first_color is not None:
                require(0<=ideal_first_color<5,'invalid ideal local color')
                require(tri.pentachoron(0).vertex(ideal_first_color).isIdeal(),
                        'chosen first corner is not an ideal vertex')
                perm=list(range(5));perm[0],perm[ideal_first_color]=perm[ideal_first_color],perm[0]
                commit(log,['R',[[i,perm] for i in range(tri.size())]])
                # Pay the cold skeleton recognition once, on the simple raw
                # source, then preserve its genuine ideal count through moves.
                require(tri.pentachoron(0).vertex(0).isIdeal(),'relabelled corner is not ideal')
            if snap_to_one:
                require(sum(v.isIdeal() for v in tri.vertices())==1,'snap phase requires exactly one ideal vertex')
                while tri.countVertices()>1 and not stopping and not stop_file.exists():
                    candidates=[e for e in tri.edges()
                                if e.vertex(0)!=e.vertex(1)
                                and (e.vertex(0).isIdeal() or e.vertex(1).isIdeal())
                                and tri.hasSnapEdge(e)]
                    require(candidates,'cannot join remaining finite vertices to ideal vertex')
                    # Prefer absorbing global vertex0 first, so the ideal
                    # vertex is first in subsequent skeleton enumeration.
                    candidates.sort(key=lambda e:(min(e.vertex(0).index(),e.vertex(1).index()),e.index()))
                    commit(log,['S',candidates[0].index()])
            while not stopping and time.monotonic()-started<seconds:
                if stop_file.exists():mode='STOP_FILE_REQUEST';break
                options=[];scan=time.monotonic()
                if growth_remaining:
                    options=eligible(tri,'P',3,ideal_first_color is not None)
                    if options:
                        growth_remaining-=1
                        commit(log,rng.choice(options))
                        continue
                    growth_remaining=0
                if mix_remaining:
                    options=eligible(tri,'P',2,ideal_first_color is not None)
                    if options:
                        mix_remaining-=1
                        commit(log,rng.choice(options))
                        continue
                    mix_remaining=0
                for kind,k in [('Z',2),('Z',1),('P',1),('P',0),('Z',0)]:
                    options=eligible(tri,kind,k,ideal_first_color is not None)
                    if options:break
                reducing=bool(options)
                if not options:
                    if neutral>=neutral_limit:
                        if growth and tri.size()+2<=best_size+growth:
                            growth_remaining=max(1,(best_size+growth-tri.size())//2)
                            mix_remaining=mix_steps;growth_cycles+=1
                            neutral=0
                            continue
                        if not options:
                            mode='LOCAL_PLATEAU';break
                    else:
                        options=eligible(tri,'P',2,ideal_first_color is not None);neutral+=1
                if not options:mode='NO_ALLOWED_MOVE';break
                op=rng.choice(options)
                if reducing:neutral=0
                commit(log,op)
                if time.monotonic()-last_update>=10:
                    print(json.dumps({'moves':n,'current':tri.size(),'best':best_size,
                                      'scan_and_move_seconds':time.monotonic()-scan,'seconds':time.monotonic()-started}),flush=True)
                    last_update=time.monotonic()
            else:mode='INTERRUPTED_AFTER_LAST_COMPLETE_MOVE' if stopping else 'WALL_BUDGET'
        except KeyboardInterrupt:
            # SIGINT normally only sets our flag.  If an embedding injects a
            # Python KeyboardInterrupt directly, do not export potentially
            # half-committed memory.  Resume replays the atomic JSON prefix.
            print('Interrupted outside the normal signal path; preserving the last committed checkpoint for guarded replay.',flush=True)
            signal.signal(signal.SIGINT,old_sigint);signal.signal(signal.SIGTERM,old_sigterm)
            return
    publish_closed_log(out/'moves.jsonl',n,current)
    result=checkpoint(mode)
    signal.signal(signal.SIGINT,old_sigint);signal.signal(signal.SIGTERM,old_sigterm)
    print(json.dumps({k:v for k,v in result.items() if k!='rng_state'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('search');s.add_argument('source',type=Path);s.add_argument('--out',type=Path,required=True)
    s.add_argument('--seconds',type=float,default=600);s.add_argument('--seed',type=int,default=0)
    s.add_argument('--neutral-limit',type=int,default=100);s.add_argument('--growth',type=int,default=0)
    s.add_argument('--mix-steps',type=int,default=200)
    s.add_argument('--ideal-first-color',type=int)
    s.add_argument('--snap-to-one-ideal',action='store_true')
    u=sub.add_parser('resume');u.add_argument('directory',type=Path);u.add_argument('--source',type=Path)
    u.add_argument('--seconds',type=float,default=600)
    u.add_argument('--growth',type=int);u.add_argument('--neutral-limit',type=int)
    u.add_argument('--mix-steps',type=int);u.add_argument('--reseed',type=int)
    r=sub.add_parser('replay');r.add_argument('directory',type=Path);r.add_argument('--source',type=Path)
    r.add_argument('--best',action='store_true');r.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command=='search':search(a.source,a.out,a.seconds,a.seed,a.neutral_limit,a.growth,a.ideal_first_color,
                                 snap_to_one=a.snap_to_one_ideal,mix_steps=a.mix_steps)
    elif a.command=='resume':
        md=json.loads((a.directory/'SOURCE.json').read_text());policy=dict(md['policy'])
        sessions=a.directory/'RUN_SESSIONS.jsonl'
        if sessions.exists():
            last=json.loads(sessions.read_text().splitlines()[-1])
            for key in ('neutral_limit','growth_ceiling','mix_steps'):
                if key in last:policy[key]=last[key]
        cp=json.loads((a.directory/'CHECKPOINT.json').read_text())
        policy.update(cp.get('effective_policy',{}))
        search(a.source or Path(md['source_file']),a.directory,a.seconds,md['seed'],
               policy['neutral_limit'] if a.neutral_limit is None else a.neutral_limit,
               policy['growth_ceiling'] if a.growth is None else a.growth,
               policy.get('ideal_first_color'),resume=True,
               mix_steps=policy.get('mix_steps',200) if a.mix_steps is None else a.mix_steps,
               reseed=a.reseed)
    else:
        source_data=json.loads((a.directory/'SOURCE.json').read_text());c=json.loads((a.directory/'CHECKPOINT.json').read_text())
        source=a.source or Path(source_data['source_file']);require(sha_file(source)==source_data['source_sha256'],'source binding differs')
        result,tri=replay(source,a.directory/'moves.jsonl',c['best_moves'] if a.best else c['moves'],
                          c['best_tight_sha256'] if a.best else c['current_tight_sha256'])
        save_json(a.out,result);print(json.dumps(result),flush=True)

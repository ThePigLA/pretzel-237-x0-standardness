#!/usr/bin/env python3
"""Bounded bidirectional pure3D Pachner BFS between small boundaries."""
import argparse,collections,json,time
from pathlib import Path
import regina
from boundary_pachner_fast import face,mapped_face_index,encoding,digest,replay,require,save

def neighbours(sig,ceiling):
    t=regina.Triangulation3.fromIsoSig(sig)
    for d in range(4):
        if t.size()+2*d-3>ceiling:continue
        for f in (t.simplices() if d==3 else t.faces(d)):
            if not t.hasPachner(f):continue
            n=regina.Triangulation3(t);i=f.index();require(n.pachner(face(n,d,i)),'guard changed')
            yield n.isoSig(),d,i

def route(parents,end):
    moves=[]
    while parents[end] is not None:
        before,d,i=parents[end];moves.append({'before':before,'dimension':d,'face':i,'after':end});end=before
    return list(reversed(moves))

def inverse(move):
    d=3-move['dimension'];t=regina.Triangulation3.fromIsoSig(move['after'])
    for f in (t.simplices() if d==3 else t.faces(d)):
        if not t.hasPachner(f):continue
        n=regina.Triangulation3(t);i=f.index();n.pachner(face(n,d,i))
        if n.isoSig()==move['before']:
            return {'before':move['after'],'dimension':d,'face':i,'after':move['before']}
    raise ValueError('inverse Pachner edge not found')

def labelled(path):
    live=regina.Triangulation3.fromIsoSig(path['start']);start=encoding(live);moves=[]
    for step,m in enumerate(path['moves']):
        original=regina.Triangulation3.fromIsoSig(m['before']);iso=original.isIsomorphicTo(live)
        require(iso is not None,'path discontinuity')
        i=mapped_face_index(original,live,iso,m['dimension'],m['face'])
        before=digest(encoding(live));n=live.size();require(live.pachner(face(live,m['dimension'],i)),'transported move fails')
        require(live.isoSig()==m['after'],'wrong transported endpoint')
        moves.append({'dimension':m['dimension'],'face':i,'before_sha256':before,
                      'after_sha256':digest(encoding(live)),'source_tetrahedra':n,'target_tetrahedra':live.size()})
    data={'schema':'motoo.boundary.labelled-pachner.v1','status':'BFS_BRIDGE_REPLAYED',
          'start_tight':start,'end_tight':encoding(live),'start_tetrahedra':regina.Triangulation3.fromIsoSig(path['start']).size(),
          'end_tetrahedra':live.size(),'moves':moves,'all_moves_guarded':True,'native_simplify_used':False}
    _,data['replay']=replay(data);return data

def search(start,goal,out,ceiling,seconds):
    out=Path(out);require(not out.exists(),'output exists');out.mkdir(parents=True)
    parents=[{start:None},{goal:None}];queues=[collections.deque([start]),collections.deque([goal])]
    now=time.monotonic();expanded=[0,0];meet=start if start==goal else None;lastprint=now
    while meet is None and queues[0] and queues[1] and time.monotonic()-now<seconds:
        side=0 if len(queues[0])<=len(queues[1]) else 1
        current=queues[side].popleft();expanded[side]+=1
        for after,d,i in neighbours(current,ceiling):
            if after in parents[side]:continue
            parents[side][after]=(current,d,i);queues[side].append(after)
            if after in parents[1-side]:meet=after;break
        if time.monotonic()-lastprint>5:
            print(json.dumps({'states':[len(p) for p in parents],'expanded':expanded,'seconds':time.monotonic()-now}),flush=True);lastprint=time.monotonic()
    report={'status':'FOUND' if meet is not None else 'NOT_FOUND_WITHIN_RESOURCE_BOUND',
            'ceiling_tetrahedra':ceiling,'seconds':time.monotonic()-now,'states':[len(p) for p in parents],
            'expanded':expanded,'start':start,'end':goal}
    if meet is not None:
        moves=route(parents[0],meet)+[inverse(m) for m in reversed(route(parents[1],meet))]
        canonical={'schema':'motoo.boundary.pachner-path.v1','status':'BFS_FOUND','start':start,'end':goal,'moves':moves}
        save(out/'canonical_path.json',canonical);save(out/'labelled_path.json',labelled(canonical));report['moves']=len(moves)
    save(out/'BFS_RESULT.json',report);print(json.dumps(report),flush=True);return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--start',required=True);p.add_argument('--goal',required=True)
    p.add_argument('--out',required=True);p.add_argument('--ceiling',type=int,default=10);p.add_argument('--seconds',type=float,default=60)
    a=p.parse_args();search(a.start,a.goal,a.out,a.ceiling,a.seconds)

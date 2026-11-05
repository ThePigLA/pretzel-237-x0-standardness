#!/usr/bin/env python3
"""PurePachner3D reduction using forced two-move4-4 excursions.

A legal2-3 on a face through a degree-four edge reduces that old edge's
degree to three. If it then admits3-2, execute that old-edge deletion,
instead of immediately undoing the newly introduced edge. Both steps
are independently guarded Pachner moves; no native4-4 operation is used.
"""
import argparse,json,random,time
from pathlib import Path
import regina
from boundary_pachner_fast import digest,encoding,face,replay,require,save

def pillow_macro(t,rng):
    """Search a four-Pachner replacement for an admissible edge2-0.

    Native2-0 is only a target oracle here. An accepted replacement is
    entirely composed of guarded Pachner moves and must end at a
    triangulation isomorphic to that independently computed target.
    """
    edges=[e for e in t.edges() if e.degree()==2 and t.has20(e)];rng.shuffle(edges)
    for edge in edges[:12]:
        target=regina.Triangulation3(t);require(target.move20(target.edge(edge.index())),'2-0 target guard changed')
        embeddings=list(edge.embeddings());options=[]
        for e in embeddings:
            for f in (e.vertices()[0],e.vertices()[1]):options.append(e.simplex().triangle(f).index())
        rng.shuffle(options)
        for index in options:
            tri=t.triangle(index)
            if not t.hasPachner(tri):continue
            removed={e.simplex().index() for e in tri.embeddings()}
            survivor=next((e for e in embeddings if e.simplex().index() not in removed),None)
            if survivor is None:continue
            vertices={survivor.vertices()[0],survivor.vertices()[1]};ti=survivor.simplex().index()
            le=next(i for i in range(6) if {survivor.simplex().edgeMapping(i)[0],survivor.simplex().edgeMapping(i)[1]}==vertices)
            a=regina.Triangulation3(t);kept=a.simplex(ti);a.pachner(a.triangle(index));oldedge=kept.edge(le)
            if not a.hasPachner(oldedge):continue
            ei=oldedge.index();a.pachner(oldedge)
            for e1 in a.edges():
                if e1.degree()!=3 or not a.hasPachner(e1):continue
                i1=e1.index();b=regina.Triangulation3(a);b.pachner(b.edge(i1))
                for e2 in b.edges():
                    if e2.degree()!=3 or not b.hasPachner(e2):continue
                    i2=e2.index();c=regina.Triangulation3(b);c.pachner(c.edge(i2))
                    if c.isIsomorphicTo(target) is not None:return [(2,index),(1,ei),(1,i1),(1,i2)]
    return None

def macro(t,rng):
    edges=[e for e in t.edges() if e.degree()==4];rng.shuffle(edges)
    for edge in edges[:40]:
        embeddings=list(edge.embeddings());triangles=set()
        for e in embeddings:
            vertices={e.vertices()[0],e.vertices()[1]}
            for f in set(range(4))-vertices:triangles.add(e.simplex().triangle(f).index())
        choices=list(triangles);rng.shuffle(choices)
        for index in choices:
            tri=t.triangle(index)
            if not t.hasPachner(tri):continue
            removed={e.simplex().index() for e in tri.embeddings()}
            survivor=next((e for e in embeddings if e.simplex().index() not in removed),None)
            if survivor is None:continue
            tet_index=survivor.simplex().index();vertices={survivor.vertices()[0],survivor.vertices()[1]}
            local_edge=next(i for i in range(6) if {survivor.simplex().edgeMapping(i)[0],survivor.simplex().edgeMapping(i)[1]}==vertices)
            probe=regina.Triangulation3(t);kept=probe.simplex(tet_index)
            require(probe.pachner(probe.triangle(index)),'macro first guard changed')
            second=kept.edge(local_edge)
            if probe.hasPachner(second):return [(2,index),(1,second.index())]
    return None

def run(source,out,seconds,seed,goal_size=7):
    out=Path(out);require(not out.exists(),'output exists');out.mkdir(parents=True)
    t=regina.Triangulation3(regina.open(str(source)));start=encoding(t);n0=t.size();rng=random.Random(seed)
    moves=[];bestmoves=[];best=start;bestscore=(t.countVertices(),t.size());now=lastsave=time.monotonic();steps=0;macros=0
    def checkpoint(status):
        endpoint=regina.Triangulation3.tightDecoding(best)
        data={'schema':'motoo.boundary.labelled-pachner.v1','status':status,'start_tight':start,'end_tight':best,
              'start_tetrahedra':n0,'end_tetrahedra':endpoint.size(),'moves':bestmoves,'seed':seed,
              'all_moves_guarded':True,'native_simplify_used':False,'forced_two_move_macros':macros,
              'attempted_pachner_moves':steps,'seconds':time.monotonic()-now}
        _,data['replay']=replay(data);save(out/'best_labelled_path.json',data)
        regina.PacketOfTriangulation3(endpoint).save(str(out/'best_boundary.rga'));return data
    while time.monotonic()-now<seconds:
        if t.size()<=goal_size:
            best,bestmoves,bestscore=encoding(t),list(moves),(t.countVertices(),t.size())
            result=checkpoint('SMALL_BOUNDARY_FOUND');print(json.dumps({'status':result['status'],'moves':len(bestmoves),'tetrahedra':t.size(),'isosig':t.isoSig(),'seconds':result['seconds']}),flush=True);return result
        choices=[(0,v.index()) for v in t.vertices() if v.degree()==4 and t.hasPachner(v)]
        choices=choices or [(1,e.index()) for e in t.edges() if e.degree()==3 and t.hasPachner(e)]
        block=[rng.choice(choices)] if choices else None
        if block is None:
            block=pillow_macro(t,rng)
            if block:macros+=1
        if block is None and macros%100!=99:
            block=macro(t,rng)
            if block:macros+=1
        if block is None:
            for attempt in range(100):
                i=rng.randrange(t.countTriangles())
                if t.hasPachner(t.triangle(i)):block=[(2,i)];break
            if block is None:block=[(3,rng.randrange(t.size()))]
            macros+=1
        for d,i in block:
            before=digest(encoding(t));n=t.size();require(t.pachner(face(t,d,i)),'selectedPachnerguard failed')
            moves.append({'dimension':d,'face':i,'before_sha256':before,'after_sha256':digest(encoding(t)),
                          'source_tetrahedra':n,'target_tetrahedra':t.size()});steps+=1
        score=(t.countVertices(),t.size())
        if score<bestscore:
            bestscore,best,bestmoves=score,encoding(t),list(moves)
            if t.size()%20==0 or t.size()<20:
                print(json.dumps({'best':score,'moves':len(moves),'seconds':time.monotonic()-now}),flush=True)
        if time.monotonic()-lastsave>20:
            checkpoint('RUNNING');lastsave=time.monotonic()
    result=checkpoint('TIME_BUDGET_EXHAUSTED');print(json.dumps({'status':result['status'],'best':bestscore,'moves':len(bestmoves),'seconds':result['seconds']}),flush=True);return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    p.add_argument('--seconds',type=float,default=60);p.add_argument('--seed',type=int,default=0);p.add_argument('--goal-size',type=int,default=7)
    a=p.parse_args();run(a.source,a.out,a.seconds,a.seed,a.goal_size)

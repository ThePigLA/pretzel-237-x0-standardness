#!/usr/bin/env python3
"""Label-preserving, logged pure 3D Pachner search for upper-handle gluing.

Exploration avoids computing an isomorphism signature at every step.
Every move has before/after tight-encoding digests and is guarded by
Regina. The winning path can be converted to the existing canonical
boundary_pachner_path schema, so the established 4D layering/checker
continues to apply unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import time

import regina


def require(ok, message):
    if not ok:
        raise ValueError(message)


def face(t, d, i):
    return t.simplex(i) if d == 3 else t.face(d, i)


def encoding(t):
    return t.tightEncoding()


def digest(s):
    return hashlib.sha256(s.encode()).hexdigest()


def save(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def replay(data):
    t = regina.Triangulation3.tightDecoding(data["start_tight"])
    h = digest(encoding(t))
    for step, m in enumerate(data["moves"]):
        require(h == m["before_sha256"], f"source mismatch at {step}")
        require(t.pachner(face(t, m["dimension"], m["face"])), f"illegal move at {step}")
        h = digest(encoding(t))
        require(h == m["after_sha256"], f"target mismatch at {step}")
    require(encoding(t) == data["end_tight"], "wrong final labelled triangulation")
    return t, {"status":"PASS_EXPLICIT_3D_PACHNER_REPLAY", "moves":len(data["moves"]),
               "source_tetrahedra":data["start_tetrahedra"], "target_tetrahedra":t.size()}


def mapped_face_index(original, canonical, relabel, d, i):
    if d == 3:
        return relabel.simpImage(i)
    e = face(original, d, i).embedding(0)
    old_tet = e.simplex().index()
    tet = canonical.simplex(relabel.simpImage(old_tet))
    perm = relabel.facetPerm(old_tet)
    vertices = {perm[e.vertices()[j]] for j in range(d+1)}
    if d == 0:
        return tet.vertex(next(iter(vertices))).index()
    if d == 2:
        return tet.triangle(next(iter(set(range(4))-vertices))).index()
    for edge in range(6):
        mapping = tet.edgeMapping(edge)
        if {mapping[0],mapping[1]} == vertices:
            return tet.edge(edge).index()
    raise ValueError("edge relabelling could not be recovered")


def canonicalise(data, output):
    """Convert a finite retained path, checking every canonical move."""
    t = regina.Triangulation3.tightDecoding(data["start_tight"])
    before, relabel = t.isoSigDetail()
    moves = []
    start = before
    then = time.monotonic()
    for step, move in enumerate(data["moves"]):
        require(digest(encoding(t)) == move["before_sha256"], f"bad source at {step}")
        canonical = regina.Triangulation3.fromIsoSig(before)
        d, i = move["dimension"], move["face"]
        ci = mapped_face_index(t, canonical, relabel, d, i)
        require(canonical.pachner(face(canonical,d,ci)), f"canonical guard failed at {step}")
        require(t.pachner(face(t,d,i)), f"labelled guard failed at {step}")
        require(digest(encoding(t)) == move["after_sha256"], f"labelled endpoint differs at {step}")
        after, relabel = t.isoSigDetail()
        require(canonical.isoSig() == after, f"canonical endpoint differs at {step}")
        moves.append({"before":before,"dimension":d,"face":ci,"after":after})
        before = after
        if (step+1)%100 == 0:
            print(json.dumps({"canonicalised_moves":step+1,"seconds":time.monotonic()-then}),flush=True)
    result = {"schema":"motoo.boundary.pachner-path.v1", "status":"CANONICALISED_AND_REPLAYED",
              "start":start,"end":before,"moves":moves,
              "source_labelled_path_sha256":digest(json.dumps(data,sort_keys=True)),
              "seconds":time.monotonic()-then}
    save(output,result)
    return result


def search(source_path, output, seconds, seed, goals, growth=12, loop_limit=100, loop_period=20):
    output = Path(output)
    require(not output.exists(), "output exists")
    output.mkdir(parents=True)
    source = regina.open(str(source_path))
    require(isinstance(source,regina.Triangulation3), "input is not a 3D triangulation")
    t = regina.Triangulation3(source)
    require(t.isValid() and t.isClosed() and t.isConnected(), "source boundary is not closed/valid/connected")
    initial = encoding(t)
    best_encoding = initial
    best_score = (t.countVertices(),t.size())
    moves, best_moves = [], []
    rng = random.Random(seed)
    started = last_save = time.monotonic()
    steps = restarts = 0
    high = t.size()
    mode = "reduce"
    growth_target = t.size()
    before = digest(initial)
    counters = [0]*4
    short_positions = {}
    erased_moves = 0

    def checkpoint(status):
        best = regina.Triangulation3.tightDecoding(best_encoding)
        result = {"schema":"motoo.boundary.labelled-pachner.v1", "status":status,
                  "source_file":str(source_path), "source_sha256":hashlib.sha256(Path(source_path).read_bytes()).hexdigest(),
                  "seed":seed,"start_tight":initial,"end_tight":best_encoding,
                  "start_tetrahedra":source.size(),"end_tetrahedra":best.size(),
                  "best_score_vertices_tetrahedra":list(best_score),"moves":best_moves,
                  "attempted_moves":steps,"restarts":restarts,"max_tetrahedra":high,
                  "move_counts_by_dimension":counters,"seconds":time.monotonic()-started,
                  "isomorphism_loop_limit":loop_limit,"erased_moves":erased_moves,
                  "isomorphism_loop_period":loop_period,
                  "native_simplify_used":False,"all_moves_guarded":True}
        save(output/"best_labelled_path.json",result)
        regina.PacketOfTriangulation3(best).save(str(output/"best_boundary.rga"))
        return result

    while time.monotonic()-started < seconds:
        if t.size() <= 10 and t.isoSig() in goals:
            best_encoding, best_moves, best_score = encoding(t),list(moves),(t.countVertices(),t.size())
            result=checkpoint("FOUND")
            _,result["replay"] = replay(result)
            save(output/"SUCCESS_LABELLED.json",result)
            print(json.dumps({"status":"FOUND","steps":steps,"moves":len(best_moves),"tetrahedra":t.size()}),flush=True)
            return result
        vertex_moves = [(0,v.index()) for v in t.vertices() if v.degree()==4 and t.hasPachner(v)]
        reductions = vertex_moves or [(1,e.index()) for e in t.edges() if e.degree()==3 and t.hasPachner(e)]
        if mode == "reduce" and reductions:
            choices = reductions
        else:
            if mode == "reduce":
                if t.size() >= best_score[1]+growth or t.countVertices() > best_score[0]+1:
                    t = regina.Triangulation3.tightDecoding(best_encoding)
                    moves = list(best_moves)
                    before = digest(best_encoding)
                    short_positions = {}
                    restarts += 1
                    continue
                mode = "grow"
                growth_target = min(best_score[1]+growth,t.size()+rng.choice((2,4,6,8)))
            if t.size() >= growth_target:
                mode = "reduce"
                choices = reductions
            else:
                choices = []
            if not choices:
                # Sampling avoids a complete legality scan over thousands
                # of triangles each time a single growth move is needed.
                choices=[]
                for trial in range(min(64,t.countTriangles())):
                    triangle=t.triangle(rng.randrange(t.countTriangles()))
                    if t.hasPachner(triangle):
                        choices.append((2,triangle.index()))
                        break
                if not choices:
                    choices=[(2,f.index()) for f in t.triangles() if t.hasPachner(f)]
                if not choices or (t.size()+3 <= growth_target and rng.random()<0.003):
                    choices=[(3,rng.randrange(t.size()))]
        d,i=rng.choice(choices)
        n=t.size()
        require(t.pachner(face(t,d,i)),"selected move failed")
        tight=encoding(t);after=digest(tight)
        moves.append({"dimension":d,"face":i,"before_sha256":before,"after_sha256":after,
                      "source_tetrahedra":n,"target_tetrahedra":t.size()})
        before=after;steps+=1;counters[d]+=1;high=max(high,t.size())
        if t.size() <= loop_limit and (t.size()<=20 or steps%loop_period==0):
            sig=t.isoSig()
            previous=short_positions.get(sig)
            if previous is not None:
                where,old_tight=previous
                erased_moves += len(moves)-where
                moves=moves[:where]
                short_positions={k:v for k,v in short_positions.items() if v[0]<=where}
                # Restart at the EXACT earlier labels, so retained move
                # indices need no unlogged relabelling bridge.
                t=regina.Triangulation3.tightDecoding(old_tight)
                tight=old_tight;before=digest(tight)
            else:
                short_positions[sig]=(len(moves),tight)
        score=(t.countVertices(),t.size())
        if score<best_score:
            best_score=score;best_encoding=tight;best_moves=list(moves);mode="reduce"
            if steps%10==0 or t.size()<100:
                print(json.dumps({"best":best_score,"steps":steps,"retained":len(moves),"seconds":time.monotonic()-started}),flush=True)
        if time.monotonic()-last_save>20:
            checkpoint("RUNNING");last_save=time.monotonic()
    result=checkpoint("TIME_BUDGET_EXHAUSTED")
    _,result["replay"] = replay(result)
    save(output/"best_labelled_path.json",result)
    print(json.dumps({"status":result["status"],"best":best_score,"steps":steps,"retained":len(best_moves),"replay":result["replay"]}),flush=True)
    return result


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source")
    p.add_argument("--out",required=True)
    p.add_argument("--seconds",type=float,default=120)
    p.add_argument("--seed",type=int,default=0)
    p.add_argument("--growth",type=int,default=12)
    p.add_argument("--loop-limit",type=int,default=100)
    p.add_argument("--loop-period",type=int,default=20)
    p.add_argument("--goal",action="append",default=[])
    p.add_argument("--canonicalise")
    a=p.parse_args()
    if a.canonicalise:
        canonicalise(json.loads(Path(a.canonicalise).read_text()),a.out)
    else:
        require(a.source,"supply --source")
        search(a.source,a.out,a.seconds,a.seed,set(a.goal or ["gLALQccbefefdpcghfc"]),a.growth,a.loop_limit,a.loop_period)

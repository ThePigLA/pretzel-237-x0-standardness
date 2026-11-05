"""Extend a labelled Pachner path by breadth-first search (on labelled representatives) to any goal isoSig.

usage: bridge_bfs.py path.json out.json goal1,goal2 [ceiling]
The input path's end_tight is the start; moves found are appended in the same labelled format.
"""
import sys, json, hashlib, collections, time, regina


def dig(s):
    return hashlib.sha256(s.encode()).hexdigest()


def face(t, d, i):
    return t.simplex(i) if d == 3 else t.face(d, i)


def count(t, d):
    return t.size() if d == 3 else t.countFaces(d)


src, out, goals = sys.argv[1], sys.argv[2], set(sys.argv[3].split(','))
ceiling = int(sys.argv[4]) if len(sys.argv) > 4 else 10
data = json.load(open(src))
start = regina.Triangulation3.tightDecoding(data['end_tight'])
t0 = time.time()
parent = {start.isoSig(): None}; rep = {start.isoSig(): start}
q = collections.deque([start.isoSig()]); found = start.isoSig() if start.isoSig() in goals else None
while q and not found:
    s = q.popleft(); cur = rep[s]
    for d in (1, 0, 2, 3):
        for i in range(count(cur, d)):
            if not cur.hasPachner(face(cur, d, i)):
                continue
            nxt = regina.Triangulation3(cur); nxt.pachner(face(nxt, d, i))
            if nxt.size() > ceiling:
                continue
            ns = nxt.isoSig()
            if ns in parent:
                continue
            parent[ns] = (s, d, i); rep[ns] = nxt; q.append(ns)
            if ns in goals:
                found = ns; break
        if found:
            break
print('states', len(parent), 'found', found, '%.1fs' % (time.time()-t0), flush=True)
if not found:
    sys.exit(1)
chain = []; s = found
while parent[s] is not None:
    p, d, i = parent[s]; chain.append((d, i)); s = p
t = regina.Triangulation3(start)
moves = list(data['moves'])
for d, i in reversed(chain):
    b = dig(t.tightEncoding()); n = t.size()
    assert t.pachner(face(t, d, i))
    moves.append(dict(dimension=d, face=i, before_sha256=b, after_sha256=dig(t.tightEncoding()),
                      source_tetrahedra=n, target_tetrahedra=t.size()))
assert t.isoSig() == found
out_data = dict(data); out_data.update(moves=moves, end_tight=t.tightEncoding(), end_tetrahedra=t.size(),
                                      bridge=dict(goal=found, moves=len(chain), states=len(parent)))
json.dump(out_data, open(out, 'w'))
print('bridged with', len(chain), 'moves; total', len(moves), 'end', found)

"""Labelled pure-Pachner path on a closed 3-manifold boundary triangulation, ending isomorphic to a goal.

Writes the path in the format used by boundary_layer_fast.layer_labelled:
  {"start_tight", "end_tight", "moves": [{"dimension", "face", "before_sha256", "after_sha256"}]}
Every move is applied with Regina's checked Triangulation3.pachner on the *labelled* current state, so the
path replays verbatim.  Phase 1: randomised descent (3-2 and 4-1 moves, with 2-3 kicks).  Phase 2: breadth-first
search over labelled representatives (keyed by isoSig) until a state isomorphic to a goal appears.
usage: boundary_path.py boundary.rga out.json goal_isosig[,goal2...] [seconds] [seed]
"""
import sys, json, hashlib, random, time, collections, regina


def enc(t):
    return t.tightEncoding()


def dig(s):
    return hashlib.sha256(s.encode()).hexdigest()


def face(t, d, i):
    return t.simplex(i) if d == 3 else t.face(d, i)


def count(t, d):
    return t.size() if d == 3 else t.countFaces(d)


def legal(t, d):
    return [i for i in range(count(t, d)) if t.hasPachner(face(t, d, i))]


src, out, goals = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
secs = float(sys.argv[4]) if len(sys.argv) > 4 else 600
rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 0)
goal_tris = [regina.Triangulation3.fromIsoSig(g) for g in goals]
goal_sigs = {g.isoSig() for g in goal_tris}
goal_size = min(g.size() for g in goal_tris)
t = regina.Triangulation3(regina.open(src))
start_tight = enc(t)
moves = []
t0 = time.time()


def do(t, d, i):
    b = dig(enc(t))
    assert t.pachner(face(t, d, i)), 'illegal move'
    moves.append(dict(dimension=d, face=i, before_sha256=b, after_sha256=dig(enc(t))))


print('start', t.size(), 'tets', t.countVertices(), 'vertices', flush=True)
# phase 1: descent with kicks, on the labelled object (moves recorded as we go)
best = t.size(); stall = 0
while time.time()-t0 < secs and t.size() > goal_size+2:
    red = [(1, i) for i in legal(t, 1)] + [(0, i) for i in legal(t, 0)]
    if red:
        d, i = rng.choice(red); do(t, d, i); stall = 0
    else:
        stall += 1
        ups = legal(t, 2)
        # prefer 2-3 moves that create a degree-3 edge or a degree-4 vertex soon: just random kick
        d, i = 2, rng.choice(ups); do(t, d, i)
    if t.size() < best:
        best = t.size()
        if best % 50 == 0 or best < 40:
            print(f'  {time.time()-t0:.0f}s size {t.size()} vertices {t.countVertices()} moves {len(moves)}', flush=True)
print('phase 1 end: size', t.size(), 'vertices', t.countVertices(), 'moves', len(moves), flush=True)
# compress 1-vertex: also allow 4-1 moves already included.  phase 2: BFS from current labelled state
ceiling = max(t.size(), goal_size)+3
start = regina.Triangulation3(t)
parent = {start.isoSig(): None}
rep = {start.isoSig(): start}
q = collections.deque([start.isoSig()])
found = None
while q and time.time()-t0 < secs*3:
    s = q.popleft(); cur = rep[s]
    if s in goal_sigs:
        found = s; break
    for d in (1, 2, 0, 3):
        for i in range(count(cur, d)):
            f = face(cur, d, i)
            if not cur.hasPachner(f):
                continue
            nxt = regina.Triangulation3(cur)
            nxt.pachner(face(nxt, d, i))
            if nxt.size() > ceiling:
                continue
            ns = nxt.isoSig()
            if ns in parent:
                continue
            parent[ns] = (s, d, i); rep[ns] = nxt; q.append(ns)
            if ns in goal_sigs:
                found = ns; q.clear(); break
        if found:
            break
print('phase 2: states', len(parent), 'found', found is not None, f'{time.time()-t0:.0f}s', flush=True)
if not found:
    json.dump(dict(status='NO_PATH', start_tight=start_tight, moves=moves, end_tight=enc(t)), open(out, 'w'))
    sys.exit(1)
chain = []
s = found
while parent[s] is not None:
    p, d, i = parent[s]; chain.append((d, i)); s = p
for d, i in reversed(chain):
    do(t, d, i)
assert t.isoSig() in goal_sigs
json.dump(dict(schema='motoo.boundary.labelled-pachner-path.v1', status='PASS', start_tight=start_tight,
               end_tight=enc(t), end_isosig=t.isoSig(), moves=moves), open(out, 'w'))
print('PATH', len(moves), 'moves; end', t.isoSig(), flush=True)

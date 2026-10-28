"""Replay a (possibly partial) simplicial move certificate with check.py's rules and export the resulting complex
as Regina gluing data.  The export's facet order and vertex order are the sorted final facets, so the Regina
triangulation is literally that simplicial complex (checked back by matching every gluing to facet labels)."""
import gzip, json, sys, hashlib
from itertools import combinations
src, cert, out = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.loads(gzip.decompress(open(src, 'rb').read()))
def _load_moves(cert):
    if cert.endswith('.jsonl') or cert.endswith('.jsonl.gz'):
        op = gzip.open if cert.endswith('.gz') else open
        with op(cert, 'rt') as fh: return {'moves': [json.loads(l) for l in fh if l.strip()]}
    return json.loads(gzip.decompress(open(cert, 'rb').read())) if cert.endswith('.gz') else json.load(open(cert))
C = _load_moves(cert)
facets = set(frozenset(f) for f in d['facets']); vstar = {}
for f in facets:
    for v in f: vstar.setdefault(v, set()).add(f)
def having(face):
    face = list(face); s = set(vstar.get(face[0], ()))
    for v in face[1:]: s &= vstar.get(v, set())
    return s
def put(f):
    assert f not in facets; facets.add(f)
    for v in f: vstar.setdefault(v, set()).add(f)
def take(f):
    facets.remove(f)
    for v in f:
        vstar[v].discard(f)
        if not vstar[v]: del vstar[v]
def lk(fs, drop):
    o = set()
    for f in fs:
        r = sorted(f - drop)
        for k in range(1, len(r) + 1): o.update(frozenset(c) for c in combinations(r, k))
    return o
for n, m in enumerate(C['moves']):
    if m[0] == 'c':
        u, v = m[1], m[2]; Su, Sv = set(vstar[u]), set(vstar[v]); Suv = Su & Sv
        assert Suv and lk(Su, {u}) & lk(Sv, {v}) == lk(Suv, {u, v}), f'move {n}'
        for f in Suv: take(f)
        for f in list(vstar.get(u, ())): take(f); put(frozenset(v if x == u else x for x in f))
    else:
        S, T = frozenset(m[1]), frozenset(m[2]); A = having(S)
        assert not (S & T) and len(S) + len(T) == 6 and A == {S | (T - {x}) for x in T}, f'move {n}'
        assert (next(iter(T)) not in vstar) if len(T) == 1 else not having(T), f'move {n}'
        for f in A: take(f)
        for y in S: put((S - {y}) | T)
F = sorted(tuple(sorted(f)) for f in facets)
ridges = {}
for i, f in enumerate(F):
    for k in range(5):
        ridges.setdefault(tuple(x for j, x in enumerate(f) if j != k), []).append((i, k))
gl = []
for r, occ in ridges.items():
    assert len(occ) == 2, 'not closed'
    (i, a), (j, b) = occ
    if (i, a) > (j, b): (i, a), (j, b) = (j, b), (i, a)
    perm = [None] * 5
    for x in range(5):
        perm[x] = b if x == a else F[j].index(F[i][x])
    gl.append([i, a, j, perm])
# check back
for i, a, j, p in gl:
    for x in range(5):
        if x != a: assert F[j][p[x]] == F[i][x]
fsha = hashlib.sha256(json.dumps(F).encode()).hexdigest()
json.dump({'simplex_count': len(F), 'gluings': gl, 'facets': F, 'facets_sha256': fsha,
           'chain': {'source': src, 'certificate': cert, 'moves': len(C['moves'])}}, gzip.open(out, 'wt'))
print(f'replayed {len(C["moves"])} moves; exported {len(F)} facets, {len(vstar)} vertices; facets sha256 {fsha}')

"""My own check of the analyst's closed edge-contraction trace (master 103,744 -> 44,728 pentachora).
Reads the master facets directly from the marked master, applies each recorded contraction u -> v, and checks
the full link condition lk(u) & lk(v) == lk(uv) (all nonempty faces) before each step.  At the end the facet set
must equal the published endpoint, and the published Regina export must be exactly the endpoint's complex.
Under the link condition, a contraction in a closed combinatorial manifold preserves PL type (Nevo, Thm 1.4)."""
import gzip, json, sys, time, hashlib
from itertools import combinations
from collections import defaultdict
import os
R = os.path.join(sys.argv[1], '') if len(sys.argv) > 1 else '../inputs/'
t0 = time.time()
ld = lambda p: json.loads(gzip.decompress(open(R + p, 'rb').read()))
master = ld('marked_control.json.gz')
mf = [tuple(sorted(f)) for f in master['facets']]
src = ld('unmarked_contractions_v1/source_facets.json.gz')
assert [tuple(sorted(f)) for f in src['facets']] == mf, 'trace source differs from the marked master facets'
cert = ld('unmarked_contractions_v1/certificate.json.gz'); end = ld('unmarked_contractions_v1/endpoint.json.gz')
print('master facets', len(mf), 'steps', len(cert['steps']), flush=True)
# master is the complex of the Regina export that native admission checked
reg = ld('cover_regina.json.gz'); F = master['facets']
for i, k, j, p in reg['gluings']:
    for x in range(5):
        if x != k: assert F[j][p[x]] == F[i][x]
assert len(set(mf)) == len(mf)
facets = {i: f for i, f in enumerate(mf)}; star = defaultdict(set)
for i, f in facets.items():
    for v in f: star[v].add(i)
nxt = len(mf)
def faces_minus(fids, drop):
    out = set()
    for i in fids:
        rest = tuple(x for x in facets[i] if x not in drop)
        for k in range(1, len(rest) + 1):
            out.update(combinations(rest, k))
    return out
for n, st in enumerate(cert['steps']):
    u, v = st['u'], st['v']
    su, sv = star[u], star[v]; suv = su & sv
    assert suv, f'step {n}: uv is not an edge'
    assert len(facets) == st['before']
    lku = faces_minus(su, {u}); lkv = faces_minus(sv, {v}); lkuv = faces_minus(suv, {u, v})
    assert lku & lkv == lkuv, f'step {n}: link condition fails for {u},{v}'
    for i in list(suv):                      # delete facets containing the edge
        for w in facets[i]: star[w].discard(i)
        del facets[i]
    for i in list(star[u]):                  # u -> v in the remaining star of u
        f = tuple(sorted(v if w == u else w for w in facets[i]))
        assert v not in facets[i]
        for w in facets[i]: star[w].discard(i)
        del facets[i]; facets[nxt] = f
        for w in f: star[w].add(nxt)
        nxt += 1
    assert not star[u]; del star[u]
    assert len(facets) == st['after'], f'step {n}: count {len(facets)} != {st["after"]}'
fs = sorted(facets.values()); es = sorted(tuple(sorted(f)) for f in end['facets'])
assert len(set(fs)) == len(fs), 'endpoint has a repeated facet'
assert fs == es, 'endpoint facets differ'
er = ld('unmarked_contractions_v1/endpoint_regina.json.gz'); E = end['facets']; seen = set()
assert er['simplex_count'] == len(E)
for i, k, j, p in er['gluings']:
    for x in range(5):
        if x != k: assert E[j][p[x]] == E[i][x]
    seen.add((i, k)); seen.add((j, p[k]))
assert len(seen) == 5 * len(E)
print(f'PASS: {len(cert["steps"])} contractions, link condition each step; endpoint facets and Regina export match',
      f'({len(fs)} pentachora) [{time.time()-t0:.0f}s]', flush=True)

"""Certify that a 4-dimensional simplicial complex is a combinatorial 4-manifold: every vertex link is reduced to
the boundary of the 4-simplex by 3-dimensional bistellar flips only (no contractions, so no manifold premise is
needed: bistellar flips are PL homeomorphisms of arbitrary simplicial complexes).
Each link's flip sequence is produced by a small randomized reducer and then replayed by a separate verifier that
checks every flip's precondition from scratch.  Standard library only."""
import gzip, json, sys, time, random
from itertools import combinations

def reduce3(tets, rng, max_steps=200000):
    """bistellar reduction of a 3-dim complex; returns the flip list (sigma, tau) or None"""
    F = set(frozenset(t) for t in tets); star = {}
    for t in F:
        for v in t: star.setdefault(v, set()).add(t)
    log = []
    def having(face):
        face = list(face); s = set(star.get(face[0], ()))
        for v in face[1:]: s &= star.get(v, set())
        return s
    def flip(S, T):
        A = having(S)
        if A != {S | (T - {x}) for x in T}: return False
        if len(T) > 1 and having(T): return False
        for t in A:
            F.discard(t)
            for v in t:
                star[v].discard(t)
                if not star[v]: del star[v]
        for y in S:
            t = (S - {y}) | T; F.add(t)
            for v in t: star.setdefault(v, set()).add(t)
        log.append((sorted(S), sorted(T))); return True
    steps = 0
    while len(F) > 5 and steps < max_steps:
        steps += 1; done = False
        for v in [v for v, s in star.items() if len(s) == 4]:          # 4-1
            if v in star and len(star[v]) == 4:
                T = frozenset(x for t in star[v] for x in t) - {v}
                if len(T) == 4 and flip(frozenset([v]), T): done = True
        if done: continue
        edges = {}
        for t in F:
            for e in combinations(sorted(t), 2): edges[e] = edges.get(e, 0) + 1
        cands = [e for e, d in edges.items() if d == 3]; rng.shuffle(cands)
        for e in cands:                                                  # 3-2
            S = frozenset(e); A = having(S)
            T = frozenset(x for t in A for x in t) - S
            if len(T) == 3 and flip(S, T): done = True; break
        if done: continue
        tl = list(F)                                                     # 2-3 heating
        for _ in range(100):
            t = rng.choice(tl); S = frozenset(rng.sample(sorted(t), 3)); A = having(S)
            if len(A) != 2: continue
            T = frozenset(x for a in A for x in a) - S
            if len(T) == 2 and flip(S, T): break
    V = set(star)
    if len(V) == 5 and F == {frozenset(c) for c in combinations(sorted(V), 4)}: return log
    return None

def verify3(tets, log):
    """independent replay: every flip precondition, final boundary of the 4-simplex"""
    F = set(frozenset(t) for t in tets)
    assert len(F) == len(tets) and all(len(t) == 4 for t in F)
    tri = {}
    for t in F:
        for r in combinations(sorted(t), 3): tri[r] = tri.get(r, 0) + 1
    assert all(c == 2 for c in tri.values()), 'link is not a closed 3-pseudomanifold'
    for S, T in log:
        S, T = frozenset(S), frozenset(T)
        assert not (S & T) and len(S) + len(T) == 5
        A = {t for t in F if S <= t}
        assert A == {S | (T - {x}) for x in T}
        if len(T) > 1: assert not any(T <= t for t in F)
        else: assert not any(T <= t for t in F)
        F -= A; F |= {(S - {y}) | T for y in S}
    V = set().union(*F)
    return len(V) == 5 and F == {frozenset(c) for c in combinations(sorted(V), 4)}

def verify_saved(src, certpath):
    d = json.loads(gzip.decompress(open(src, 'rb').read()))
    c = json.loads(gzip.decompress(open(certpath, 'rb').read())) if certpath.endswith('.gz') else json.load(open(certpath))
    facets = [tuple(sorted(f)) for f in d['facets']]; assert len(set(facets)) == len(facets)
    links = {}
    for f in facets:
        for v in f: links.setdefault(v, []).append(tuple(x for x in f if x != v))
    assert len(links) == c['vertices'] == len(c['links'])
    t0 = time.time(); n = 0
    for v, L in links.items():
        assert verify3(L, c['links'][str(v)]), f'vertex {v}: link replay failed'
        n += len(c['links'][str(v)])
    print(f'PASS: all {len(links)} vertex links replay to the boundary of the 4-simplex ({n} flips) [{time.time()-t0:.0f}s]')

if __name__ == '__main__' and sys.argv[1] == 'verify':
    verify_saved(sys.argv[2], sys.argv[3]); sys.exit(0)
if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    d = json.loads(gzip.decompress(open(src, 'rb').read()))
    facets = [tuple(sorted(f)) for f in d['facets']]
    assert len(set(facets)) == len(facets)
    links = {}
    for f in facets:
        for v in f: links.setdefault(v, []).append(tuple(x for x in f if x != v))
    t0 = time.time(); res = {}; flips = 0; worst = 0
    for n, (v, L) in enumerate(sorted(links.items())):
        rng = random.Random(v); log = None
        for attempt in range(20):
            log = reduce3(L, rng)
            if log is not None: break
        assert log is not None, f'vertex {v}: link not reduced'
        assert verify3(L, log), f'vertex {v}: replay failed'
        res[v] = log; flips += len(log); worst = max(worst, len(L))
        if n % 500 == 0: print(f'  {n}/{len(links)} links certified [{time.time()-t0:.0f}s]', flush=True)
    json.dump({'vertices': len(links), 'flips_total': flips, 'largest_link_tets': worst, 'links': res}, open(out, 'w'))
    print(f'PASS: all {len(links)} vertex links are PL 3-spheres (bistellar to the boundary of the 4-simplex);'
          f' {flips} flips, largest link {worst} tetrahedra [{time.time()-t0:.0f}s]', flush=True)

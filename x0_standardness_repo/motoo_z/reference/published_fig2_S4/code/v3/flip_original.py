"""Crossing changes of a PD link (singles and pairs) whose exterior volume matches a target volume."""
import sys, json, time, itertools, numpy as np, snappy
fn, target = sys.argv[1], float(sys.argv[2]); maxk = int(sys.argv[3]) if len(sys.argv) > 3 else 2
d = json.load(open(fn))
pd = [list(x) for x in d['pd_code']]
blocks = d['component_labels'] if 'component_labels' in d else [c['arc_labels'] for c in d['components']]
comp = {l: i for i, b in enumerate(blocks) for l in b}
nxt = {}
for b in blocks:
    for i, l in enumerate(b):
        nxt[l] = b[(i+1) % len(b)]
def flip(x):
    a, b, c, dd = x
    if nxt[dd] == b:
        return [dd, a, b, c]
    if nxt[b] == dd:
        return [b, c, dd, a]
    raise ValueError(x)
def vol(p):
    try:
        M = snappy.Link(p).exterior()
        return M.volume()
    except Exception:
        return float('nan')
v0 = vol(pd); print('base volume', v0, 'target', target, 'crossings', len(pd), flush=True)
kinds = [(comp[x[0]], comp[x[1]]) for x in pd]
hits = []; t0 = time.time(); n = 0
for k in range(1, maxk+1):
    for S in itertools.combinations(range(len(pd)), k):
        q = [flip(x) if i in S else x for i, x in enumerate(pd)]
        v = vol(q); n += 1
        if abs(v-target) < 1e-5:
            hits.append((S, [kinds[i] for i in S], v)); print('  HIT', hits[-1], flush=True)
        if n % 300 == 0:
            print(f'  {n} tested, {time.time()-t0:.0f}s', flush=True)
print('tested', n, 'hits', hits)

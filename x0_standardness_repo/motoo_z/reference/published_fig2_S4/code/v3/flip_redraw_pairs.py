"""In a redraw PD, find pairs of crossings (framed comp F with dotted D1, D2) whose change gives a target volume."""
import sys, json, itertools, numpy as np, snappy
fn, target = sys.argv[1], float(sys.argv[2]); F = sys.argv[3]; Ds = sys.argv[4].split(',')
d = json.load(open(fn)); order = d['order']
pd = [list(x) for x in d['pd_code']]; blocks = d['component_labels']
comp = {l: order[i] for i, b in enumerate(blocks) for l in b}
nxt = {}
for b in blocks:
    for i, l in enumerate(b):
        nxt[l] = b[(i+1) % len(b)]
def flip(x):
    a, b, c, dd = x
    return [dd, a, b, c] if nxt[dd] == b else [b, c, dd, a]
def vol(p):
    try:
        return snappy.Link(p).exterior().volume()
    except Exception:
        return float('nan')
cands = {D: [k for k, x in enumerate(pd) if {comp[x[0]], comp[x[1]]} == {F, D}] for D in Ds}
print({D: len(v) for D, v in cands.items()})
hits = []
for i in cands[Ds[0]]:
    for j in cands[Ds[1]]:
        q = [flip(x) if k in (i, j) else x for k, x in enumerate(pd)]
        v = vol(q)
        if abs(v-target) < 1e-5:
            hits.append((i, j)); print('HIT', i, j, d['crossing_points'][i], d['crossing_points'][j], 'over before:', comp[pd[i][1]], comp[pd[j][1]], v, flush=True)
            d2 = dict(d); d2['pd_code'] = q; json.dump(d2, open(sys.argv[5].replace('.json', f'_{i}_{j}.json'), 'w'))
print('hits', hits)

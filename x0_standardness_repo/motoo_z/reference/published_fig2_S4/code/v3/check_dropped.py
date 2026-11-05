"""For each piece the linker dropped: does any curve cross it, and how far is it from its own colour's curve?"""
import sys, json, pickle, numpy as np
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from link import link, seg_hits
from scipy.spatial import cKDTree
trace_pkl, build_pkl, cfg_fn = sys.argv[1:4]
T = pickle.load(open(trace_pkl, 'rb')); B = pickle.load(open(build_pkl, 'rb')); cfg = json.load(open(cfg_fn))
res = T['res']; curves = B['curves']
A = []; Bs = []
for n, r in res.items():
    for p in r['pieces']:
        A.append(p['P'][:-1]); Bs.append(p['P'][1:])
segs = (np.vstack(A), np.vstack(Bs))
for nm, r in res.items():
    kw = dict(cfg.get('link_kw', {}).get(nm, {}))
    for key in ('forced', 'forbidden'):
        if key in kw:
            kw[key] = set(tuple(x) for x in kw[key])
    ends, mate, comps, links, dropped = link(r['pieces'], segs, **kw)
    own = cKDTree(curves[nm]['P'])
    for i in sorted(dropped):
        P = r['pieces'][i]['P']
        hits = {}
        for on, cv in curves.items():
            Q = cv['P']; Qa = Q; Qb = np.roll(Q, -1, axis=0)
            if on == nm:
                continue
            h = sum(seg_hits(P[k], P[k+1], Qa, Qb) for k in range(len(P)-1))
            if h:
                hits[on] = h
        d, _ = own.query(P)
        print(f'{nm} dropped L={r["pieces"][i]["L"]:.1f} at {P[len(P)//2].round(1).tolist()}: crossed by {hits or "nothing"}; max dist to own curve {d.max():.1f}')

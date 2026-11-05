"""Median rendered stroke width (pt) along each traced piece, to flag faint pencil touches."""
import sys, json, pickle, numpy as np
from scipy.ndimage import distance_transform_edt
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from masks import render_masks
fn, trace_pkl, cfg_fn, out = sys.argv[1:5]
T = pickle.load(open(trace_pkl, 'rb')); cfg = json.load(open(cfg_fn)); res = T['res']
m, H, s = render_masks(fn, ppt=6.0, exclude=set(cfg['labels']), recolour={int(k): v for k, v in cfg.get('recolour', {}).items()})
W = {}
for nm, r in res.items():
    dt = distance_transform_edt(m[nm])
    ws = []
    for p in r['pieces']:
        P = p['P']; c = np.clip((P[:, 0]*s).round().astype(int), 0, dt.shape[1]-1); rr = np.clip(((H-P[:, 1])*s).round().astype(int), 0, dt.shape[0]-1)
        w = 2*dt[rr, c]/s
        ws.append(float(np.median(w)))
    W[nm] = ws
    a = np.array(ws)
    print(nm, 'pieces', len(ws), 'median width quantiles', np.percentile(a, [0, 5, 25, 50, 75, 100]).round(2))
    for i, (p, w) in enumerate(zip(r['pieces'], ws)):
        if w < 0.45*np.median(a):
            print(f'   thin piece {i}: width {w:.2f} L={p["L"]:.1f} from {p["P"][0].round(1).tolist()} to {p["P"][-1].round(1).tolist()}')
json.dump(W, open(out, 'w'))

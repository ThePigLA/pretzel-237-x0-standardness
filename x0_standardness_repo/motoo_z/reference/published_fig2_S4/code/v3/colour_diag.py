"""Diagnose the linking of one colour: pieces, matching links (per component), dropped pieces."""
import sys, json, pickle, numpy as np, pymupdf, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from link import link
from masks import COLS
fn, trace_pkl, cfg_fn, nm, out = sys.argv[1:6]
win = json.loads(sys.argv[6]) if len(sys.argv) > 6 else None
T = pickle.load(open(trace_pkl, 'rb')); cfg = json.load(open(cfg_fn)); res = T['res']
A = []; B = []
for n, r in res.items():
    for p in r['pieces']:
        A.append(p['P'][:-1]); B.append(p['P'][1:])
segs = (np.vstack(A), np.vstack(B))
pcs = res[nm]['pieces']
kw = dict(cfg.get('link_kw', {}).get(nm, {}))
for key in ('forced', 'forbidden'):
    if key in kw:
        kw[key] = set(tuple(x) for x in kw[key])
ends, mate, comps, links, dropped = link(pcs, segs, **kw)
comp_of = {}
for ci, c in enumerate(comps):
    for n_ in c['nodes']:
        comp_of[n_//2] = ci
print(nm, 'pieces', len(pcs), 'components', [(c['npieces'], c['closed']) for c in comps])
for i in sorted(dropped):
    print('  dropped piece', i, 'L=%.1f' % pcs[i]['L'], pcs[i]['P'][0].round(1).tolist(), pcs[i]['P'][-1].round(1).tolist())
for ci, c in enumerate(comps):
    if c['npieces'] < 6:
        for n_ in c['nodes'][::2]:
            i = n_//2
            print(f'  comp {ci} piece {i} L={pcs[i]["L"]:.1f}', pcs[i]['P'][0].round(1).tolist(), pcs[i]['P'][-1].round(1).tolist())
doc = pymupdf.open(fn); page = doc[0]; H = page.rect.height
if win:
    x0, y0, x1, y1 = win
else:
    allP = np.vstack([p['P'] for p in pcs]); x0, y0 = allP.min(0)-10; x1, y1 = allP.max(0)+10
clip = pymupdf.Rect(x0, H-y1, x1, H-y0)
pix = page.get_pixmap(matrix=pymupdf.Matrix(3, 3), clip=clip, annots=True)
img = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3]
fig, ax = plt.subplots(figsize=(22, 22*(y1-y0)/(x1-x0)))
ax.imshow(img, extent=[x0, x1, y0, y1], alpha=0.35)
cmap = plt.get_cmap('tab10')
for i, p in enumerate(pcs):
    col = 'k' if i in dropped else cmap(comp_of.get(i, 9) % 10)
    ax.plot(p['P'][:, 0], p['P'][:, 1], color=col, lw=1.4)
    m = p['P'][len(p['P'])//2]
    if x0 <= m[0] <= x1 and y0 <= m[1] <= y1:
        ax.text(m[0], m[1], str(i), fontsize=7, color=col)
for a, b in mate.items():
    if a < b and ends[a]['i'] != ends[b]['i']:
        pa, pb = ends[a]['p'], ends[b]['p']
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color='m', lw=0.9, ls='--')
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect('equal')
plt.tight_layout(); plt.savefig(out, dpi=70)

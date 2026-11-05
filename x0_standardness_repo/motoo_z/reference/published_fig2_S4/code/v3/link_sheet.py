"""Sheet of windows: page crop with one colour's pieces (thin black, numbered) and matching links (magenta)."""
import sys, json, pickle, numpy as np, pymupdf, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from link import link
fn, trace_pkl, cfg_fn, out = sys.argv[1:5]; spots = json.loads(sys.argv[5])
T = pickle.load(open(trace_pkl, 'rb')); cfg = json.load(open(cfg_fn)); res = T['res']
A = []; B = []
for n, r in res.items():
    for p in r['pieces']:
        A.append(p['P'][:-1]); B.append(p['P'][1:])
segs = (np.vstack(A), np.vstack(B))
cache = {}
def L(nm):
    if nm not in cache:
        kw = dict(cfg.get('link_kw', {}).get(nm, {}))
        for key in ('forced', 'forbidden'):
            if key in kw:
                kw[key] = set(tuple(x) for x in kw[key])
        cache[nm] = link(res[nm]['pieces'], segs, **kw)
    return cache[nm]
doc = pymupdf.open(fn); page = doc[0]; H = page.rect.height
cols = min(3, len(spots)); rows = int(np.ceil(len(spots)/cols))
fig, axes = plt.subplots(rows, cols, figsize=(7*cols, 7*rows)); axes = np.atleast_1d(axes).ravel()
for ax, (nm, x, y, h) in zip(axes, spots):
    clip = pymupdf.Rect(x-h, H-y-h, x+h, H-y+h)
    pix = page.get_pixmap(matrix=pymupdf.Matrix(8, 8), clip=clip, annots=True)
    img = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3]
    ax.imshow(img, extent=[x-h, x+h, y-h, y+h], alpha=0.6)
    ends, mate, comps, links, dropped = L(nm)
    for i, p in enumerate(res[nm]['pieces']):
        P = p['P']; m = (np.abs(P[:, 0]-x) < h*1.2) & (np.abs(P[:, 1]-y) < h*1.2)
        if not m.any():
            continue
        ax.plot(P[:, 0], P[:, 1], color='r' if i in dropped else 'k', lw=0.8)
        ax.plot(*P[0], 'k.', ms=4); ax.plot(*P[-1], 'kx', ms=5)
        k = np.nonzero(m)[0][len(np.nonzero(m)[0])//2]; ax.text(P[k, 0], P[k, 1], str(i), fontsize=8)
    for a, b in mate.items():
        if a < b and ends[a]['i'] != ends[b]['i']:
            pa, pb = ends[a]['p'], ends[b]['p']
            ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color='m', lw=1.2, ls='--')
    ax.set_xlim(x-h, x+h); ax.set_ylim(y-h, y+h); ax.set_title(f'{nm} ({x},{y})', fontsize=10)
for ax in axes[len(spots):]:
    ax.axis('off')
plt.tight_layout(); plt.savefig(out, dpi=55)

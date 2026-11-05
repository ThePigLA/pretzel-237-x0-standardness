"""Contact sheet of page crops (full render) at given y-up points: crops.py pdf out.png '[[x,y,half],...]'"""
import sys, json, numpy as np, pymupdf, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
fn, out = sys.argv[1], sys.argv[2]; spots = json.loads(sys.argv[3])
doc = pymupdf.open(fn); page = doc[0]; H = page.rect.height
cols = min(4, len(spots)); rows = int(np.ceil(len(spots)/cols))
fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 5*rows)); axes = np.atleast_1d(axes).ravel()
for ax, (x, y, h) in zip(axes, spots):
    clip = pymupdf.Rect(x-h, H-y-h, x+h, H-y+h)
    pix = page.get_pixmap(matrix=pymupdf.Matrix(8, 8), clip=clip, annots=True)
    img = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3]
    ax.imshow(img, extent=[x-h, x+h, y-h, y+h]); ax.set_title(f'({x},{y})', fontsize=9)
    ax.tick_params(labelsize=6)
for ax in axes[len(spots):]:
    ax.axis('off')
plt.tight_layout(); plt.savefig(out, dpi=60)

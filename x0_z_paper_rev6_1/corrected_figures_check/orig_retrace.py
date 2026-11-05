"""Independent re-transcription of Motoo's original vector Figure (strokes only) with the raster pipeline."""
import sys, json, pickle, numpy as np, pymupdf, networkx as nx
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from link import link
from build_pd import find_crossings, pd_code, linking_matrix, dotted_report, trim_ends
from pd_genus import genus

pdf, arch_json, out = sys.argv[1], sys.argv[2], sys.argv[3]
y0, y1 = float(sys.argv[4]), float(sys.argv[5])          # vertical band of the figure (PDF y-down)
src = pymupdf.open(pdf); page = src[0]
W, H = page.rect.width, page.rect.height
doc = pymupdf.open(); np_ = doc.new_page(width=W, height=H)
sh = np_.new_shape()
n = 0
for d in page.get_drawings():
    if d['type'] != 's' or d['rect'].y0 < y0 or d['rect'].y1 > y1:
        continue
    for it in d['items']:
        if it[0] == 'l':
            sh.draw_line(it[1], it[2])
        elif it[0] == 'c':
            sh.draw_bezier(it[1], it[2], it[3], it[4])
        elif it[0] == 're':
            sh.draw_rect(it[1])
    sh.finish(color=(0, 0, 0), width=d['width'], closePath=False, lineCap=1)
    n += 1
sh.commit()
def bez(p, m=24):
    p = np.asarray(p, float); t = np.linspace(0, 1, m)[:, None]
    return (1-t)**3*p[0]+3*(1-t)**2*t*p[1]+3*(1-t)*t**2*p[2]+t**3*p[3]
pieces = []
for d in page.get_drawings():
    if d['type'] != 's' or d['rect'].y0 < y0 or d['rect'].y1 > y1:
        continue
    Q = []
    for it in d['items']:
        if it[0] == 'l':
            seg = np.array([[it[1].x, it[1].y], [it[2].x, it[2].y]])
        elif it[0] == 'c':
            seg = bez([[q.x, q.y] for q in it[1:5]])
        else:
            continue
        if Q and np.linalg.norm(Q[-1][-1]-seg[0]) < 1e-6:
            seg = seg[1:]
        Q.append(seg)
    P = np.vstack(Q)
    P[:, 1] = H-P[:, 1]   # to y-up like the rest of the pipeline
    L = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
    # a closed stroke (e.g. an unbroken circle) is kept as is
    pieces.append(dict(P=P, L=L, closed=np.linalg.norm(P[0]-P[-1]) < 1e-3))
print('vector strokes', len(pieces), 'closed', sum(p['closed'] for p in pieces))
A = np.vstack([p['P'][:-1] for p in pieces]); B = np.vstack([p['P'][1:] for p in pieces])
forced_xy = json.loads(sys.argv[6]) if len(sys.argv) > 6 else []
from link import tangent
endpts = []
for i, p in enumerate(pieces):
    endpts += [p['P'][0], p['P'][-1]]
endpts = np.array(endpts)
forced = set()
for (xa, ya), (xb, yb) in forced_xy:
    a = int(np.argmin(np.linalg.norm(endpts-np.array([xa, H-ya]), axis=1)))
    b = int(np.argmin(np.linalg.norm(endpts-np.array([xb, H-yb]), axis=1)))
    forced.add((min(a, b), max(a, b)))
print('forced pairs', forced)
ends, mate, comps, links, dropped = link(pieces, (A, B), forced=forced)
print('components', [(c['npieces'], c['closed']) for c in comps], 'dropped', len(dropped))
curves = {}
for ci, c in enumerate(comps):
    if not c['closed']:
        continue
    nodes = set(c['nodes']); start = min(nodes); cur = start; pts = []; ink = []; seen = 0
    while True:
        i, e = ends[cur]['i'], ends[cur]['e']
        P = trim_ends(pieces[i]['P'], 0.3); P = P[::-1] if e == 1 else P
        for k in range(len(P)):
            pts.append(P[k]); ink.append(k < len(P)-1)
        nxt = mate[2*i+(1-e)]
        seen += 1
        if nxt == start or seen > 10000:
            break
        cur = nxt
    Pa = np.array(pts) + np.random.default_rng(ci).uniform(-0.02, 0.02, (len(pts), 2))
    curves[f'k{ci}'] = dict(P=Pa, ink=np.array(ink, bool))
print('closed curves', len(curves))
# identify with archive components (archive geometry is PDF y-down)
arch = json.load(open(arch_json))
strokes = {st['id']: st for st in arch['geometry']['strokes']}
def bez(p, m=10):
    p = np.asarray(p, float); t = np.linspace(0, 1, m)[:, None]
    return (1-t)**3*p[0]+3*(1-t)**2*t*p[1]+3*(1-t)*t**2*p[2]+t**3*p[3]
apts = {}
for comp in arch['components']:
    Q = []
    for sc in comp['oriented_stroke_cycle']:
        for itm in strokes[sc['stroke']]['items']:
            Q.append(bez(itm['points']) if itm['type'] == 'c' else np.asarray(itm['points'], float))
    apts[comp['index']] = np.vstack(Q)
from scipy.spatial import cKDTree
trees = {k: cKDTree(v) for k, v in apts.items()}
ident = {}
for nm, cv in curves.items():
    Pd = np.column_stack([cv['P'][:, 0], H-cv['P'][:, 1]])
    scores = {k: float(np.mean(t.query(Pd)[0])) for k, t in trees.items()}
    best = min(scores, key=scores.get)
    ident[nm] = best
    print(nm, '-> archive C%d' % best, 'mean dist %.2f' % scores[best])
# per-piece archive component, and links that join different archive components
piece_comp = []
for p in pieces:
    Pd = np.column_stack([p['P'][:, 0], H-p['P'][:, 1]])
    sc = {k: float(np.mean(t.query(Pd)[0])) for k, t in trees.items()}
    piece_comp.append(min(sc, key=sc.get))
for l, a, b in links:
    if l.get('drop'): continue
    ca, cb = piece_comp[ends[a]['i']], piece_comp[ends[b]['i']]
    if ca != cb:
        pa = ends[a]['p']; pb = ends[b]['p']
        print(f"CROSS-COMPONENT LINK: archive C{ca} end at ({pa[0]:.1f},{H-pa[1]:.1f}) <-> C{cb} end at ({pb[0]:.1f},{H-pb[1]:.1f}) cost {l['cost']:.1f} d {l['d']:.1f}")
        alt = []
        for e2 in range(len(ends)):
            if e2 in (a, b): continue
            if piece_comp[ends[e2]['i']] in (ca, cb):
                d2 = np.linalg.norm(ends[e2]['p']-pa)
                if d2 < 15: alt.append((round(float(d2), 1), f"C{piece_comp[ends[e2]['i']]}", (round(ends[e2]['p'][0], 1), round(H-ends[e2]['p'][1], 1))))
        print('    other nearby ends:', sorted(alt)[:6])
assert len(set(ident.values())) == len(ident) == len(arch['components']), 'component identification failed'
X = find_crossings(curves)
for x in X:
    if x['status'] == 'ok':
        continue
    if x['a'] == x['b'] and x['status'] == 'ink-ink':
        x['over'] = 'a' if x['ia']+x['ta'] < x['ib']+x['tb'] else 'b'
    elif x['status'] == 'gap-gap':
        Pa = curves[x['a']]['P']; Pb = curves[x['b']]['P']
        la = np.linalg.norm(Pa[(x['ia']+1) % len(Pa)]-Pa[x['ia']]); lb = np.linalg.norm(Pb[(x['ib']+1) % len(Pb)]-Pb[x['ib']])
        x['over'] = 'a' if la < lb else 'b'
    else:
        raise SystemExit(f"unresolved {x['status']} {x['a']} {x['b']} {x['pt']}")
print('crossings', len(X), 'ambiguous', sum(x['status'] != 'ok' for x in X))
order = sorted(curves, key=lambda nm: ident[nm])
pd, passes = pd_code(curves, X, order)
print('genus', genus(pd))
dec = [arch['decorations'][ident[nm]] for nm in order]
comp_labels = []; base = 0
for nm in order:
    m = len(passes[nm]); comp_labels.append(list(range(base, base+m))); base += m
json.dump(dict(pd_code=pd, decorations=dec, order=[f'C{ident[nm]}' for nm in order], component_labels=comp_labels),
          open(out, 'w'))
print('linking matrix (archive order):'); print(linking_matrix(curves, X, order))
print('archive linking matrix:'); print(np.array(arch['linking_matrix']))

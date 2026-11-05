"""Resolve ambiguous crossings by fixed rules, build PD + decorations, report Katie position and linking."""
import sys, json, pickle, numpy as np
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from build_pd import pd_code, linking_matrix, dotted_report

build_pkl, cfg_fn, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
B = pickle.load(open(build_pkl, 'rb')); cfg = json.load(open(cfg_fn))
curves, X = B['curves'], B['crossings']
overrides = {tuple(map(float, k.split(','))): v for k, v in cfg.get('crossing_overrides', {}).items()}


def seglen(name, i):
    P = curves[name]['P']; return float(np.linalg.norm(P[(i+1) % len(P)]-P[i]))


_comp = None
def visible_colour(pt, a, b):
    global _comp
    import pymupdf
    from masks import COLS
    if _comp is None:
        doc = pymupdf.open(cfg['pdf']); page = doc[0]; Hh = page.rect.height
        pix = page.get_pixmap(matrix=pymupdf.Matrix(8, 8), annots=True, alpha=False)
        _comp = (np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3].astype(float)/255, Hh)
    img, Hh = _comp
    r = int(round((Hh-pt[1])*8)); c = int(round(pt[0]*8))
    patch = img[max(r-2, 0):r+3, max(c-2, 0):c+3].reshape(-1, 3)
    score = {}
    for nm in (a, b):
        col = np.array(COLS[nm])
        score[nm] = float(np.min(np.linalg.norm(patch-col, axis=1)))
    return min(score, key=score.get), score

resolved = []
for x in X:
    key = None
    for k, v in overrides.items():
        if np.linalg.norm(np.array(k)-x['pt']) < 1.5:
            key = k
    if key is not None:
        v = overrides[key]
        x['over'] = v if v in ('a', 'b') else ('a' if x['a'] == v else 'b')
        assert v in ('a', 'b') or v in (x['a'], x['b']), (v, x['a'], x['b'])
        x['rule'] = 'override (%s over)' % v
    elif x['status'] == 'ink-ink' and x['a'] == x['b']:
        # descending rule: earlier pass along the curve goes over
        pa, pb = x['ia']+x['ta'], x['ib']+x['tb']
        x['over'] = 'a' if pa < pb else 'b'; x['rule'] = 'descending'
    elif x['status'] == 'gap-gap':
        la, lb = seglen(x['a'], x['ia']), seglen(x['b'], x['ib'])
        x['over'] = 'a' if la < lb else 'b'; x['rule'] = f'shorter-gap-over ({la:.1f} vs {lb:.1f})'
    elif x['status'] == 'ink-ink':
        top, sc = visible_colour(x['pt'], x['a'], x['b'])
        x['over'] = 'a' if top == x['a'] else 'b'
        x['rule'] = f"visible-on-top={top} (colour distances {', '.join(f'{k}:{v:.2f}' for k, v in sc.items())})"
    else:
        x['rule'] = 'ink-vs-gap'
    if x['status'] != 'ok':
        resolved.append((x['status'], x['a'], x['b'], np.round(x['pt'], 1).tolist(), x['rule'], x['over']))
unres = [r for r in resolved if r[5] is None]
print('resolved ambiguous crossings:')
for r in resolved:
    print('   ', r)
if unres:
    print('UNRESOLVED', unres); sys.exit(1)

order = cfg['order']
dec = ['x' if n in cfg['dotted'] else int(cfg['framings'][n]) for n in order]

# remove Reidemeister-I curls: self-crossings whose two passes are consecutive along the component
def r1_cleanup(X):
    removed = []
    while True:
        passes = {}
        for k, x in enumerate(X):
            passes.setdefault(x['a'], []).append((x['ia']+x['ta'], k))
            passes.setdefault(x['b'], []).append((x['ib']+x['tb'], k))
        hit = None
        for n, p in passes.items():
            p.sort(); ks = [k for _, k in p]; m = len(ks)
            for i in range(m):
                if ks[i] == ks[(i+1) % m] and X[ks[i]]['a'] == X[ks[i]]['b'] and m > 2:
                    hit = ks[i]; break
            if hit is not None:
                break
        if hit is None:
            return X, removed
        removed.append(np.round(X[hit]['pt'], 1).tolist())
        X = X[:hit]+X[hit+1:]

try:
    import regina
    pd0, _ = pd_code(curves, X, order)
    R0 = regina.Link.fromPD([[l+1 for l in x] for x in pd0])
    sys.path.insert(0, '/home/claude/motoo_followup/v2'); from pd_genus import genus
    print('BEFORE curl removal: crossings', len(pd0), 'planar (classical):', R0.isClassical(), genus(pd0))
except Exception as e:
    print('pre-check failed', e)
X, r1 = r1_cleanup(X)
print(f'removed {len(r1)} Reidemeister-I curls at', r1)
pd, passes = pd_code(curves, X, order)
pickle.dump(X, open(out_json.replace('.json', '_X.pkl'), 'wb'))
LM = linking_matrix(curves, X, order)
rep = dotted_report(passes, X, order, cfg['dotted'])
print('order', order); print('decorations', dec)
print('crossings', len(pd))
print('linking matrix (component order above):'); print(LM)
admissible = True
for n, r in rep.items():
    ok = r['self'] == 0 and r['mutual'] == 0 and r['transitions'] == 2
    admissible &= ok
    print(f"dotted {n}: {r['n']} crossings, self {r['self']}, with other dotted {r['mutual']}, over/under transitions {r['transitions']}  {r['seq']}  {'OK' if ok else 'NOT in Katie position'}")
# writhe per component (diagnostic)
comp_labels = []
base = 0
for n in order:
    m = len(passes[n]); comp_labels.append(list(range(base, base+m))); base += m
json.dump(dict(schema='motoo.redraw.v1', source=cfg.get('source'), order=order, pd_code=pd, decorations=dec,
               component_labels=comp_labels,
               pd_component_to_source=list(range(len(order))), linking_matrix=LM.tolist(),
               dotted_report=rep, katie_position=bool(admissible), ambiguous_resolved=resolved,
               crossing_points=[np.round(x['pt'], 2).tolist() for x in X]), open(out_json, 'w'), indent=1)
try:
    import regina
    RL = regina.Link.fromPD([[l+1 for l in x] for x in pd])
    print('planar (classical) diagram:', RL.isClassical(), '| regina crossings', RL.size(), 'components', RL.countComponents())
except Exception as e:
    print('regina check failed', e)
print('wrote', out_json, 'katie_position', admissible)

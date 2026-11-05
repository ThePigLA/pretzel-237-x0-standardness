"""Finger moves on a Katie-position diagram (geometric curves + resolved crossings).

A finger is a thin band from a framed component F (the 'base', default blue) along a centreline
    p (on F) --hop--> q (beside dotted circle D) --along D--> t (an over/under transition of D) --poke--> across D.
Its tip clasps D at t: the band crosses D once just on D's over-side of t (band under D) and once just on
D's under-side (band over D).  So D keeps one over-block and one under-block (Katie position is preserved),
and the move is a crossing change between F and D, slid along D to t.
Heights along D are forced (the band stays at D's height: over strands that pass under D, under strands that
pass over D).  Heights in the hop are free choices.  Other dotted circles are always crossed by the rule
'over an under-block, under an over-block'.
"""
import sys, itertools, numpy as np
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from build_pd import pd_code, linking_matrix, dotted_report

W = 0.2      # band half-width
D0 = 0.8     # centreline offset from D while sliding along it
D1 = 1.0     # how far the tip pokes through to the other side


def rot(v):
    return np.array([-v[..., 1], v[..., 0]]).T if v.ndim > 1 else np.array([-v[1], v[0]])


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v/np.where(n == 0, 1, n)


class Grid:
    def __init__(self, curves, cell=4.0):
        self.cell = cell; self.cells = {}; self.curves = curves
        for name, c in curves.items():
            P = c['P']; A = P; B = np.roll(P, -1, axis=0)
            lo = np.floor(np.minimum(A, B)/cell).astype(int); hi = np.floor(np.maximum(A, B)/cell).astype(int)
            for j in range(len(P)):
                for gx in range(lo[j, 0], hi[j, 0]+1):
                    for gy in range(lo[j, 1], hi[j, 1]+1):
                        self.cells.setdefault((gx, gy), []).append((name, j))

    def candidates(self, a0, a1):
        lo = np.floor(np.minimum(a0, a1)/self.cell).astype(int); hi = np.floor(np.maximum(a0, a1)/self.cell).astype(int)
        out = set()
        for gx in range(lo[0], hi[0]+1):
            for gy in range(lo[1], hi[1]+1):
                out.update(self.cells.get((gx, gy), ()))
        return out


def seg_x(a0, a1, b0, b1):
    r = a1-a0; s = b1-b0; den = r[0]*s[1]-r[1]*s[0]
    if abs(den) < 1e-12:
        return None
    w = b0-a0
    t = (w[0]*s[1]-w[1]*s[0])/den; u = (w[0]*r[1]-w[1]*r[0])/den
    if 0 <= t < 1 and 0 <= u < 1:
        return t, u
    return None


def polyline_hits(Q, grid, skip=None):
    """Crossings of open polyline Q with the indexed curves: list of (k, t, name, j, u, pt)."""
    hits = []
    for k in range(len(Q)-1):
        a0, a1 = Q[k], Q[k+1]
        for name, j in grid.candidates(a0, a1):
            if skip and skip(name, j):
                continue
            P = grid.curves[name]['P']; b0 = P[j]; b1 = P[(j+1) % len(P)]
            r = seg_x(a0, a1, b0, b1)
            if r:
                hits.append((k, r[0], name, j, r[1], a0+r[0]*(a1-a0)))
    hits.sort(key=lambda h: h[0]+h[1])
    return hits


class Base:
    def __init__(self, curves, X, order, dotted, base='blue'):
        self.curves, self.X, self.order, self.dotted, self.base = curves, X, order, dotted, base
        self.grid = Grid(curves)
        pd, passes = pd_code(curves, X, order)
        self.passes = passes
        self.lk0 = linking_matrix(curves, X, order)
        # blue segments that carry a crossing
        busy = set()
        for x in X:
            if x['a'] == base:
                busy.add(x['ia'])
            if x['b'] == base:
                busy.add(x['ib'])
        self.busy = busy
        # pass intervals ('edges') of each curve, for de-duplication
        self.edge_pos = {n: np.array([p for p, _, _ in passes[n]]) for n in order}

    def N(self, n):
        return len(self.curves[n]['P'])

    def point(self, n, pos):
        P = self.curves[n]['P']; N = len(P); i = int(np.floor(pos)) % N; t = pos-np.floor(pos)
        return P[i]+t*(P[(i+1) % N]-P[i])

    def tangent(self, n, pos, h=1.0):
        # forward tangent by arc-length-ish finite difference in index space
        P = self.curves[n]['P']; N = len(P)
        i = int(np.floor(pos)) % N
        a = P[(i-3) % N]; b = P[(i+4) % N]
        return unit(b-a)

    def edge_id(self, n, pos):
        e = self.edge_pos[n]
        return int(np.searchsorted(e, pos % self.N(n)))

    def block(self, n, pos):
        """role of curve n at pos ('over'/'under') if both neighbouring passes agree, else None."""
        p = self.passes[n]; e = self.edge_pos[n]; m = len(p)
        k = int(np.searchsorted(e, pos % self.N(n)))
        r1 = p[(k-1) % m][2]; r2 = p[k % m][2]
        return r1 if r1 == r2 else None

    def transitions(self, n):
        p = self.passes[n]; m = len(p); N = self.N(n); out = []
        for k in range(m):
            a, b = p[k], p[(k+1) % m]
            if a[2] != b[2]:
                pa, pb = a[0], b[0]
                if pb < pa:
                    pb += N
                out.append(dict(pos=((pa+pb)/2) % N, before=a[2], after=b[2], lo=pa, hi=pb))
        return out


def along_path(B, D, pos_q, pos_t, direction, side, step=5):
    """Centreline beside D from pos_q to pos_t (moving in index direction 'direction'), offset to 'side'."""
    N = B.N(D)
    if direction > 0:
        span = (pos_t-pos_q) % N
        poss = [pos_q] + [pos_q+k for k in np.arange(np.ceil(0) + 1, span, step)] + [pos_q+span]
    else:
        span = (pos_q-pos_t) % N
        poss = [pos_q] + [pos_q-k for k in np.arange(1, span, step)] + [pos_q-span]
    pts = []
    for s in poss:
        pts.append(B.point(D, s) + side*D0*rot(B.tangent(D, s)))
    return np.array(pts), poss


def offset_open(Q, d):
    T = np.zeros_like(Q)
    T[1:-1] = unit(Q[2:]-Q[:-2]); T[0] = unit(Q[1]-Q[0]); T[-1] = unit(Q[-1]-Q[-2])
    return Q + d*rot(T)


def build_band(G, s_out):
    out = offset_open(G, s_out*W)
    back = offset_open(G, -s_out*W)[::-1]
    T = G[-1]; g = unit(G[-1]-G[-2]); n = rot(g)
    a0 = np.arctan2(*(s_out*n)[::-1]); a1 = np.arctan2(*(-s_out*n)[::-1]); ag = np.arctan2(g[1], g[0])
    # semicircle from a0 to a1 passing through ag
    def ang_between(a, b, c):
        d1 = (c-a) % (2*np.pi); d2 = (b-a) % (2*np.pi)
        return d1 < d2
    if ang_between(a0, a1, ag):
        angs = a0 + np.linspace(0, (a1-a0) % (2*np.pi), 7)
    else:
        angs = a0 - np.linspace(0, (a0-a1) % (2*np.pi), 7)
    cap = np.array([T + W*np.array([np.cos(a), np.sin(a)]) for a in angs[1:-1]])
    return np.vstack([out, cap, back])


def apply_finger(B, spec, heights=None, want_events=False):
    """spec: dict(D, t (transition dict), side, direction, pos_q, base_vertex i).
    Returns (curves', X', info) or (None, None, reason)."""
    D = spec['D']; base = B.base; P = B.curves[base]['P']; N = len(P); i = spec['i']
    if i in B.busy or (i-1) % N in B.busy:
        return None, None, 'busy base vertex'
    p = P[i]
    along, poss = along_path(B, D, spec['pos_q'], spec['t']['pos'], spec['direction'], spec['side'])
    tpt = B.point(D, spec['t']['pos']); tn = rot(B.tangent(D, spec['t']['pos']))
    tip = tpt - spec['side']*D1*tn
    G = np.vstack([p[None], along, tip[None]])
    # section labels by centreline vertex: 0 = hop (segment 0), 1..len(along)-1 = along, last = poke
    n_along = len(along)
    g0 = unit(G[1]-G[0]); tb = unit(P[(i+1) % N]-P[(i-1) % N])
    if abs(np.dot(g0, tb)) > 0.85:
        return None, None, 'hop tangent to base'
    s_out = 1 if np.dot(rot(g0), -tb) > 0 else -1
    band = build_band(G, s_out)
    # centreline crossings (for heights and event grouping)
    chits = polyline_hits(G, B.grid)
    events = []
    for (k, t, name, j, u, pt) in chits:
        sec = 'hop' if k == 0 else ('poke' if k == len(G)-2 else 'along')
        events.append(dict(k=k, t=t, name=name, j=j, u=u, pt=pt, sec=sec))
    # the poke must cross only D, exactly once
    pk = [e for e in events if e['sec'] == 'poke']
    if len(pk) != 1 or pk[0]['name'] != D:
        return None, None, 'poke crosses %s' % [e['name'] for e in pk]
    free = [e for e in events if e['sec'] == 'hop' and e['name'] not in B.dotted]
    if want_events:
        return free
    if heights is None:
        heights = [1]*len(free)
    if len(heights) != len(free):
        return None, None, 'heights length'
    for e, h in zip(free, heights):
        e['band_over'] = bool(h)
    # forced heights
    for e in events:
        if e['sec'] == 'hop' and e['name'] not in B.dotted:
            continue
        if e['name'] in B.dotted and not (e['sec'] == 'poke'):
            r = B.block(e['name'], e['j']+e['u'])
            if r is None:
                return None, None, 'crosses dotted at a transition'
            e['band_over'] = (r == 'under')
        elif e['sec'] == 'along':
            # nearest existing crossing between D and this strand
            best = None
            for x in B.X:
                if {x['a'], x['b']} == {D, e['name']} and x['a'] != x['b']:
                    dd = np.linalg.norm(x['pt']-e['pt'])
                    if best is None or dd < best[0]:
                        Dover = (x['over'] == 'a') == (x['a'] == D)
                        best = (dd, Dover)
            if best is None or best[0] > 3.0:
                return None, None, 'along-D crossing with %s not near a D crossing' % e['name']
            e['band_over'] = not best[1]
    # new base polyline
    ins = band
    m = len(ins)
    newP = np.vstack([P[:i], ins, P[i+1:]])
    sec_lo, sec_hi = i-1, i+m-1          # new segments forming the inserted section (inclusive)
    def remap(j):
        if j <= i-2:
            return j
        if j >= i+1:
            return j+m-1
        raise ValueError('crossing on replaced segment')
    sectionQ = np.vstack([P[(i-1) % N][None], ins, P[(i+1) % N][None]])
    def skip(name, j):
        return name == base and j in ((i-1) % N, i)
    bh = polyline_hits(sectionQ, B.grid, skip)
    # self-crossings of the section
    for a in range(len(sectionQ)-1):
        for b in range(a+2, len(sectionQ)-1):
            if seg_x(sectionQ[a], sectionQ[a+1], sectionQ[b], sectionQ[b+1]):
                return None, None, 'band self-crossing'
    newX = []
    for x in B.X:
        y = dict(x)
        if y['a'] == base:
            y['ia'] = remap(y['ia'])
        if y['b'] == base:
            y['ib'] = remap(y['ib'])
        newX.append(y)
    for (k, t, name, j, u, pt) in bh:
        # match to a centreline event of the same curve
        cands = [e for e in events if e['name'] == name and np.linalg.norm(e['pt']-pt) < 1.6]
        if not cands:
            return None, None, 'band crossing with %s at %s has no centreline event' % (name, np.round(pt, 1).tolist())
        e = min(cands, key=lambda e: np.linalg.norm(e['pt']-pt))
        if e['sec'] == 'poke':
            # which side of the transition along D
            tt = spec['t']; pos = j+u; N_D = B.N(D)
            lo, hi = tt['lo'], tt['hi']; mid = (lo+hi)/2
            pp = pos if pos >= lo-1e-9 else pos+N_D
            side_role = tt['before'] if pp < mid else tt['after']
            band_over = (side_role == 'under')
        else:
            band_over = e['band_over']
        jj = j if name != base else remap(j)
        newX.append(dict(a=base, ia=sec_lo+k, ta=float(t), b=name, ib=jj, tb=float(u), pt=pt,
                         over='a' if band_over else 'b', rule='finger'))
    curves2 = {n: (dict(P=newP, ink=np.ones(len(newP), bool)) if n == base else c) for n, c in B.curves.items()}
    return curves2, newX, dict(n_new=len(bh), s_out=s_out, events=len(events))


def evaluate(B, curves2, X2):
    import snappy
    pd, passes = pd_code(curves2, X2, B.order)
    lk = linking_matrix(curves2, X2, B.order)
    rep = dotted_report(passes, X2, B.order, B.dotted)
    try:
        v = snappy.Link(pd).exterior().volume()
    except Exception:
        v = float('nan')
    return pd, lk, rep, v

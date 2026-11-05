"""Search finger moves (base component clasping dotted circle D) that turn a diagram into a target link.

usage: finger_search.py curves.pkl X.pkl diagram.json D target_volume target_lk out_prefix [Hmax] [max_free]
Writes out_prefix_hit.pkl (curves, X, spec, heights) for the first hit and prints progress.
"""
import sys, json, pickle, time, itertools, numpy as np
sys.path.insert(0, '/home/claude/motoo_followup/v3')
from finger import Base, apply_finger, evaluate, along_path

cpk, xpk, djs, D, tvol, tlk, outp = sys.argv[1:8]
tvol = float(tvol); tlk = int(tlk)
Hmax = float(sys.argv[8]) if len(sys.argv) > 8 else 60.0
max_free = int(sys.argv[9]) if len(sys.argv) > 9 else 5
base_name = 'blue'
Bk = pickle.load(open(cpk, 'rb')); curves = Bk['curves'] if 'curves' in Bk else Bk
X = pickle.load(open(xpk, 'rb')); dj = json.load(open(djs))
order = dj['order']; dotted = [n for n, d in zip(order, dj['decorations']) if d == 'x']
t0 = time.time()
B = Base(curves, X, order, dotted, base_name)
iD, ib = order.index(D), order.index(base_name)
print('base lk(blue,%s) =' % D, B.lk0[ib, iD], ' target', tlk, ' grid built %.1fs' % (time.time()-t0), flush=True)
P = curves[base_name]['P']; NB = len(P)
cand_base = np.array([i for i in range(0, NB, 6) if i not in B.busy and (i-1) % NB not in B.busy])
specs = []
for t in B.transitions(D):
    p = B.passes[D]; m = len(p); ND = B.N(D)
    for k in range(m):
        lo = p[k][0]; hi = p[(k+1) % m][0]
        if hi <= lo:
            hi += ND
        pos_q = ((lo+hi)/2) % ND
        for direction in (1, -1):
            for side in (1, -1):
                along, poss = along_path(B, D, pos_q, t['pos'], direction, side)
                alen = float(np.linalg.norm(np.diff(along, axis=0), axis=1).sum())
                qa = along[0]
                d = np.linalg.norm(P[cand_base]-qa, axis=1)
                for i, h in zip(cand_base[d < Hmax], d[d < Hmax]):
                    specs.append((alen+float(h), dict(D=D, t=t, side=side, direction=direction, pos_q=pos_q, i=int(i), qedge=k)))
specs.sort(key=lambda s: s[0])
print('specs', len(specs), '(%.1fs)' % (time.time()-t0), flush=True)
seen = set(); n_eval = 0; n_sig = 0; reasons = {}
for cost, spec in specs:
    free = apply_finger(B, spec, want_events=True)
    if isinstance(free, tuple):
        reasons[free[2].split(' ')[0]] = reasons.get(free[2].split(' ')[0], 0)+1
        continue
    if len(free) > max_free:
        continue
    sig = (spec['t']['pos'], spec['qedge'], spec['side'], spec['direction'], B.edge_id(base_name, spec['i']),
           tuple((e['name'], B.edge_id(e['name'], e['j']+e['u'])) for e in free))
    if sig in seen:
        continue
    seen.add(sig); n_sig += 1
    for hts in itertools.product([1, 0], repeat=len(free)):
        c2, X2, info = apply_finger(B, spec, list(hts))
        if c2 is None:
            reasons[info.split(' ')[0]] = reasons.get(info.split(' ')[0], 0)+1
            break
        pd, lk, rep, v = evaluate(B, c2, X2)
        n_eval += 1
        if lk[ib, iD] != tlk:
            continue
        katie = all(r['transitions'] <= 2 for r in rep.values()) if isinstance(rep, dict) else None
        if abs(v-tvol) < 1e-6:
            print(f'HIT cost {cost:.1f} spec {dict((k, spec[k]) for k in ("side", "direction", "pos_q", "i", "qedge"))} '
                  f't {spec["t"]["pos"]:.1f} heights {hts} volume {v} crossings {len(pd)} katie {katie}', flush=True)
            pickle.dump(dict(curves=c2, X=X2, spec=spec, heights=hts, pd=pd, order=order), open(outp+'_hit.pkl', 'wb'))
            print('evaluated', n_eval, 'signatures', n_sig, '%.0fs' % (time.time()-t0)); sys.exit(0)
    if n_sig % 200 == 0:
        print(f'  signatures {n_sig}, evaluated {n_eval}, cost {cost:.1f}, {time.time()-t0:.0f}s, rejects {reasons}', flush=True)
print('NO HIT; evaluated', n_eval, 'signatures', n_sig, 'rejects', reasons)

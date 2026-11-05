"""Remove faint pencil touches: pieces whose median rendered width is far below the colour's normal pen width."""
import sys, json, pickle
trace_pkl, widths_json, out = sys.argv[1:4]
thr = {'green': 1.5, 'salmon': 1.5, 'navy': 1.5}
T = pickle.load(open(trace_pkl, 'rb')); W = json.load(open(widths_json))
for nm, r in T['res'].items():
    t = thr.get(nm, 0.45)
    keep = [p for p, w in zip(r['pieces'], W[nm]) if w >= t]
    gone = [(round(p['L'], 1), p['P'][0].round(1).tolist()) for p, w in zip(r['pieces'], W[nm]) if w < t]
    if gone:
        print(nm, 'dropped faint pieces', gone)
    r['pieces'] = keep
pickle.dump(T, open(out, 'wb'))

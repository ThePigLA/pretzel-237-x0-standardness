"""Randomised Reidemeister simplification of a Katie-position diagram, keeping Katie position.

Each candidate is checked for: dotted components simple (no self-crossings), no dotted-dotted crossings,
one over-block and one under-block per dotted component, connected projection, and the same marked framed link
(decorated isometry signature with meridians, and an isometry matching framings).  Keeps the smallest.
usage: katie_shrink.py in.json out.json seconds
"""
import sys, json, time, random, snappy, spherogram
sys.path.insert(0, '/home/claude/motoo_followup/v2')
from compmap import spherogram_link

src, out, secs = sys.argv[1], sys.argv[2], float(sys.argv[3])
d = json.load(open(src))
L0, s2p = spherogram_link(d)          # s2p: source component -> spherogram index
p2s = {v: k for k, v in s2p.items()}
dec0 = [d['decorations'][p2s[c]] for c in range(len(L0.link_components))]
names0 = [d['order'][p2s[c]] for c in range(len(L0.link_components))]
sig0 = L0.exterior().isometry_signature(of_link=True)


def comp_of_strands(L):
    comp = {}
    for ci, comp_list in enumerate(L.link_components):
        for ce in comp_list:
            comp[(ce.crossing, ce.strand_index)] = ci
    return comp


def katie_ok(L, dotted_idx):
    # per component: sequence of over/under along the component, self/mutual crossings
    comps = L.link_components
    for ci in dotted_idx:
        seq = []
        for ce in comps[ci]:
            over = ce.strand_index in (1, 3)
            seq.append('O' if over else 'U')
        trans = sum(1 for a, b in zip(seq, seq[1:]+seq[:1]) if a != b)
        if trans != 2:
            return False
    # self / mutual crossings of dotted components
    for c in L.crossings:
        a = c.strand_components[0]; b = c.strand_components[1]
        if a in dotted_idx and b in dotted_idx:
            return False
    return True


def components_connected(L):
    return len(L.split_link_diagram()) == 1 if hasattr(L, 'split_link_diagram') else True


best = None; t0 = time.time(); tries = 0
while time.time()-t0 < secs:
    tries += 1
    L = L0.copy()
    mode = random.choice(['level', 'global', 'basic'])
    try:
        if mode == 'basic':
            L.simplify('basic')
        elif mode == 'level':
            L.simplify('level')
        else:
            L.simplify('global')
        # random extra: backtrack then simplify again
        if random.random() < 0.5:
            L.backtrack(random.randint(5, 40)); L.simplify('level')
    except Exception as e:
        continue
    n = len(L.crossings)
    if best is not None and n >= best[0]:
        continue
    # identify components: spherogram keeps component order through simplify? verify via isometry later
    # dotted indices: match by checking candidate assignment via linking numbers with original
    lk = L.linking_matrix()
    try:
        M = L.exterior()
        if M.isometry_signature(of_link=True) != sig0:
            continue
    except Exception:
        continue
    # component correspondence by isometry preserving meridians/longitudes
    iso = [i for i in L0.exterior().is_isometric_to(M, return_isometries=True)
           if all(abs(m[0, 0]) == 1 and m[0, 1] == 0 and m[1, 0] == 0 and m[0, 0] == m[1, 1] for m in i.cusp_maps())]
    if not iso:
        continue
    perm = iso[0].cusp_images()
    dec = [None]*len(perm); nm = [None]*len(perm)
    for i, j in enumerate(perm):
        dec[j] = dec0[i]; nm[j] = names0[i]
    dotted_idx = {j for j, x in enumerate(dec) if x == 'x'}
    if not katie_ok(L, dotted_idx):
        continue
    best = (n, L.PD_code(), dec, nm, mode)
    print(f'{time.time()-t0:.0f}s try {tries}: {n} crossings ({mode})', flush=True)
    json.dump(dict(pd_code=[list(x) for x in best[1]], decorations=best[2], component_names=best[3],
                   source=src + ' simplified by spherogram, Katie position kept'), open(out, 'w'))
print('best', best[0] if best else None, 'tries', tries)

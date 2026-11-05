"""Exact obstruction to a homeomorphism of #3(S^1 x S^2) carrying one family of 2-handle attaching circles
to the other, for the first drawing and the corrected drawing of a figure.

usage: python3 cover_obstruction.py OLD.json NEW.json --degrees 2,3 --output RESULT.json

Each record is a labelled PD code (0-based) with component_labels and decorations ('x' = dotted circle).
The link exterior is built with SnapPy; each dotted cusp is filled along its 0-slope (Seifert longitude),
which gives the exterior of the framed attaching circles in #3(S^1 x S^2).  The connected n-fold covers are
enumerated combinatorially (permutation representations of the filled fundamental group) and their first
homology is computed exactly (Smith normal form).  Different multisets show that the two exteriors are not
homeomorphic, so no homeomorphism of #3(S^1 x S^2) carries one family of attaching circles to the other.
Requires SnapPy (tested with 3.3.2).
"""
import argparse, hashlib, json
from collections import Counter
import snappy


def link_with_map(d):
    L = snappy.Link(d['pd_code'])
    comp_of = {l: i for i, b in enumerate(d['component_labels']) for l in b}
    pairs = {}
    for c in L.crossings:
        for pos in range(4):
            pairs.setdefault(c.strand_components[pos], set()).add(comp_of[d['pd_code'][c.label][pos]])
    if any(len(v) != 1 for v in pairs.values()):
        raise SystemExit('inconsistent component mapping')
    return L, {k: next(iter(v)) for k, v in pairs.items()}       # SnapPy cusp -> source component


def filled_exterior(path):
    d = json.load(open(path))
    L, cusp_to_src = link_with_map(d)
    M = L.exterior()
    dec = [d['decorations'][cusp_to_src[c]] for c in range(M.num_cusps())]
    M.dehn_fill([(0, 1) if x == 'x' else (0, 0) for x in dec])
    N = M.filled_triangulation()
    N.simplify()
    return d, N, sum(1 for x in dec if x == 'x')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('old'); ap.add_argument('new')
    ap.add_argument('--degrees', default='2,3'); ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = dict(schema='x0_z.cover_obstruction.v1', snappy=snappy.__version__, inputs={})
    Ns = {}
    for tag, p in (('old', a.old), ('new', a.new)):
        d, N, ndot = filled_exterior(p)
        Ns[tag] = N
        out['inputs'][tag] = dict(path=p, sha256=hashlib.sha256(open(p, 'rb').read()).hexdigest(),
                                  dotted_filled=ndot, cusps_after_filling=N.num_cusps(), H1=str(N.homology()))
    out['degrees'] = []
    out['status'] = 'NO_DIFFERENCE_FOUND'
    for n in [int(x) for x in a.degrees.split(',')]:
        po = Counter(str(C.homology()) for C in Ns['old'].covers(n))
        pn = Counter(str(C.homology()) for C in Ns['new'].covers(n))
        diff = {h: [po[h], pn[h]] for h in sorted(set(po) | set(pn)) if po[h] != pn[h]}
        out['degrees'].append(dict(degree=n, covers=[sum(po.values()), sum(pn.values())],
                                   differing_H1_counts_old_new=diff))
        print(f'{n}-fold covers: {sum(po.values())} vs {sum(pn.values())}; '
              f'{"multisets differ in %d classes" % len(diff) if diff else "same multiset"}', flush=True)
        if diff:
            out['status'] = 'DIFFERENT_EXTERIORS'
            break
    json.dump(out, open(a.output, 'w'), indent=1)
    print(out['status'])


if __name__ == '__main__':
    main()

"""Native Regina simplification of a closed 4-manifold triangulation, with checkpoints and S^4 recognition.

usage: closed_simplify.py closed.rga outdir rounds [known_standard_isosig ...]
Saves outdir/best.rga after every stage (and best.isosig once small).  Reports PL-standard S^4 if the result is
isomorphic to Regina's standard sphere (double of a pentachoron) or to a listed known-standard signature.
"""
import sys, os, time, regina
src, outdir, rounds = sys.argv[1], sys.argv[2], int(sys.argv[3])
known = set(sys.argv[4:])
os.makedirs(outdir, exist_ok=True)
T = regina.Triangulation4(regina.open(src))
S4 = regina.Example4.sphere()
known.add(S4.isoSig())
t00 = time.time()


def report(tag, t0):
    regina.PacketOfTriangulation4(T).save(os.path.join(outdir, 'best.rga'))
    std = False; sig = ''
    if T.size() <= 400:
        sig = T.isoSig(); open(os.path.join(outdir, 'best.isosig'), 'w').write(sig + '\n')
        std = sig in known or (T.size() == S4.size() and T.isIsomorphicTo(S4) is not None)
    print(f'{tag}: {T.size()} pentachora, {T.countVertices()} vertices, {time.time()-t0:.0f}s '
          f'(total {time.time()-t00:.0f}s) {sig if T.size() <= 40 else ""}{"  STANDARD_S4" if std else ""}', flush=True)
    if std:
        open(os.path.join(outdir, 'STANDARD.txt'), 'w').write(sig + '\n'); sys.exit(0)


t0 = time.time(); report('start', t0)
t0 = time.time(); T.simplify(); report('simplify', t0)
for r in range(rounds):
    t0 = time.time()
    U = regina.Triangulation4(T)
    U.simplifyUpDown(); U.simplify()
    if U.size() <= T.size():
        T = U
    report(f'round {r}', t0)

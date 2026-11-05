"""Time Regina's native simplification of the ideal Katie triangulation; save checkpoints."""
import sys, time, regina
sig = open(sys.argv[1]).read().split()[0]; out = sys.argv[2]
T = regina.Triangulation4.fromIsoSig(sig)
t = time.time(); nv = T.countVertices(); print('skeleton: vertices', nv, '%.1fs' % (time.time()-t), flush=True)
t = time.time(); T.simplify(); print('simplify ->', T.size(), 'pentachora, %d vertices, %.1fs' % (T.countVertices(), time.time()-t), flush=True)
open(out, 'w').write(T.isoSig() + '\n')
for r in range(int(sys.argv[3]) if len(sys.argv) > 3 else 0):
    t = time.time(); T.simplifyUpDown() if hasattr(T, 'simplifyUpDown') else None; T.simplify()
    print('round', r, '->', T.size(), '%.1fs' % (time.time()-t), flush=True)
    open(out, 'w').write(T.isoSig() + '\n')

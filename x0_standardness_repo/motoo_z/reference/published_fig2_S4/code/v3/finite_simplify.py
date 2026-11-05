"""Regina native simplification of a finite 4-manifold triangulation with boundary; report boundary size."""
import sys, time, regina
src, out = sys.argv[1], sys.argv[2]; rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 2
T = regina.Triangulation4(regina.open(src))
def rep(tag, t0):
    b = [x for x in T.boundaryComponents() if x.isReal()]
    print(f'{tag}: {T.size()} pentachora, {T.countVertices()} vertices, boundary tets {[x.size() for x in b]}, {time.time()-t0:.0f}s', flush=True)
t0 = time.time(); rep('start', t0)
t0 = time.time(); T.simplify(); rep('simplify', t0)
regina.PacketOfTriangulation4(T).save(out)
for r in range(rounds):
    t0 = time.time(); T.simplifyUpDown(); T.simplify(); rep(f'round {r}', t0)
    regina.PacketOfTriangulation4(T).save(out)

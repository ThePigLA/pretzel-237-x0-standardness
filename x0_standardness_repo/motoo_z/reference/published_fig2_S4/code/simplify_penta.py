#!/usr/bin/env python3
"""
Certified parallel simplification search for small closed 4-manifold triangulations (Regina 7.4.x).

Default start: the 8-pentachoron endpoint of the Figure-2-redraw reduction
    iLLLAQQkbcfdegegfghhhUagbdayaSazawaob8aca1atbba
Goal: reach the standard 2-pentachoron 4-sphere (or any triangulation listed in a --known file).

Once a triangulation with <= 4 pentachora is reached, a bounded bidirectional breadth-first search
connects it to Regina's standard S^4 (Example4.sphere(), the double of a pentachoron) by explicit moves;
segments found from the sphere's side are recorded as 'backward' steps (PL equivalence is symmetric).

Every accepted step is an explicit list of Regina moves applied to the canonically labelled
triangulation Triangulation4.fromIsoSig(parent).  Each move is recorded as (kind, face-dimension,
face-index) and performed with Regina's own checked move routines (pachner, move20, move44,
collapseEdge), which refuse any move that would change the PL type.  The `verify` sub-command
replays the whole chain from the start signature, so a reported result is a replayable certificate,
not just a smaller isomorphism signature.

Search: each episode starts from a state in a pool of the best triangulations found so far,
makes a random excursion upward (2-4 and 1-5 Pachner moves, mixed with neutral 3-3 and 4-4 moves),
then descends greedily with explicit reducing moves (edge collapses, 5-1, 4-2, 2-0), using a
bounded number of neutral 3-3 / 4-4 moves to escape plateaus.  Excursion size grows automatically
while the search is stuck.

Usage
  python simplify_penta.py search  --out runs/r1 [--hours 8] [--workers N] [--seed 1]
                                   [--start ISOSIG | --start-file T.rga] [--known known.txt]
  python simplify_penta.py search  --resume runs/r1 [--hours 8]
  python simplify_penta.py verify  runs/r1/CERTIFICATE.json [--known known.txt]
  python simplify_penta.py verify  runs/r1/BEST_CHAIN.json          (checks a partial chain)
  python simplify_penta.py exhaustive [--height 4] [--threads N] [--start ...]
        (Regina's breadth-first Pachner search: discovery only, no move path is recorded)

On a Mac:  caffeinate -i python -u simplify_penta.py search --out runs/r1 --hours 12
"""
import argparse, json, os, random, sys, time, hashlib
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

DEFAULT_START = 'iLLLAQQkbcfdegegfghhhUagbdayaSazawaob8aca1atbba'


# ----------------------------------------------------------------------------- moves
def _face(T, dim, idx):
    return (T.vertex, T.edge, T.triangle, T.tetrahedron, T.pentachoron)[dim](idx)


def _count(T, dim):
    return (T.countVertices, T.countEdges, T.countTriangles, T.countTetrahedra, T.size)[dim]()


def apply_move(T, mv):
    """Apply one recorded move; returns True iff Regina accepted and performed it."""
    kind, dim, idx = mv
    if idx >= _count(T, dim):
        return False
    f = _face(T, dim, idx)
    if kind == 'P':
        return bool(T.pachner(f))
    if kind == '20':
        return bool(T.move20(f))
    if kind == '44':
        return bool(T.move44(f))
    if kind == 'C':
        return bool(T.collapseEdge(f))
    raise ValueError(f'unknown move kind {kind}')


# reducing moves, in priority order: (kind, dim)
REDUCERS = [('C', 1), ('P', 0), ('P', 1), ('20', 0), ('20', 1), ('20', 2)]
NEUTRAL = [('P', 2), ('44', 1)]


def try_reduce_once(T, rng, moves):
    """Perform one reducing move if any exists; record it. Returns True if a move was made."""
    for kind, dim in REDUCERS:
        idxs = list(range(_count(T, dim)))
        rng.shuffle(idxs)
        for i in idxs:
            if apply_move(T, (kind, dim, i)):
                moves.append((kind, dim, i))
                return True
    return False


def descend(T, rng, moves, max_neutral):
    """Greedy explicit descent with bounded neutral moves. Size never increases."""
    neutral_left = max_neutral
    while True:
        if try_reduce_once(T, rng, moves):
            neutral_left = max_neutral
            continue
        if neutral_left <= 0:
            return
        # one random neutral move (3-3 about a triangle, or 4-4 about an edge)
        kind, dim = NEUTRAL[0] if rng.random() < 0.7 else NEUTRAL[1]
        n = _count(T, dim)
        for _ in range(12):
            i = rng.randrange(n)
            if apply_move(T, (kind, dim, i)):
                moves.append((kind, dim, i))
                break
        neutral_left -= 1


def excursion(T, rng, moves, k_up, p15):
    """Random upward excursion: k_up size-increasing moves interleaved with neutral moves."""
    for _ in range(k_up):
        if rng.random() < p15:
            kind, dim = 'P', 4                       # 1-5: +4 pentachora, +1 vertex
        else:
            kind, dim = 'P', 3                       # 2-4: +2 pentachora
        n = _count(T, dim)
        for _ in range(12):
            i = rng.randrange(n)
            if apply_move(T, (kind, dim, i)):
                moves.append((kind, dim, i)); break
        for _ in range(rng.randint(0, 3)):           # a few neutral shuffles
            kind, dim = NEUTRAL[0] if rng.random() < 0.7 else NEUTRAL[1]
            n = _count(T, dim)
            for _ in range(8):
                i = rng.randrange(n)
                if apply_move(T, (kind, dim, i)):
                    moves.append((kind, dim, i)); break


def episode(task):
    """Worker: one excursion + descent from a pool state. Returns the explicit move list."""
    import regina
    sig, seed, k_up, p15, max_neutral = task
    rng = random.Random(seed)
    T = regina.Triangulation4.fromIsoSig(sig)
    moves = []
    excursion(T, rng, moves, k_up, p15)
    descend(T, rng, moves, max_neutral)
    return dict(start=sig, moves=moves, end=T.isoSig(), size=T.size(), vertices=T.countVertices())


# ----------------------------------------------------------------------------- connector
ALL_MOVES = [('P', 0), ('P', 1), ('P', 2), ('P', 3), ('P', 4), ('20', 0), ('20', 1), ('20', 2), ('44', 1), ('C', 1)]


def neighbours(sig, max_size):
    import regina
    T = regina.Triangulation4.fromIsoSig(sig)
    out = []
    for kind, dim in ALL_MOVES:
        if kind == 'P' and dim == 4 and T.size()+4 > max_size:
            continue
        if kind == 'P' and dim == 3 and T.size()+2 > max_size:
            continue
        for i in range(_count(T, dim)):
            U = regina.Triangulation4(T)
            if apply_move(U, (kind, dim, i)):
                out.append(((kind, dim, i), U.isoSig()))
    return out


def connect(a, b, max_size, seconds):
    """Bidirectional BFS in the graph of Regina moves (sizes <= max_size).
    Returns a list of steps (each one move) from a to b, some recorded backwards, or None."""
    if a == b:
        return []
    par = [{a: None}, {b: None}]            # sig -> (parent sig, move applied to parent)
    front = [[a], [b]]
    t_end = time.time()+seconds
    while front[0] and front[1] and time.time() < t_end:
        side = 0 if len(front[0]) <= len(front[1]) else 1
        nxt = []
        for sig in front[side]:
            for mv, child in neighbours(sig, max_size):
                if child in par[side]:
                    continue
                par[side][child] = (sig, mv)
                if child in par[1-side]:
                    return _join_paths(par, child)
                nxt.append(child)
            if time.time() > t_end:
                break
        front[side] = nxt
    return None


def _join_paths(par, meet):
    steps = []
    cur = meet; fw = []
    while par[0][cur] is not None:          # a ... -> meet  (forward steps)
        p, mv = par[0][cur]; fw.append(dict(**{'from': p, 'to': cur, 'moves': [mv]})); cur = p
    steps = list(reversed(fw))
    cur = meet
    while par[1][cur] is not None:          # meet -> ... b : each recorded backwards
        p, mv = par[1][cur]                 # applying mv to p gives cur
        steps.append(dict(**{'from': cur, 'to': p, 'moves': [mv], 'direction': 'backward'}))
        cur = p
    return steps


# ----------------------------------------------------------------------------- bookkeeping
def chain_to(parents, start, sig):
    steps = []
    cur = sig
    while cur != start:
        p, mv = parents[cur]
        steps.append(dict(**{'from': p, 'moves': mv, 'to': cur}))
        cur = p
    return list(reversed(steps))


def save_json(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f)
    os.replace(tmp, path)


def load_known(path):
    if not path:
        return set()
    txt = open(path).read().split()
    return set(t.strip() for t in txt if t.strip())


def is_standard_s4(sig):
    import regina
    T = regina.Triangulation4.fromIsoSig(sig)
    S = regina.Example4.sphere()
    return T.size() == S.size() and T.isIsomorphicTo(S) is not None


# ----------------------------------------------------------------------------- search
def cmd_search(a):
    import regina
    os.makedirs(a.out or a.resume, exist_ok=True)
    out = a.out or a.resume
    state_path = os.path.join(out, 'STATE.json')
    log_path = os.path.join(out, 'log.txt')
    known = load_known(a.known)
    if a.resume:
        st = json.load(open(state_path))
        start = st['start']; parents = {k: tuple(v) for k, v in st['parents'].items()}
        pool = st['pool']; best = st['best_size']; seed = st['seed'] + 1000003
        k_up = st.get('k_up', 2); episodes = st.get('episodes', 0)
    else:
        if a.start_file:
            start = regina.Triangulation4(regina.open(a.start_file).firstChild()).isoSig() \
                if not a.start_file.endswith('.sig') else open(a.start_file).read().strip()
        else:
            start = a.start
        T0 = regina.Triangulation4.fromIsoSig(start)
        start = T0.isoSig()
        if not (T0.isValid() and T0.isClosed() and T0.isOrientable() and T0.isConnected()):
            sys.exit('start triangulation must be valid, closed, orientable and connected')
        parents = {}; pool = [start]; best = T0.size(); seed = a.seed; k_up = 2; episodes = 0
    rng = random.Random(seed)
    workers = a.workers or max(1, (os.cpu_count() or 2)-1)
    t_end = time.time()+3600*a.hours
    last_save = time.time(); last_improve_ep = episodes; t0 = time.time(); ep0 = episodes
    logf = open(log_path, 'a')

    def log(msg):
        line = time.strftime('%Y-%m-%d %H:%M:%S ')+msg
        print(line, flush=True); logf.write(line+'\n'); logf.flush()

    def checkpoint(final=False):
        save_json(state_path, dict(start=start, parents=parents, pool=pool, best_size=best, seed=seed,
                                   k_up=k_up, episodes=episodes, regina=regina.versionString()))
        bs = pool[0]
        save_json(os.path.join(out, 'BEST_CHAIN.json'),
                  dict(start=start, end=bs, end_size=best, steps=chain_to(parents, start, bs),
                       regina=regina.versionString(), complete=False))

    log(f'start {start} size {regina.Triangulation4.fromIsoSig(start).size()}; workers {workers}; '
        f'{"resumed" if a.resume else "new run"}; best {best}; pool {len(pool)}')
    done = False
    tried_connect = set()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        inflight = set()
        while not done and time.time() < t_end:
            while len(inflight) < 2*workers:
                sig = rng.choice(pool)
                p15 = 0.15 if k_up < 6 else 0.3
                task = (sig, rng.getrandbits(63), rng.randint(1, k_up), p15, a.neutral)
                inflight.add(ex.submit(episode, task))
            fin, inflight = wait(inflight, timeout=5, return_when=FIRST_COMPLETED)
            for fu in fin:
                r = fu.result(); episodes += 1
                sig, sz = r['end'], r['size']
                if sig == r['start'] or sig in parents or sig == start:
                    continue
                if sz < best:
                    parents[sig] = (r['start'], r['moves'])
                    best = sz; pool = [sig]; last_improve_ep = episodes; k_up = 2
                    log(f'NEW BEST size {sz} (vertices {r["vertices"]}) after {episodes} episodes: {sig}')
                    checkpoint()
                elif sz == best and len(pool) < a.pool:
                    parents[sig] = (r['start'], r['moves']); pool.append(sig)
                target = None; tail = []
                if sig in known:
                    target = 'triangulation listed in --known'
                elif sz <= a.connect_at and sig not in tried_connect:
                    tried_connect.add(sig)
                    std = regina.Example4.sphere().isoSig()
                    log(f'size {sz} reached ({sig}); trying to connect to the standard S^4 {std} '
                        f'(moves within size <= {sz+a.connect_extra}, {a.connect_seconds}s)')
                    tail = connect(sig, std, sz+a.connect_extra, a.connect_seconds)
                    if tail is not None:
                        target = 'standard S^4 (Regina Example4.sphere(), the double of a pentachoron)'
                    else:
                        log('   not connected within the bounds; continuing the search')
                if target:
                    if sig not in parents:
                        parents[sig] = (r['start'], r['moves'])
                    steps = chain_to(parents, start, sig)+tail
                    sig = steps[-1]['to'] if tail else sig
                    save_json(os.path.join(out, 'CERTIFICATE.json'),
                              dict(start=start, end=sig, end_size=regina.Triangulation4.fromIsoSig(sig).size(),
                                   target=target, steps=steps,
                                   total_moves=sum(len(s['moves']) for s in steps),
                                   regina=regina.versionString(), complete=True))
                    log(f'*** REACHED {target}: {sig}. Certificate written; run the verify command. ***')
                    done = True; break
            # adapt excursion size while stuck
            if episodes-last_improve_ep > a.stuck and k_up < a.kmax:
                k_up += 1; last_improve_ep = episodes
                log(f'stuck at size {best}: raising excursion size to {k_up} (pool {len(pool)})')
            if time.time()-last_save > 60:
                rate = (episodes-ep0)/max(time.time()-t0, 1e-9)
                log(f'episodes {episodes} ({rate:.0f}/s), best {best}, pool {len(pool)}, excursion<= {k_up}')
                checkpoint(); last_save = time.time()
        for fu in inflight:
            fu.cancel()
    checkpoint(final=True)
    log(f'finished: best size {best}, pool {len(pool)}, episodes {episodes}. Resume with --resume {out}')


# ----------------------------------------------------------------------------- verify
def cmd_verify(a):
    import regina
    c = json.load(open(a.cert))
    known = load_known(a.known)
    start = c['start']
    T = regina.Triangulation4.fromIsoSig(start)
    if T.isoSig() != start:
        sys.exit('FAIL: start signature is not canonical')
    print(f'Regina {regina.versionString()} (certificate made with {c.get("regina")})')
    print(f'start: size {T.size()}, f-vector {list(T.fVector())}, valid {T.isValid()}, closed {T.isClosed()}, '
          f'orientable {T.isOrientable()}')
    if a.expect_start and start != a.expect_start:
        sys.exit(f'FAIL: start is not the expected triangulation {a.expect_start}')
    cur = start; nmoves = 0
    for k, st in enumerate(c['steps']):
        if st['from'] != cur:
            sys.exit(f'FAIL: step {k} does not continue the chain')
        back = st.get('direction') == 'backward'
        src, dst = (st['to'], st['from']) if back else (st['from'], st['to'])
        T = regina.Triangulation4.fromIsoSig(src)
        for j, mv in enumerate(st['moves']):
            if not apply_move(T, tuple(mv)):
                sys.exit(f'FAIL: step {k} move {j} {mv} was refused by Regina')
            nmoves += 1
        if T.isoSig() != dst:
            sys.exit(f'FAIL: step {k} ends at {T.isoSig()} instead of {dst}')
        cur = st['to']
    if cur != c['end']:
        sys.exit('FAIL: chain does not end at the recorded endpoint')
    E = regina.Triangulation4.fromIsoSig(cur)
    print(f'replayed {len(c["steps"])} steps, {nmoves} moves; endpoint size {E.size()}, '
          f'f-vector {list(E.fVector())}: {cur}')
    if is_standard_s4(cur):
        print('PASS_STANDARD_S4: the start triangulation is PL-homeomorphic to the standard 4-sphere '
              '(via the replayed Regina moves).')
    elif cur in known:
        print('PASS_KNOWN_TARGET: the chain reaches a triangulation listed in the --known file.')
    else:
        print('PASS_PARTIAL_PATH: the recorded chain is valid; the endpoint is not a recognised target.')


# ----------------------------------------------------------------------------- exhaustive
def cmd_exhaustive(a):
    import regina
    T = regina.Triangulation4.fromIsoSig(a.start)
    n0 = T.size()
    threads = a.threads or (os.cpu_count() or 2)
    print(f'exhaustive Pachner search from size {n0}, height {a.height}, {threads} threads '
          f'(discovery only: no move path is recorded)', flush=True)
    found = []
    def action(sig, tri):
        if tri.size() < n0:
            found.append(sig); return True
        return False
    t = time.time()
    T.retriangulate(a.height, threads, action)
    if found:
        print(f'FOUND smaller triangulation (size {regina.Triangulation4.fromIsoSig(found[0]).size()}) '
              f'in {time.time()-t:.0f}s: {found[0]}')
        print('Seed a certified search from it with:  search --start <that signature> ; note that the link '
              'back to the original start is then only "found by Regina retriangulate", not a recorded path.')
    else:
        print(f'no smaller triangulation within height {a.height} ({time.time()-t:.0f}s)')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('search')
    s.add_argument('--out'); s.add_argument('--resume')
    s.add_argument('--start', default=DEFAULT_START); s.add_argument('--start-file')
    s.add_argument('--hours', type=float, default=8.0); s.add_argument('--workers', type=int, default=0)
    s.add_argument('--seed', type=int, default=1); s.add_argument('--known')
    s.add_argument('--neutral', type=int, default=40, help='neutral moves allowed per plateau in descent')
    s.add_argument('--pool', type=int, default=400, help='max distinct states kept at the best size')
    s.add_argument('--stuck', type=int, default=20000, help='episodes without improvement before enlarging excursions')
    s.add_argument('--kmax', type=int, default=14, help='largest excursion (number of 2-4/1-5 moves)')
    s.add_argument('--connect-at', type=int, default=4, help='size at which to try connecting to the standard S^4')
    s.add_argument('--connect-extra', type=int, default=6, help='connector may go this many pentachora above the start')
    s.add_argument('--connect-seconds', type=int, default=900)
    v = sub.add_parser('verify'); v.add_argument('cert'); v.add_argument('--known')
    v.add_argument('--expect-start', default=None)
    e = sub.add_parser('exhaustive'); e.add_argument('--start', default=DEFAULT_START)
    e.add_argument('--height', type=int, default=4); e.add_argument('--threads', type=int, default=0)
    a = ap.parse_args()
    if a.cmd == 'search':
        if not (a.out or a.resume):
            sys.exit('give --out DIR for a new run or --resume DIR')
        cmd_search(a)
    elif a.cmd == 'verify':
        cmd_verify(a)
    else:
        cmd_exhaustive(a)


if __name__ == '__main__':
    main()

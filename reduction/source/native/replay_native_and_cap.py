#!/usr/bin/env python3
"""Independent, bounded replay of X0's native macros and an explicit PL capping.

No reviewer code is imported. The native simplify macros remain a trust
dependency on the recorded Regina engine. Capping is constructed from an
explicit boundary triangulation isomorphism, whose face-pairing equations and
quotient embeddings are checked, rather than from an Euler-characteristic test.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import importlib.metadata
import itertools
import json
import platform
from pathlib import Path
import resource
import signal
import sys
import time
import warnings

warnings.filterwarnings('ignore', category=DeprecationWarning)
import regina

G6_SIG = 'gAMLMQaabcdeffff2a+a2avbmaoaxaOb0aUa'
T5_SIG = 'fHMLQcabcdeeeecanbpaoazaobWbGb'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    return sha(Path(path).read_bytes())


def integer(x):
    return type(x) is int


def numbered_rows(t):
    return [[(p.adjacentPentachoron(f).index(), str(p.adjacentGluing(f)))
             if p.adjacentPentachoron(f) is not None else None
             for f in range(5)] for p in t.pentachora()]


def numbered_hash(t):
    return sha(json.dumps(numbered_rows(t)).encode())


def verify_iso(a, b, iso, dim):
    """Check the permutation-conjugacy equations on every paired half-facet."""
    need(iso is not None, 'isomorphism absent')
    need(a.size() == b.size(), 'isomorphism sizes differ')
    images = [iso.simpImage(i) for i in range(a.size())]
    need(sorted(images) == list(range(b.size())), 'simplex map not bijective')
    maps = [iso.facetPerm(i) for i in range(a.size())]
    equations = 0
    for i, p in enumerate(maps):
        need(sorted(p[k] for k in range(dim + 1)) == list(range(dim + 1)),
             'vertex map not a permutation')
        s, u = a.simplex(i), b.simplex(images[i])
        for f in range(dim + 1):
            t, v = s.adjacentSimplex(f), u.adjacentSimplex(p[f])
            need((t is None) == (v is None), 'boundary pairing mismatch')
            if t is None:
                continue
            j = t.index()
            need(v.index() == images[j], 'adjacent simplex mismatch')
            g, h = s.adjacentGluing(f), u.adjacentGluing(p[f])
            need(all(maps[j][g[k]] == h[p[k]] for k in range(dim + 1)),
                 'face-pairing conjugacy equation failed')
            equations += 1
    return {'simplex_images': images,
            'vertex_permutations': [[p[k] for k in range(dim + 1)] for p in maps],
            'paired_half_facet_equations': equations}


def local_face_key(tri, dim, simplex, ordered_vertices):
    """A global face index plus its ordered affine vertex frame."""
    k = len(ordered_vertices) - 1
    if k == dim:
        return (simplex, tuple(ordered_vertices))
    vertices = tuple(ordered_vertices)
    perm = getattr(regina, f'Perm{dim + 1}')(
        *(vertices + tuple(v for v in range(dim + 1) if v not in vertices)))
    local = getattr(regina, f'Face{dim}_{k}').faceNumber(perm)
    s = tri.simplex(simplex)
    frame = s.faceMapping(k, local)
    inverse = {frame[i]: i for i in range(dim + 1)}
    return (s.face(k, local).index(), tuple(inverse[v] for v in vertices))


def boundary_embedding_check(tri):
    """Check exact point identifications of every ordered boundary subface.

    Ordered affine frames distinguish identifications with different vertex
    permutations. Equality in the boundary quotient must agree in both
    directions with equality in the ambient triangulation.
    """
    bc = tri.boundaryComponent(0)
    b = bc.build()
    rows, counts = [], []
    for k in range(4):
        forward, inverse = {}, {}
        checks = 0
        for i in range(b.size()):
            e = bc.tetrahedron(i).embedding(0)
            emb = e.vertices()
            for ordered in itertools.permutations(range(4), k + 1):
                bk = local_face_key(b, 3, i, ordered)
                ak = local_face_key(tri, 4, e.simplex().index(),
                                    tuple(emb[x] for x in ordered))
                need(bk not in forward or forward[bk] == ak,
                     'boundary embedding not well-defined')
                need(ak not in inverse or inverse[ak] == bk,
                     'boundary embedding is not injective')
                forward[bk], inverse[ak] = ak, bk
                checks += 1
        counts.append(len({q[0] for q in forward}))
        rows.append({'dimension': k, 'ordered_local_frames_checked': checks,
                     'quotient_ordered_frames': len(forward),
                     'ambient_face_count': len({q[0] for q in inverse})})
    need(counts == list(b.fVector()), 'boundary face coverage mismatch')
    return {'f_vector': counts, 'exact_affine_quotient_embedding': True,
            'by_dimension': rows}


def native_export(t):
    gluings = []
    for i in range(t.size()):
        p = t.simplex(i)
        for f in range(5):
            q = p.adjacentSimplex(f)
            if q is None:
                continue
            perm = p.adjacentGluing(f)
            if (i, f) <= (q.index(), perm[f]):
                gluings.append([i, f, q.index(), [perm[k] for k in range(5)]])
    return {'simplex_count': t.size(), 'gluings': gluings,
            'isoSig': t.isoSig(), 'f_vector': list(t.fVector()),
            'numbered_sha256': numbered_hash(t), 'numbered_rows': numbered_rows(t)}


def source_triangulation(data):
    need(set(data) == {'simplex_count', 'facets', 'gluings'}, 'source schema')
    facets = data['facets']
    need(integer(data['simplex_count']) and data['simplex_count'] == len(facets) == 1332,
         'wrong source facet count')
    need(all(type(f) is list and len(f) == 5 and all(integer(v) for v in f)
             and f == sorted(set(f)) for f in facets), 'invalid abstract facets')
    need(len(set(map(tuple, facets))) == len(facets), 'duplicate source facet')
    faces = defaultdict(list)
    for i, f in enumerate(facets):
        for k in range(5):
            faces[tuple(f[:k] + f[k + 1:])].append((i, k))
    need(all(len(x) == 2 for x in faces.values()), 'source is not closed pseudomanifold')
    expected, actual, used = {}, {}, set()
    for side in faces.values():
        (i, k), (j, l) = sorted(side)
        p = [l if q == k else facets[j].index(v) for q, v in enumerate(facets[i])]
        expected[(i, k)] = (j, tuple(p))
        inv = [p.index(q) for q in range(5)]
        expected[(j, l)] = (i, tuple(inv))
    for record in data['gluings']:
        need(type(record) is list and len(record) == 4, 'gluing schema')
        i, f, j, p = record
        need(all(integer(x) for x in [i, f, j]) and 0 <= i < len(facets)
             and 0 <= j < len(facets) and 0 <= f < 5 and type(p) is list
             and all(integer(x) for x in p) and sorted(p) == list(range(5)),
             'invalid gluing indices or permutation')
        a, b = (i, f), (j, p[f])
        need(a != b and a not in used and b not in used, 'duplicate/self-paired side')
        used.update((a, b))
        actual[a] = (j, tuple(p))
        actual[b] = (i, tuple(p.index(q) for q in range(5)))
    need(actual == expected, 'native gluing differs from abstract source')
    t = regina.Triangulation4.fromGluings(
        len(facets), [(i, f, j, regina.Perm5(*p)) for i, f, j, p in data['gluings']])
    need(t.isValid() and t.isConnected() and t.isClosed() and t.isOrientable()
         and not t.isIdeal(), 'native source admission failed')
    fvec = [len({tuple(q) for f in facets for q in itertools.combinations(f, k)})
            for k in range(1, 6)]
    need(fvec == list(t.fVector()), 'native and abstract face count mismatch')
    return t, {'facets': len(facets), 'vertices': fvec[0], 'f_vector': fvec,
               'all_half_facets_checked': len(actual), 'gluing_pairs': len(data['gluings']),
               'facet_sha256': sha(json.dumps(facets, separators=(',', ':')).encode())}


def face_pushout(t, c, u):
    out = []
    for k in range(5):
        images = []
        for tri, offset in [(t, 0), (c, t.size())]:
            im = []
            for f in (tri.simplices() if k == 4 else tri.faces(k)):
                if k == 4:
                    im.append(f.index() + offset)
                else:
                    emb = f.embedding(0)
                    q = emb.vertices()
                    im.append(local_face_key(u, 4, emb.simplex().index() + offset,
                                             tuple(q[i] for i in range(k + 1)))[0])
            need(len(set(im)) == (tri.size() if k == 4 else tri.countFaces(k)),
                 'capping collapses a face within a side')
            images.append(set(im))
        common = len(images[0] & images[1])
        expected_common = t.boundaryComponent(0).build().countFaces(k) if k < 4 else 0
        need(common == expected_common, 'capping side intersection is not just boundary')
        union_faces = u.size() if k == 4 else u.countFaces(k)
        need(len(images[0] | images[1]) == union_faces, 'capping misses faces')
        out.append({'dimension': k, 'left_faces': len(images[0]),
                    'cone_faces': len(images[1]), 'intersection': common,
                    'union_faces': union_faces, 'both_side_maps_injective': True})
    return out


def cap_by_explicit_boundary_map(t):
    need(t.size() == 5 and t.isoSig() == T5_SIG, 'unexpected native endpoint')
    need(t.isValid() and t.isConnected() and t.isOrientable() and not t.isIdeal()
         and t.countBoundaryComponents() == 1 and t.countBoundaryTetrahedra() == 1,
         'endpoint is not admitted manifold with one real boundary')
    b = t.boundaryComponent(0).build()
    ref = regina.Example3.threeSphere()
    sphere_iso = verify_iso(b, ref, b.isIsomorphicTo(ref), 3)
    need(b.isoSig() == ref.isoSig() == 'bkaagj', 'unexpected boundary sphere')
    t_embedding = boundary_embedding_check(t)
    c = regina.Example4.singleCone(b)
    need(c.size() == 1 and c.isValid() and not c.isIdeal()
         and c.countBoundaryComponents() == 1 and c.countBoundaryTetrahedra() == 1,
         'cone over standard sphere is not admitted compact ball')
    cb = c.boundaryComponent(0).build()
    iso = b.isIsomorphicTo(cb)
    boundary_iso = verify_iso(b, cb, iso, 3)
    c_embedding = boundary_embedding_check(c)
    u = regina.Triangulation4(t)
    u.insertTriangulation(c)
    maps = []
    for i in range(b.size()):
        ei = t.boundaryComponent(0).tetrahedron(i).embedding(0)
        ej = c.boundaryComponent(0).tetrahedron(iso.simpImage(i)).embedding(0)
        a, z, p = ei.vertices(), ej.vertices(), iso.facetPerm(i)
        perm = [None] * 5
        for k in range(4):
            perm[a[k]] = z[p[k]]
        perm[a[4]] = z[4]
        u.simplex(ei.simplex().index()).join(
            a[4], u.simplex(t.size() + ej.simplex().index()), regina.Perm5(*perm))
        maps.append({'left_pentachoron': ei.simplex().index(), 'left_facet': a[4],
                     'cone_pentachoron': ej.simplex().index(), 'cone_facet': z[4],
                     'perm5': perm})
    need(u.isValid() and u.isConnected() and u.isClosed() and u.isOrientable()
         and not u.isIdeal(), 'closed capping admission failed')
    pushout = face_pushout(t, c, u)
    expected = [a + z - q for a, z, q in zip(t.fVector(), c.fVector(), list(b.fVector()) + [0])]
    need(list(u.fVector()) == expected == [2, 3, 12, 15, 6], 'f-vector pushout equation')
    g = regina.Triangulation4.fromIsoSig(G6_SIG)
    target_iso = verify_iso(u, g, u.isIsomorphicTo(g), 4)
    need(u.size() == 6 and u.isoSig() == G6_SIG, 'wrong closed endpoint')
    return u, {'status': 'PASS', 'no_simplification_performed_after_capping': True,
               'boundary_isoSig': b.isoSig(), 'boundary_reference': 'regina.Example3.threeSphere()',
               'boundary_to_reference_isomorphism': sphere_iso,
               'left_boundary_embedding': t_embedding, 'cone_boundary_embedding': c_embedding,
               'boundary_identification_isomorphism': boundary_iso, 'attaching_facet_maps': maps,
               'face_pushout': pushout, 'left_f_vector': list(t.fVector()),
               'boundary_f_vector': list(b.fVector()), 'cone_f_vector': list(c.fVector()),
               'closed_f_vector': list(u.fVector()), 'closed_chi': u.eulerCharTri(),
               'target_isoSig': G6_SIG, 'target_isomorphism': target_iso}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--certificate', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path,
                        help='new or empty directory for replayed endpoint and receipts')
    args = parser.parse_args()
    need(sys.flags.optimize == 0, 'optimized Python forbidden')
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (900, 900))
    signal.alarm(1200)
    need(not args.output.exists() or not any(args.output.iterdir()), 'output must be new or empty')
    args.output.mkdir(parents=True, exist_ok=True)
    start_wall, start_cpu = time.monotonic(), time.process_time()
    raw, cert_raw = args.source.read_bytes(), args.certificate.read_bytes()
    d, cert = json.loads(gzip.decompress(raw)), json.loads(cert_raw)
    need(set(cert) == {'source', 'source_sha256', 'puncture', 'regina_version', 'blocks',
                       'final_size', 'final_isoSig', 'ball', 'producer'}, 'certificate schema')
    need(sha(raw) == cert['source_sha256'], 'source hash mismatch')
    need(cert['final_size'] == 5 and cert['final_isoSig'] == T5_SIG and cert['ball'] is False,
         'certificate does not describe the expected five-pentachoron candidate')
    need(integer(cert['puncture']) and cert['puncture'] == 25, 'puncture index mismatch')
    need(type(cert['blocks']) is list and len(cert['blocks']) == 5173, 'macro count mismatch')
    for k, block in enumerate(cert['blocks']):
        need(type(block) is dict and set(block) == {'kicks', 'hash'}
             and type(block['kicks']) is list and all(integer(x) and x >= 0 for x in block['kicks'])
             and type(block['hash']) is str and len(block['hash']) == 64
             and all(x in '0123456789abcdef' for x in block['hash']), f'block {k} schema')
    need(cert['blocks'][0]['kicks'] == [], 'first macro unexpectedly has kicks')
    t, source = source_triangulation(d)
    source['file_sha256'], source['native_numbered_sha256'] = sha(raw), numbered_hash(t)
    p = t.simplex(cert['puncture'])
    unique_faces = [len({p.face(k, j).index() for j in range([5, 10, 10, 5][k])})
                    for k in range(4)]
    need(unique_faces == [5, 10, 10, 5] and all(p.adjacentSimplex(f) != p for f in range(5)),
         'puncture is not an embedded four-simplex')
    puncture = {'index': cert['puncture'], 'vertices': d['facets'][cert['puncture']],
                'distinct_faces_by_dimension': unique_faces, 'embedded': True}
    t.removePentachoron(p)
    puncture['punctured_numbered_sha256'] = numbered_hash(t)
    puncture['punctured_size'] = t.size()
    need(t.size() == 1331 and t.countBoundaryComponents() == 1
         and t.countBoundaryTetrahedra() == 5, 'puncture did not create simplex sphere boundary')
    sizes, kick_count, explicit_max = [], 0, t.size()
    last_report = time.monotonic()
    for k, block in enumerate(cert['blocks']):
        for index in block['kicks']:
            need(index < t.countTetrahedra(), f'kick index out of range at block {k}')
            tetrahedron = t.tetrahedron(index)
            need(not tetrahedron.isBoundary(), f'boundary kick at block {k}')
            before = t.size()
            need(t.pachner(tetrahedron, True, True), f'illegal guarded kick at block {k}')
            need(t.size() == before + 2, f'kick was not a 2-4 move at block {k}')
            kick_count += 1
            explicit_max = max(explicit_max, t.size())
        regina.RandomEngine.reseedWithDefault()
        t.simplify()
        need(numbered_hash(t) == block['hash'], f'numbered gluing hash mismatch at block {k}')
        sizes.append(t.size())
        if time.monotonic() - last_report > 25:
            print(json.dumps({'blocks_replayed': k + 1, 'pentachora': t.size(),
                              'guarded_kicks': kick_count,
                              'wall_seconds': time.monotonic() - start_wall}), flush=True)
            last_report = time.monotonic()
    need(t.size() == cert['final_size'] and t.isoSig() == cert['final_isoSig'], 'wrong macro endpoint')
    native_complete_wall, native_complete_cpu = time.monotonic(), time.process_time()
    u, cap = cap_by_explicit_boundary_map(t)
    outputs = {'replayed_T5.json': native_export(t), 'capped_G6.json': native_export(u),
               'macro_block_sizes.json': {'sizes_after_blocks': sizes}, 'capping_receipt.json': cap}
    for name, obj in outputs.items():
        (args.output / name).write_text(json.dumps(obj, indent=2) + '\n')
    engine = Path(regina.engine.__file__)
    receipt = {'schema': 'x0-independent-native-capping-audit-v1', 'status': 'PASS',
               'claim': 'The 1332-facet source is PL homeomorphic to G6, conditional on Regina move implementations.',
               'does_not_prove_G6_standard': True, 'source': source,
               'certificate_file_sha256': sha(cert_raw), 'verifier_file_sha256': file_sha(__file__),
               'puncture': puncture, 'macros': {'blocks_replayed': len(sizes),
                  'numbered_hashes_checked': len(sizes), 'guarded_2_4_moves': kick_count,
                  'first_simplify_size': sizes[0], 'final_size': t.size(), 'final_isoSig': t.isoSig(),
                  'final_numbered_sha256': numbered_hash(t), 'minimum_recorded_size': min(sizes),
                  'maximum_recorded_size': max(sizes), 'maximum_explicit_size_including_source': explicit_max,
                  'macro_internal_maximum_not_observed': True},
               'capping': cap, 'environment': {'python': platform.python_version(),
                  'platform': platform.platform(), 'machine': platform.machine(),
                  'regina_api_version': regina.versionString(),
                  'regina_distribution_version': importlib.metadata.version('regina'),
                  'certificate_regina_api_version': cert['regina_version'],
                  'engine_filename': engine.name, 'engine_sha256': file_sha(engine),
                  'python_optimization': sys.flags.optimize},
               'resources': {'cpu_seconds': time.process_time() - start_cpu,
                  'wall_seconds': time.monotonic() - start_wall,
                  'through_macros_cpu_seconds': native_complete_cpu - start_cpu,
                  'through_macros_wall_seconds': native_complete_wall - start_wall,
                  'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                  'address_space_limit_bytes': 8 * 1024**3, 'cpu_limit_seconds': 900,
                  'wall_alarm_seconds': 1200},
               'outputs': {name: {'sha256': file_sha(args.output / name),
                                  'bytes': (args.output / name).stat().st_size} for name in outputs},
               'trust_dependencies': ['Regina guarded Pachner moves and simplify implementations',
                                     'Regina skeleton/face mappings and Example3.threeSphere reference'],
               'input_files_unchanged': file_sha(args.source) == sha(raw)
                                       and file_sha(args.certificate) == sha(cert_raw)}
    need(receipt['input_files_unchanged'], 'inputs mutated during replay')
    (args.output / 'independent_replay.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'blocks': len(sizes), 'kicks': kick_count,
                      'final_T5': t.isoSig(), 'closed_G6': u.isoSig(),
                      'resources': receipt['resources']}), flush=True)


if __name__ == '__main__':
    main()

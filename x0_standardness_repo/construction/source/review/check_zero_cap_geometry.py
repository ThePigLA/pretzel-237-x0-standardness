#!/usr/bin/env python3
"""Independent m=0 full-cell identity and cap cone/product checks.

No producer or previous verifier is imported. This only reads frozen data.
The optional candidate comparison reads an existing closed model; it never
assembles one or changes any model. Standard-library Python only.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import Counter
import argparse
import gzip
import hashlib
import json
import time


def require(x, why):
    if not x:
        raise ValueError(why)


def load(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def tupled(value):
    return tuple(map(tupled, value)) if isinstance(value, list) else value


def mod1(x):
    return x - x.numerator // x.denominator


def determinant(rows):
    a = [list(map(Q, r)) for r in rows]
    result = Q(1)
    for j in range(len(a)):
        pivot = next((k for k in range(j, len(a)) if a[k][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            result = -result
        q = a[j][j]
        result *= q
        for k in range(j+1, len(a)):
            scale = a[k][j] / q
            for c in range(j, len(a)):
                a[k][c] -= scale * a[j][c]
    return result


def diamond(theta):
    q = 4*mod1(theta)
    side = q.numerator // q.denominator
    a = q-side
    corners = ((Q(1),Q(0)), (Q(0),Q(1)), (Q(-1),Q(0)), (Q(0),Q(-1)))
    return tuple((1-a)*corners[side][k] + a*corners[(side+1)%4][k] for k in (0,1))


def physical_tetrahedron(data, f, times):
    tags = data['vertices']
    centers = [tags[v] for v in f if tags[v][0] == 'C']
    require(len(centers) == 1, 'tetrahedron must have one hemisphere/slab center')
    _, hemisphere, slab = centers[0]
    points = []
    for v in f:
        tag = tags[v]
        if tag[0] == 'E':
            xy, t = diamond(Q(tag[1])), Q(tag[2])
        elif tag[0] == 'P':
            require(tag[1] == hemisphere, 'wrong pole hemisphere')
            xy, t = (Q(0),Q(0)), Q(tag[2], times)
        elif tag[0] == 'C':
            require(tag[1:] == [hemisphere, slab], 'wrong slab center')
            xy, t = (Q(0),Q(0)), Q(2*slab+1,2*times)
        else:
            raise ValueError('unexpected spherical boundary vertex')
        if slab == times-1 and t == 0:
            t = Q(1)
        require(Q(slab,times) <= t <= Q(slab+1,times), 'vertex lies outside slab')
        points.append((*xy,t))
    require(determinant([(*p,Q(1)) for p in points]) != 0, 'degenerate spherical tetrahedron')
    return hemisphere, slab, points


def faces(fs):
    out = set()
    for f in fs:
        for k in range(1,len(f)+1):
            out.update(combinations(f,k))
    return out


def zero_cap_check(root):
    paths = {k: root/name for k,name in {
        'shared':'cap/data/shared_boundary.json.gz',
        'm0':'cap/data/cap_m0.json.gz',
        'm2':'cap/data/cap_m2.json.gz',
        'central_certificate':'central_roll/certificate.json.gz',
        'central_receipt':'central_roll/independent_verification.json',
        'bridge':'cap/data/source_cap_bridge.json.gz',
    }.items()}
    b,c0,c2 = (load(paths[k]) for k in ('shared','m0','m2'))
    require(c0['twist'] == 0 and c2['twist'] == 2, 'cap parameters')
    require(c0['vertices'] == c2['vertices'] and c0['deck'] == c2['deck'], 'cap label/deck identity')
    require(c0['boundary'] == c2['boundary'] == b['facets'], 'literal shared boundary')
    require(c0['shared_boundary_sha256'] == c2['shared_boundary_sha256'] == digest(b), 'shared boundary binding')
    require(c0['beta']['cap_full_turns'] == c0['beta']['downstairs_winding'] == 0,
            'zero full lift, not only parity')
    require(c0['beta']['normalized_displacement'] == [0,1] and
            c0['beta']['unwrapped_displacement'] == [0,4], 'zero lifted path')
    coarse = b['coarse_source_boundary']
    rows = c0['boundary_to_source_refinement']
    require(len(rows) == len(b['facets']) == 1792, 'complete boundary chart list')
    identity_rows = []
    for i,(f,row) in enumerate(zip(b['facets'],rows)):
        target_hem,target_slab,target = physical_tetrahedron(b,f,b['times'])
        parent = row['parent']
        require(type(parent) is int and 0 <= parent < len(coarse['facets']), 'parent range')
        source_hem,source_slab,parent_points = physical_tetrahedron(coarse,coarse['facets'][parent],b['times'])
        require((target_hem,target_slab) == (source_hem,source_slab), 'hemisphere/slab mismatch')
        matrix = [list(map(Q,r)) for r in row['barycentric']]
        require(len(matrix) == 4 and all(len(r) == 4 and min(r) >= 0 and sum(r) == 1 for r in matrix),
                'barycentric simplex map')
        require(determinant(matrix) != 0, 'singular inverse chart')
        source = [tuple(sum(r[j]*parent_points[j][k] for j in range(4)) for k in range(3)) for r in matrix]
        require(source == target, 'm0 spherical chart is not pointwise identity')
        identity_rows.append([i,parent,source_hem,source_slab])
    # The four noncoplanar vertex equations determine the entire affine map,
    # so this checks the interiors of every hemisphere tetrahedron as well.
    tags = list(map(tupled,c0['vertices']))
    index = {v:i for i,v in enumerate(tags)}
    filling = set(map(tuple,c0['filling_torus']))
    collar = {tuple(sorted(r['vertices'])) for r in c0['collar_charts']}
    require(len(collar) == len(c0['collar_charts']) and collar <= filling, 'collar inclusion')
    n,t = c0['angles'],c0['times']
    def v(i,j):
        return index[('H',i%n)] if j == -1 else index[('I',i%n,j%t)]
    product = set()
    for i in range(n):
        for j in range(t):
            for path in (
                ((i,-1),(i+1,-1),(i+1,j),(i+1,j+1)),
                ((i,-1),(i,j),(i+1,j),(i+1,j+1)),
                ((i,-1),(i,j),(i,j+1),(i+1,j+1)),
            ):
                product.add(tuple(sorted(v(a,k) for a,k in path)))
    require(len(product) == 96 and not product & collar and product | collar == filling,
            'filling torus is not the expected circle times cone-circle plus its refinement collar')
    cap_facets = set(map(tuple,c0['facets']))
    halves, bases = [],[]
    for p in (0,1):
        base = set(map(tuple,b['hemispheres'][str(p)])) | filling
        require(not (set(map(tuple,b['hemispheres'][str(p)])) & filling), 'hemisphere/filling 3-facet overlap')
        cone_vertex = index[('A',p)]
        expected = {tuple(sorted((*f,cone_vertex))) for f in base}
        require(expected == set(map(tuple,c0['hemisphere_cones'][str(p)])), 'cap half is not the asserted cone')
        require(all(cone_vertex not in f for f in base), 'cone apex lies in its base')
        # Every triangle of this 3D cone base is paired. Its geometric type
        # follows from the checked D2 x S1 and S1 x D2 attachment by identity.
        ridges = Counter(tr for tet in base for tr in combinations(tet,3))
        require(set(ridges.values()) == {2}, 'cone base is not closed')
        halves.append(expected)
        bases.append(base)
    require(halves[0].isdisjoint(halves[1]) and halves[0]|halves[1] == cap_facets,
            'cap is not exactly its two cone halves')
    require(faces(bases[0]) & faces(bases[1]) == faces(filling), 'cone halves meet outside filling torus')
    f2 = set(map(tuple,c2['facets']))
    require(len(cap_facets) == len(f2) == 5824 and len(cap_facets-f2) == len(f2-cap_facets) == 3008,
            'unexpected cap delta')
    central = load(paths['central_certificate'])
    receipt = load(paths['central_receipt'])
    require(receipt['passed'] and receipt['certificate_sha256'] == digest(central), 'central receipt binding')
    require(central['source_profile_transfer']['map'] == '(u,v,r,t)->(u+J(r,t),v,r,t)' and
            central['time_seam'] == 'identity on N and E products and on straightened central coordinates; no meridian shift',
            'central transfer is not u-only')
    require(receipt['straightening']['full_integer_winding'] == 1 and
            receipt['source_seam_witnesses']['full_unwrapped_displacement'] == ['1','0'], 'source roll winding')
    return {
        'result':'PASS', 'input_file_sha256':{k:sha(p) for k,p in paths.items()},
        'identity_spherical_tetrahedra':len(identity_rows),
        'full_cell_identity_argument':'Each nondegenerate tetrahedron has four exact physical source=target vertex equations; affine extension is identity on its entire interior and faces.',
        'identity_chart_index_digest':digest(identity_rows),
        'zero_lift_displacement':[0,4], 'upstairs_full_turns':0, 'downstairs_full_turns':0,
        'F0':'identity', 'Phi0':'identity',
        'filling_product_tetrahedra':len(product), 'filling_refinement_collar_tetrahedra':len(collar),
        'complete_cap_cone_halves':2, 'each_half_pentachora':len(halves[0]),
        'halves_intersection_exactly_filling_torus':True,
        'same_cap_vertex_labels_and_deck':True, 'same_cap_boundary':True,
        'cap_pentachora':5824, 'common_cap_pentachora':len(cap_facets&f2),
        'removed_m2_cap_pentachora':3008, 'inserted_m0_cap_pentachora':3008,
        'unchanged_central_transfer':'u-only; full source winding (1,0)',
        'scope':'Additional finite checks on the frozen open m0 cap and its exact identity attachment. Existing cap/bridge/source receipts remain required; no closed X0 assembled here.'
    }


def compare_candidate(root,candidate):
    oldpath = root/'control/marked_control.json.gz'
    m2,m0 = load(oldpath),load(candidate)
    require((m0['twist'],m0['roll']) == (0,1), 'candidate parameters')
    require(m0['vertices'] == m2['vertices'] and m0['deck'] == m2['deck'], 'candidate labels/deck changed')
    f0,f2 = (set(map(tuple,m['facets'])) for m in (m0,m2))
    unchanged = set()
    for key in ('N','E','central','bridge'):
        require(m0['parts'][key] == m2['parts'][key], 'noncap piece changed: '+key)
        unchanged.update(map(tuple,m0['parts'][key]))
    c0,c2 = (set(map(tuple,m['parts']['cap'])) for m in (m0,m2))
    require(len(unchanged) == 97920, 'unchanged noncap count')
    require(f0 == unchanged|c0 and f2 == unchanged|c2 and not(unchanged & (c0|c2)), 'global part decomposition')
    require(len(f0-f2) == len(f2-f0) == 3008 and f0-f2 == c0-c2 and f2-f0 == c2-c0,
            'global difference is not exactly cap replacement')
    signs0 = dict(zip(map(tuple,m0['facets']),m0['orientation']))
    signs2 = dict(zip(map(tuple,m2['facets']),m2['orientation']))
    require(all(signs0[f] == signs2[f] for f in f0&f2), 'orientation on shared simplex changed')
    for key in ('bulk_boundary','source_port_map','cap_port_map','orientation_anchor','ramification_triangles','twin_triangles'):
        require(m0[key] == m2[key], 'shared marking changed: '+key)
    return {
        'result':'PASS','candidate_file_sha256':sha(candidate),'original_master_file_sha256':sha(oldpath),
        'candidate_facets_sha256':digest(m0['facets']), 'unchanged_noncap_pentachora':len(unchanged),
        'common_cap_pentachora':len(c0&c2), 'common_total_pentachora':len(f0&f2),
        'removed_m2_cap_pentachora':len(f2-f0),'inserted_m0_cap_pentachora':len(f0-f2),
        'identical_vertex_labels_and_deck':True,'identical_source_and_cap_port_maps':True,
        'identical_orientation_on_all_shared_cells':True, 'identical_marked_surface_labels':True,
    }


def main():
    base = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frozen-root',type=Path,default=base/'frozen_release')
    parser.add_argument('--candidate',type=Path)
    parser.add_argument('--output',type=Path,default=base/'review/zero_cap_geometry.json')
    args = parser.parse_args()
    start = time.monotonic()
    result = {'schema':'independent-X0-zero-cap-geometry-review-v1','result':'PASS',
              'checker_sha256':sha(Path(__file__)),'zero_cap':zero_cap_check(args.frozen_root)}
    if args.candidate:
        result['candidate_comparison'] = compare_candidate(args.frozen_root,args.candidate)
    result['wall_seconds'] = time.monotonic()-start
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    main()

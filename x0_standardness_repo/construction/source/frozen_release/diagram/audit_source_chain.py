#!/usr/bin/env python3
"""Read-only audit of signed source, normal orientation and bound interfaces.

Writes only diagram/source_audit_certificate.json.  The complete collar and
ball-replacement verifiers are separate gates; their current, hash-bound
receipts are included rather than described as checks performed here.
"""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import gzip
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def req(x, message):
    if not x:
        raise ValueError(message)


def load(name):
    path = ROOT / name
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


def digest(x):
    return sha256(json.dumps(x, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def filehash(name):
    return sha256((ROOT / name).read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def parity(seq):
    return (-1) ** sum(seq[i] > seq[j] for i in range(len(seq)) for j in range(i + 1, len(seq)))


def det(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def audit():
    diagram = load('diagram/pretzel_237_diagram.json')
    diagram_check = module('independent_diagram', 'diagram/verify_diagram.py').verify(diagram)
    native_check = module('independent_cusp', 'diagram/verify_marked_isomorphism.py')
    exact = native_check.build_certificate(ROOT / 'peripheral')
    req(json.loads(json.dumps(exact)) == load('diagram/marked_isomorphism_certificate.json'),
        'Exact cusp certificate is stale')
    cut = load('peripheral/cut_source.json.gz')
    raw = load('peripheral/peripheral_chart.json.gz')
    chart = load('peripheral/diagram_chart.json.gz')
    ori = load('peripheral/cut_orientation.json')
    imported = load('peripheral/marked_coned_exterior.json.gz')
    quotient = load('peripheral/E_quotient_regina.json.gz')
    iface = load('peripheral/puncture_product_interface.json')
    req(ori['passed'] and ori['cut_sha256'] == digest(cut), 'Cut orientation source mismatch')
    req(chart['raw_chart_sha256'] == digest(raw), 'Corrected chart has another raw source')
    req(chart['native_transport_sha256'] == digest(load('peripheral/marked_exterior_transport.json')),
        'Corrected chart has another native transport')
    req(chart['N_cocycle'] == raw['N_cocycle'], 'Longitude cocycle changed')
    req(chart['E_cocycle'] == [-x for x in raw['E_cocycle']], 'Meridian cocycle not negated')
    for side in ('downstairs', 'upstairs'):
        req(chart[side + '_periods'] == [[u, -v] for u, v in raw[side + '_periods']], 'Period transition changed')
        req([(a, Q(u), Q(v)) for a, u, v in chart[side + '_vertex_coordinates']]
            == [(a, Q(u), -Q(v)) for a, u, v in raw[side + '_vertex_coordinates']], 'Coordinate transition changed')
    req(chart['oriented_core_point_cycle'] == raw['oriented_core_point_cycle'], 'Core orientation changed')
    for key in ('torus_edge_loop', 'raw_exterior_word', 'reduced_exterior_word'):
        req(chart['peripheral_loops']['meridian'][key] == [-x for x in reversed(raw['peripheral_loops']['meridian'][key])],
            'Meridian path inversion missing')
    req(chart['peripheral_loops']['longitude'] == raw['peripheral_loops']['longitude'], 'Longitude path changed')

    # Recheck all exact affine chart determinants, both deck actions, the
    # two-to-one meridional projection, and both peripheral loop periods.
    geometry = module('independent_source_geometry', 'verify_source_geometry.py')
    boundary_N = dict(zip(map(tuple, cut['T']), ori['T_boundary_N_orientations']))
    chart_proof = geometry.chart_check(cut, chart, boundary_N)
    req(all(z['orientation_relative_to_boundary_N'] == -1 for z in chart_proof['charts'].values()),
        'Diagram chart should be negative on outward tube boundary')

    N_orient = dict(zip(map(tuple, cut['N']), ori['N_tetrahedron_orientations']))
    E_orient = dict(zip(map(tuple, cut['E']), ori['E_tetrahedron_orientations']))
    core = chart['oriented_core_point_cycle']
    cocycle = dict(zip(map(tuple, chart['N_skeleton']['edge_representatives']), chart['N_cocycle']))
    req(core[0] == core[-1] and all(cut['deck'][v] == v for v in core), 'Core is not a closed fixed circle')
    degree = sum((1 if a < b else -1) * cocycle[tuple(sorted((a, b)))] for a, b in zip(core, core[1:]))
    req(degree == 1, 'Chosen oriented core is not positive longitude')
    edge = cut['puncture']['fixed_edge_cut_vertices']
    req(edge == core[:2] == iface['poles'], 'Puncture pole order disagrees with positive core')
    raw_cycle = cut['puncture']['equator_cycle']
    positive_cycle = [raw_cycle[0], *reversed(raw_cycle[1:])]
    req(iface['equator_cycle_source'] == raw_cycle and iface['equator_cycle_diagram_positive'] == positive_cycle,
        'Puncture meridian was not reversed exactly once')
    joins = []
    for name, cycle in [('raw', raw_cycle), ('diagram_positive', positive_cycle)]:
        signs = []
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            ordered = (*edge, a, b)
            signs.append(N_orient[tuple(sorted(ordered))] * parity(ordered))
        req(signs == ([-1] * 4 if name == 'raw' else [1] * 4), 'Puncture normal orientation failed')
        joins.append({'cycle_name': name, 'cycle': cycle, 'join_orientation_signs_vs_master': signs})
    tet_multipliers = [imported['tetrahedron_orientations_before_relabeling'][i] * E_orient[tuple(t)]
                       for i, t in enumerate(quotient['tetrahedron_representatives'])]
    req(len(tet_multipliers) == 584 and set(tet_multipliers) == {-1}, 'Imported simplex orientation multiplier')

    # The standard positive affine simplex has determinant +1.  At every
    # truncated vertex the kernel's declared outward corner order has
    # determinant -1/64 against an outward vector.  Therefore kernel
    # right_handed means NEGATIVE standard ordered-simplex orientation.
    points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    calibrations = []
    for v, corner_order in enumerate(native_check.CORNERS):
        small = [tuple(Q(points[v][k]) + Q(points[w][k] - points[v][k], 4) for k in range(3))
                 for w in corner_order]
        centre = tuple(sum(p[k] for p in small) / 3 for k in range(3))
        outward = tuple(Q(points[v][k]) - centre[k] for k in range(3))
        a, b, c = small
        determinant = det(outward, tuple(b[k] - a[k] for k in range(3)), tuple(c[k] - a[k] for k in range(3)))
        req(determinant == Q(-1, 64), 'Kernel orientation calibration failed')
        calibrations.append({'vertex': v, 'outward_corner_order': corner_order,
                             'outward_first_determinant_in_standard_simplex': str(determinant)})
    positive_diagram_basis = diagram['marking']['positive_normal_basis']
    req(det(diagram['marking']['oriented_tangent'], *positive_diagram_basis) == 1, 'Diagram meridian is not right-hand positive')
    wanted_signs = {
        'ambient_orientation_multiplier_vs_Seifert_master': 1,
        'imported_ordered_simplex_multiplier_vs_master': -1,
        'snappy_right_handed_vs_standard_ordered_simplex': -1,
        'boundary_N_vs_du_dv': -1,
        'boundary_E_vs_du_dv': 1,
    }
    req(all(chart['diagram_sign_transport'][k] == v for k, v in wanted_signs.items()), 'Orientation metadata not calibrated')
    req(iface['diagram_chart_sha256'] == digest(chart) and iface['source_cut_sha256'] == digest(cut), 'Puncture interface binding')
    req(iface['ordered_triangles'] == cut['puncture']['boundary_triangles'], 'Puncture triangles changed')
    req(all(iface['deck'][str(positive_cycle[i])] == positive_cycle[(i + 2) % 4] for i in range(4)),
        'Puncture deck is not angular half-turn')

    # Read and bind the independently checked collar/adaptor and cap bridge.
    adaptor = load('adaptor/adaptor.json.gz')
    model = load('adaptor/model.json.gz')
    adaptor_receipt = load('adaptor/independent_verification.json')
    req(adaptor['chart_sha256'] == digest(chart) and adaptor['cut_sha256'] == digest(cut), 'Adaptor source binding')
    req(adaptor['diagram_sign_transport'] == chart['diagram_sign_transport'], 'Adaptor sign metadata stale')
    req(model['diagram_chart_sha256'] == digest(chart) and model['adaptor_sha256'] == digest(adaptor), 'Adaptor model binding')
    req(model['source_cut_sha256'] == digest(cut), 'Adaptor model cut binding')
    req(adaptor_receipt['passed'] and adaptor_receipt['model_sha256'] == digest(model), 'Adaptor verification stale')
    bridge = load('cap/data/source_cap_bridge.json.gz')
    bridge_receipt = load('cap/bridge_verification.json')
    req(bridge['source_interface'] == iface and bridge['source_interface_file_sha256'] == filehash('peripheral/puncture_product_interface.json'), 'Cap bridge source binding stale')
    req(load('cap/data/bridge_generation.json')['canonical_sha256'] == digest(bridge), 'Bridge generation digest stale')
    req(bridge_receipt['result'] == 'PASS' and bridge_receipt['source_angle_cycle'] == positive_cycle
        and bridge_receipt['source_poles'] == edge, 'Bridge orientation receipt mismatch')
    cap_receipt = load('cap/verification.json')
    req(cap_receipt['result'] == 'PASS', 'Cap gate failed')
    req(cap_receipt['cases']['0']['unwrapped_displacement'] == [0, 4]
        and cap_receipt['cases']['2']['unwrapped_displacement'] == [-8, 4], 'Integer cap winding missing')

    paths = [
        'diagram/pretzel_237_diagram.json', 'diagram/verify_diagram.py',
        'diagram/verify_marked_isomorphism.py', 'diagram/marked_isomorphism_certificate.json',
        'diagram/audit_source_chain.py', 'seifert/data/master_orientation.json',
        'coarsened/master_source.json.gz', 'coarsened/source.json.gz',
        'coarsened/coarsening_certificate.json.gz', 'coarsened/coarsening_verified.json',
        'coarsened/source_orientation.json', 'peripheral/cut_source.json.gz',
        'peripheral/cut_orientation.json', 'peripheral/peripheral_chart.json.gz',
        'peripheral/diagram_chart.json.gz', 'peripheral/marked_coned_exterior.json.gz',
        'peripheral/native_import_peripheral_arrays.json.gz',
        'peripheral/marked_coned_exterior.tri', 'peripheral/marked_ideal_exterior.tri',
        'peripheral/matched_marked_exterior.tri', 'peripheral/diagram_exterior.tri',
        'peripheral/marked_exterior_transport.json', 'peripheral/puncture_product_interface.json',
        'transport_chart.py', 'verify_source_geometry.py',
        'adaptor/adaptor.json.gz', 'adaptor/model.json.gz', 'adaptor/initial_stages.json.gz',
        'adaptor/coarsening_certificate.json.gz', 'adaptor/independent_verification.json',
        'verify_torus_adaptor.py', 'cap/data/shared_boundary.json.gz',
        'cap/data/cap_m0.json.gz', 'cap/data/cap_m2.json.gz',
        'cap/data/source_cap_bridge.json.gz', 'cap/verification.json',
        'cap/bridge_verification.json', 'cap/verify_bridge.py',
    ]
    conventions = {
        'kernel/kernel_code/link_complement.c': '1943e1d1695666ecbd7502fe741d84ca12f385911363c0ade6f3f99080f2ce7d',
        'kernel/headers/link_projection.h': '64717c8e23f5973a6d146733c746a9bc5c991fa7a4e751f2db1fed10877f749e',
        'python/snap/peripheral/link.py': '6718cee9889336807cafa76ac0459d7946d787f5d563c275be014f41fc7d2eb0',
        'python/__init__.py': 'ec32a0e700665926426bc3a178680463f40c5a50695776b653caa99c2c112277',
        'cython/core/tail.pyx': 'daa934a228f7e86c970824ae7882592dd330fee43f6ff7165ce1839041465ff2',
    }
    return {
        'schema': 'signed-source-interface-audit-v1', 'passed': True,
        'scope': 'Signed 3D source and its geometric/peripheral/cap interfaces; no faithful X2 construction or standardization is asserted here.',
        'diagram_verification': diagram_check,
        'core': {'oriented_cycle': core, 'N_cocycle_degree': degree, 'puncture_oriented_edge': edge},
        'puncture_normal_orientation': joins,
        'kernel_orientation_calibration': calibrations,
        'orientation_conventions': wanted_signs,
        'imported_compact_tetrahedra_checked': len(tet_multipliers),
        'native_import': exact['native_import_audit'],
        'exact_endpoint_map': exact['all_isomorphisms'][exact['selected_identity_peripheral_isomorphism']],
        'diagram_chart_verification': chart_proof,
        'canonical_digest_links': {
            'cut': digest(cut), 'raw_chart': digest(raw), 'diagram_chart': digest(chart),
            'puncture_interface': digest(iface), 'torus_adaptor': digest(adaptor),
            'torus_model': digest(model), 'cap_bridge': digest(bridge),
        },
        'bound_independent_gates': {
            'torus_adaptor': {'passed': True, 'tetrahedra': adaptor_receipt['tetrahedra'],
                              'model_sha256': digest(model)},
            'cap_bridge': {'passed': True, 'pentachora': bridge_receipt['pentachora'],
                           'ports': bridge_receipt['ports']},
            'cap_winding': {'m0': [0, 4], 'm2': [-8, 4]},
        },
        'primary_sources': [{'url': f'https://raw.githubusercontent.com/3-manifolds/SnapPy/{native_check.PINNED_COMMIT}/{p}',
                             'sha256': h} for p, h in conventions.items()],
        'trust_boundaries': [
            'Spherogram signed-diagram conversion and SnapPy link-complement construction use their pinned native semantics; L follows link direction, M obeys the right-hand rule, and the writhe correction gives the zero longitude.',
            'Source-to-endpoint finite-vertex removal and simplification use trusted native PL moves; this package does not export a complete elementary-move ancestry with all intermediate peripheral arrays.',
            'The exact endpoint identifies oriented peripheral slopes and extends over meridional filling. It does not export a based Wirtinger-generator isomorphism or a pointwise map from the source puncture to the chosen diagram cube; these markings agree up to oriented tubular isotopy.',
            'The torus adaptor and cap bridge are source-bound PL collars. Faithful X2 still requires an actual certified global roll/twist map, its mapping torus and cap attachment, closed native admission, and the requested standardization.',
        ],
        'file_sha256': {p: filehash(p) for p in paths},
    }


if __name__ == '__main__':
    result = audit()
    path = HERE / 'source_audit_certificate.json'
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'passed': result['passed'], 'orientation_conventions': result['orientation_conventions'],
                      'core': result['core'], 'puncture_normal_orientation': result['puncture_normal_orientation'],
                      'canonical_digest_links': result['canonical_digest_links']}, indent=2))

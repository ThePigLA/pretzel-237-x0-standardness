#!/usr/bin/env python3
"""Bind the complete PL-construction derivation to its saved evidence.

This checks hashes, endpoints and receipt consistency. It does not rerun
the individual heavy recognition or reduction programs. Fresh replay is
performed by those programs separately. No native recognition receipt is
created or simulated here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import regina


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bind(root, output):
    root, output = Path(root).resolve(), Path(output)
    base = Path('results/figure2_working_20261102')
    branch = base / 'completion/finite_candidate_phase3_a'
    staircase = base / 'completion/raw_staircase'
    prefix = base / 'native_diagnosis/candidate2136_verified'
    upper = Path('results/boundary/standard_cap_six_boundary_v1')
    upper_reference = Path('results/boundary/upper_cap_selftest_v1')
    upper_guarded = base / 'upper_cap_guarded_v1'
    paths = {
        'source_finite': staircase / 'finite_staircase.rga',
        'source_truncation_receipt': staircase / 'FINITE_TRUNCATION_VERIFICATION.json',
        'source_finite_links_receipt': staircase / 'FINITE_VERTEX_LINKS.json',
        'source_prism_receipt': staircase / 'INDEPENDENT_PRISM_CHECK.json',
        'lower_2136': branch / 'lower_candidate.rga',
        'lower_prefix_receipt': prefix / 'STATUS.json',
        'lower_binding_receipt': prefix / 'CANDIDATE_BINDING.json',
        'lower_move_log': prefix / 'from_original/moves.jsonl',
        'lower_replayed_endpoint': prefix / 'lower_reference.rga',
        'boundary_path': branch / 'short_boundary_path.json',
        'layered_lower': branch / 'layered_lower_fast/layered_piece.rga',
        'layered_lower_events': branch / 'layered_lower_fast/layering_events.jsonl.gz',
        'collar_receipt': branch / 'layered_lower_fast/COLLAR_ATTACHMENT_AUDIT.json',
        'collar_lemma': branch / 'COLLAR_LEMMA.md',
        'collar_checker': Path('code/audit_layering_collar.py'),
        'upper_reference': upper_reference / 'reference_cap.rga',
        'upper_reference_construction': upper_reference / 'REFERENCE_CONSTRUCTION.json',
        'upper_explicit_24': upper_guarded / 'constructed_upper24.rga',
        'upper_guarded_construction': upper_guarded / 'SOURCE_CONSTRUCTION.json',
        'upper_guarded_log': upper_guarded / 'moves.jsonl',
        'upper_guarded_endpoint': upper_guarded / 'guarded_upper9.rga',
        'upper_guarded_result': upper_guarded / 'RESULT.json',
        'upper_guarded_independent_replay': upper_guarded / 'INDEPENDENT_REPLAY.json',
        'upper_474': upper / 'cap_with_matched_boundary.rga',
        'upper_layering_receipt': upper / 'LAYERING_VERIFICATION.json',
        'upper_shelling_receipt': upper / 'SHELLING_BACKCHECK.json',
        'upper_483': branch / 'matched_upper_cap/layered_piece.rga',
        'upper_nine_shellings': branch / 'matched_upper_cap/INVERSE_SHELLING_VERIFICATION.json',
        'closed': branch / 'closed_fast/closed_candidate.rga',
        'closed_gluing': branch / 'closed_fast/boundary_gluing.json',
        'closed_gluing_receipt': branch / 'closed_fast/INDEPENDENT_BOUNDARY_GLUING.json',
        'closed_f_vector': branch / 'closed_fast/RAW_F_VECTOR.json',
        'finite_primal_admission': base / 'layered_primal_admission/FINITE_MANIFOLD_ADMISSION.json',
        'closed_primal_admission': base / 'closed_primal_admission/CLOSED_MANIFOLD_ADMISSION.json',
        'actual_inverse_shellings': branch / 'layered_lower_fast/RAW_INVERSE_SHELLING_VERIFICATION.json',
        'historical_conditional_receipt': branch / 'CONDITIONAL_PL_COMPLETION.json',
    }
    evidence = {key: {'path': str(path), 'sha256': sha(root / path),
                      'bytes': (root / path).stat().st_size}
                for key, path in paths.items()}

    def read(key):
        return json.loads((root / paths[key]).read_text())

    def matches(receipt, field, key):
        require(receipt[field] == evidence[key]['sha256'], field + ': evidence hash differs')

    def tri(key, expected_size):
        obj = regina.open(str(root / paths[key]))
        require(isinstance(obj, regina.Triangulation4) and obj.size() == expected_size,
                key + ': wrong triangulation type or size')
        return obj

    truncation = read('source_truncation_receipt')
    require(truncation['status'] == 'PASS_EXPLICIT_FINITE_TRUNCATION', 'truncation did not pass')
    matches(truncation, 'finite_sha256', 'source_finite')
    require(truncation['finite_manifold_follows_from_explicit_truncation'], 'missing finite premise')
    finite_links = read('source_finite_links_receipt')
    require(finite_links['status'] == 'PASS_FINITE_MANIFOLD_LINK_AUDIT' and
            len(finite_links['vertices']) == 29 and
            all(v['recognized_ball_or_sphere'] for v in finite_links['vertices']),
            'the recorded source link audit did not pass')
    prefix_receipt = read('lower_prefix_receipt')
    require(prefix_receipt['status'] == 'PASS_FROZEN_PREFIX_AND_BOUNDARY_EXPORT',
            'the recorded source reduction did not pass')
    matches(prefix_receipt, 'original_source_sha256', 'source_finite')
    require(prefix_receipt['prefix_moves'] + prefix_receipt['phase_moves'] == 4901,
            'wrong source reduction length')
    lower = tri('lower_2136', 2136)
    require(lower == tri('lower_replayed_endpoint', 2136), 'selected lower differs from replay')
    lower_tight = hashlib.sha256(lower.tightEncoding().encode()).hexdigest()
    require(lower_tight == prefix_receipt['lower_tight_sha256'], 'lower prefix endpoint differs')
    binding = read('lower_binding_receipt')
    require(binding['status'] == 'PASS_CANDIDATE_SOURCE_BINDING_AND_BOUNDARY_ISOMORPHISM',
            'candidate source binding did not pass')
    matches(binding, 'candidate_sha256', 'lower_2136')
    require(binding['candidate_tight_sha256'] == lower_tight, 'candidate tight encoding differs')

    collar = read('collar_receipt')
    require(collar['status'] == 'PASS_ALL_LOCAL_COLLAR_ATTACHMENT_HYPOTHESES' and
            collar['checked_additions'] == 22581 and
            collar['all_embeddings_of_each_interior_support_face_exhausted'] and
            collar['support_open_interior_embedded_and_disjoint_from_boundary_image'],
            'actual collar hypotheses were not all checked')
    for field, key in [('source_sha256', 'lower_2136'), ('path_sha256', 'boundary_path'),
                       ('layered_sha256', 'layered_lower'), ('events_sha256', 'layered_lower_events'),
                       ('checker_sha256', 'collar_checker')]:
        matches(collar, field, key)
    require(tri('layered_lower', 24717).size() == lower.size() + collar['checked_additions'],
            'layered lower size differs')

    # The historical reference was originally found using native simplify.
    # Its identification is now proved by this independently reconstructed
    # and replayed path of 15 guarded edge collapses to the exact same file.
    upper_result = read('upper_guarded_result')
    upper_replay = read('upper_guarded_independent_replay')
    require(upper_replay['status'] == 'PASS_INDEPENDENT_GUARDED_STANDARD_UPPER_CAP' and
            upper_replay['source_pentachora'] == 24 and upper_replay['target_pentachora'] == 9 and
            upper_replay['moves'] == upper_replay['guarded_edge_collapses'] == 15 and
            upper_replay['independent_coordinate_product_reconstruction'] and
            upper_replay['embedded_boundary_sum_tetrahedron_checked'] and
            upper_replay['all_move_guards_checked'] and
            upper_replay['exact_labelled_saved_target_reached'] and
            not upper_replay['native_4d_simplify_used'], 'explicit guarded upper reference proof failed')
    for field, key in [('construction', 'upper_guarded_construction'),
                       ('result', 'upper_guarded_result'), ('moves', 'upper_guarded_log'),
                       ('source', 'upper_explicit_24'), ('endpoint', 'upper_guarded_endpoint'),
                       ('target', 'upper_reference')]:
        require(upper_replay['input_sha256'][field] == evidence[key]['sha256'],
                'guarded upper proof input differs: ' + field)
    require(upper_result['moves_sha256'] == evidence['upper_guarded_log']['sha256'] and
            upper_result['source_sha256'] == evidence['upper_explicit_24']['sha256'] and
            upper_result['endpoint_sha256'] == evidence['upper_reference']['sha256'] ==
            evidence['upper_guarded_endpoint']['sha256'], 'guarded upper result binding differs')
    explicit_upper = tri('upper_explicit_24', 24)
    reference_upper = tri('upper_reference', 9)
    require(hashlib.sha256(explicit_upper.tightEncoding().encode()).hexdigest() ==
            upper_replay['source_tight_sha256'] and
            hashlib.sha256(reference_upper.tightEncoding().encode()).hexdigest() ==
            upper_replay['target_tight_sha256'], 'guarded upper labelled endpoint differs')
    upper_layers = read('upper_layering_receipt')
    require(upper_layers['status'] == 'PASS_EXPLICIT_4D_LAYERING' and upper_layers['moves'] == 465,
            'recorded upper layering did not pass')
    matches(upper_layers, 'source_cap_sha256', 'upper_reference')
    upper_audit = upper_layers['final_cap_audit']
    require(upper_audit['status'] == 'PASS' and upper_audit['pentachora'] == 474 and
            all(v['recognised'] for v in upper_audit['vertex_links']), 'recorded upper finite audit failed')
    upper_back = read('upper_shelling_receipt')
    require(upper_back['status'] == 'PASS_GUARDED_4D_BOUNDARY_SHELLING' and
            upper_back['native_result_exactly_equals_reference'] and
            upper_back['shelling_simplex_indices'] == list(range(473, 8, -1)),
            'recorded upper shelling backcheck failed')
    matches(upper_back, 'source_cap_sha256', 'upper_474')
    matches(upper_back, 'reference_sha256', 'upper_reference')
    upper_nine = read('upper_nine_shellings')
    require(upper_nine['status'] == 'PASS_ALL_INVERSE_4D_SHELLING_GUARDS' and
            upper_nine['moves'] == 9 and upper_nine['source_recovered_exactly'],
            'recorded last nine upper shellings failed')
    matches(upper_nine, 'source_sha256', 'upper_474')
    matches(upper_nine, 'layered_sha256', 'upper_483')
    tri('upper_474', 474); tri('upper_483', 483)

    gluing = read('closed_gluing_receipt')
    require(gluing['status'] == 'PASS_INDEPENDENT_RAW_CLOSED_GLUING' and
            gluing['attaching_map_is_explicit_combinatorial_boundary_isomorphism'] and
            gluing['all_closed_facet_pairings_checked'], 'closed raw gluing audit failed')
    for field, key in [('first_sha256', 'layered_lower'), ('second_sha256', 'upper_483'),
                       ('closed_rga_sha256', 'closed')]:
        matches(gluing, field, key)
    closed = tri('closed', 25200)
    require(hashlib.sha256(closed.tightEncoding().encode()).hexdigest() == gluing['closed_tight_sha256'],
            'closed endpoint tight encoding differs')

    # Bind each individual primal-group proof to its actual saved vertex
    # link, without re-running a recognition algorithm inside this binder.
    admissions = {}
    for label, key, source_key, expected_vertices, required_status in [
            ('finite', 'finite_primal_admission', 'layered_lower', 13,
             'PASS_FINITE_ORIENTABLE_PL_4_MANIFOLD'),
            ('closed', 'closed_primal_admission', 'closed', 19,
             'PASS_CLOSED_ORIENTABLE_PL_4_MANIFOLD')]:
        admission = read(key)
        require(admission['status'] == required_status and
                admission['recognition_method'] ==
                'primal_triangle_presentation_signed_elimination_with_independent_replay' and
                admission['all_face_stabilizers_trivial'] and admission['facet_pairings_reciprocal'] and
                len(admission['vertex_links']) == expected_vertices,
                label + ': raw manifold admission did not pass')
        matches(admission, 'source_sha256', source_key)
        kinds = {'sphere': 0, 'ball': 0}
        for row in admission['vertex_links']:
            require(row['recognised'] and row['independent_elimination_replay_passed'] and
                    row['independent_triangle_relators_regenerated'] and
                    row['final_group_generators'] == 0 and
                    row['recognised_kind'] in kinds,
                    label + ': an individual vertex link was not certified')
            kind = row['recognised_kind']; kinds[kind] += 1
            directory = Path(row['group_proof_directory'])
            if directory.is_absolute():
                directory = directory.resolve().relative_to(root)
            require((root / directory).resolve().is_relative_to(root), 'unsafe group proof directory')
            files = {
                'source': directory / 'explicit_link.rga',
                'presentation': directory / 'elimination/SOURCE.json',
                'result': directory / 'elimination/RESULT.json',
                'log': directory / 'elimination/eliminations.jsonl',
                'replay': directory / 'INDEPENDENT_REPLAY.json',
            }
            hashes = {}
            for part, path in files.items():
                hashes[part] = sha(root / path)
                evidence['%s_vertex_%04d_%s' % (label, row['vertex'], part)] = {
                    'path': str(path), 'sha256': hashes[part], 'bytes': (root / path).stat().st_size}
            require(hashes['source'] == row['link_rga_sha256'] and
                    hashes['result'] == row['primal_result_sha256'] and
                    hashes['log'] == row['primal_elimination_log_sha256'] and
                    hashes['replay'] == row['independent_group_replay_sha256'],
                    label + ': vertex proof artifact checksum differs')
            result = json.loads((root / files['result']).read_text())
            replay = json.loads((root / files['replay']).read_text())
            presentation = json.loads((root / files['presentation']).read_text())
            require(replay['status'] == 'PASS_INDEPENDENT_TRIVIAL_GROUP_PROOF' and
                    replay['recognised_link_kind_by_poincare_and_3d_schoenflies'] == kind and
                    replay['source_sha256'] == result['source_sha256'] ==
                    presentation['source_sha256'] == hashes['source'] and
                    replay['proof_result_sha256'] == hashes['result'] and
                    replay['elimination_log_sha256'] == result['elimination_log_sha256'] == hashes['log'] and
                    result['all_original_generators_proved_identity'] and
                    replay['remaining_generator_classes'] == result['remaining_generator_classes'] == 0,
                    label + ': vertex proof semantics or source binding differs')
        if label == 'closed':
            require(kinds == {'sphere': 19, 'ball': 0} and admission['euler_characteristic'] == 2,
                    'closed links or Euler characteristic differ')
        admissions[label] = {'receipt': key, 'vertex_link_kinds': kinds,
                             'euler_characteristic': admission['euler_characteristic'],
                             'all_individual_group_proof_artifacts_bound': True}

    shellings = read('actual_inverse_shellings')
    require(shellings['status'] == 'PASS_ALL_RAW_INVERSE_SHELLING_GUARDS' and
            shellings['moves'] == 22581 and shellings['source_recovered_exactly'] and
            shellings['initial_finite_PL_manifold_admitted'] and
            shellings['all_face_validity_preserved_inductively'],
            'actual inverse shelling replay did not pass')
    for field, key in [('source_sha256', 'lower_2136'), ('layered_sha256', 'layered_lower'),
                       ('admission_sha256', 'finite_primal_admission')]:
        matches(shellings, field, key)
    report = {
        'status': 'PASS_COMPLETE_PL_CONSTRUCTION_CHAIN',
        'supersedes': str(paths['historical_conditional_receipt']),
        'working_diagram_assumption': 'Figure 2 red passes over dashed green',
        'derivation': [
            'The finite5952 source has an explicit truncation construction and a recorded29-link PL-manifold audit.',
            'The independently guarded4901-move source path ends at the exact selected finite2136 lower piece.',
            'All22581 actual4D attachments satisfy the collar lemma, so layered24717 is PL homeomorphic to finite2136.',
            'The explicit standard24-pentachoron upper reference reaches the exact saved9-pentachoron reference by15 guarded edge collapses, independently replayed from a fresh coordinate product construction.',
            'The483 upper piece comes from that independently certified reference by474 checked layers.',
            'Independent raw manifold admissions certify all13 finite-lower vertex links and all19 closed vertex links through regenerated and independently replayed exact primal group proofs.',
            'All22581 raw inverse shelling guards pass and recover the exact2136 source, independently of the forward collar-lemma certificate.',
            'The recorded seven-tetrahedron boundary bijection is an exact combinatorial isomorphism of the two boundaries.',
            'Gluing the two finite PL pieces along that isomorphism gives the saved closed25200 PL manifold.',
        ],
        'source_finite_pl_premise': {
            'explicitly_required': True,
            'supplied_by': ['source_truncation_receipt', 'source_finite_links_receipt',
                            'source_prism_receipt', 'lower_prefix_receipt', 'lower_binding_receipt'],
            'regina_guarded_move_semantics_and_3d_link_recognition_trusted': True,
        },
        'upper_reference_certificate': {
            'construction': 'boundary connected sum of two explicit Delta3 times a three-edge circle products',
            'raw_pentachora': 24, 'native_simplified_pentachora': 9,
            'guarded_edge_collapses_24_to_9': 15,
            'independently_reconstructed_and_replayed': True,
            'unlogged_native_regina_simplification_24_to_9_is_a_proof_premise': False,
            'historical_native_simplification_receipt_preserved_but_superseded':
                str(paths['upper_reference_construction']),
            'guarded_layer_additions_after_reference': 474,
            'guarded_regina_edge_collapse_semantics_trusted': True,
        },
        'conclusions_from_the_bound_certificates': {
            'layered_lower_is_finite_pl_and_pl_homeomorphic_to_source_2136': True,
            'saved_gluing_is_a_closed_pl_four_manifold': True,
        },
        'independent_raw_manifold_admissions': admissions,
        'actual_raw_inverse_shelling_guards_passed': 22581,
        'new_vertex_link_recognition_claimed': True,
        'vertex_link_recognition_basis': 'Verified finite3-manifold hypotheses and independently replayed trivial fundamental group proofs; Poincare, dimension3 PL compatibility, and PL3 Schoenflies for ball links.',
        'native_recognition_caches_populated_or_modified': False,
        'giant_vertex_link_recognition_receipt_synthesized': False,
        'forward_collar_lemma_and_actual_inverse_shelling_certificates_both_pass': True,
        'individual_historical_recognition_and_reduction_programs_rerun_by_this_binder': False,
        'standardness_asserted': False,
        'evidence': evidence,
        'binder_sha256': sha(__file__),
    }
    temporary = output.with_suffix(output.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(output)
    print(json.dumps({'status': report['status'], 'evidence_files_bound': len(evidence),
                      'output': str(output), 'output_sha256': sha(output)}), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=Path(__file__).resolve().parents[1], type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    bind(args.root, args.out)

#!/usr/bin/env python3
"""Try a direct trivial-group proof on one explicitly assembled 3D link.

No geometric simplification is used. This certifies only the selected link,
not the source 4-manifold. Run under an external time limit for a bounded probe.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import regina
from audit_raw_closed4 import raw_pairings, explicit_vertex_links, require


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def probe(source, out, vertex=0, finite=False, extract_only=False):
    require(importlib.metadata.version('regina') == '7.4.1',
            'requires Regina distribution 7.4.1')
    source, out = Path(source), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    source_sha = sha(source)
    tri = regina.Triangulation4(regina.open(str(source)))
    require(sha(source) == source_sha, 'source changed while loading')
    pairs, _ = raw_pairings(tri, closed=not finite)
    link = None
    for index, corners, candidate in explicit_vertex_links(tri, pairs):
        if index == vertex:
            link = candidate
            break
    require(link is not None, 'vertex index does not exist')
    require(link.isValid() and not link.isIdeal() and link.isConnected()
            and link.isOrientable(),
            'requires a valid finite connected orientable 3-manifold link')
    boundary = link.hasBoundaryFacets()
    if boundary:
        require(finite, 'closed source has a boundary vertex link')
        require(link.countBoundaryComponents() == 1,
                'ball criterion requires exactly one boundary component')
        bc = link.boundaryComponent(0)
        require(bc.isOrientable() and bc.eulerChar() == 2,
                'ball criterion requires an orientable S2 boundary')
        boundary_info = {'components': 1, 'orientable': True,
                         'euler_characteristic': 2,
                         'triangles': bc.countTriangles()}
    else:
        require(link.isClosed(), 'sphere criterion requires a closed link')
        boundary_info = {'components': 0}
    link_path = out / 'explicit_link.rga'
    require(regina.PacketOfTriangulation3(link).save(str(link_path)),
            'could not save explicit link')
    result = {'source_sha256': source_sha, 'vertex': vertex,
              'first_local_corner': list(corners[0]),
              'tetrahedra': link.size(), 'link_rga_sha256': sha(link_path),
              'valid_finite_connected_orientable_3_manifold': True,
              'closed': not boundary, 'boundary': boundary_info,
              'geometric_simplification_used': False,
              'four_dimensional_skeleton_used': False,
              'cached_flags_modified': False,
              'source_four_manifold_admission_asserted': False,
              'standardness_of_source_asserted': False}
    def progress(stage):
        result.update(stage=stage, seconds=time.monotonic() - started)
        atomic_json(out / 'PROGRESS.json', result)
        print(json.dumps(result), flush=True)
    progress('LINK_VALIDATED_EXPLICIT_EXTRACTION' if extract_only
             else 'LINK_VALIDATED_BEGIN_FUNDAMENTAL_GROUP')
    if extract_only:
        result['status'] = 'PASS_EXPLICIT_FINITE_3_MANIFOLD_LINK_EXTRACTION'
        result['stage'] = 'EXTRACTION_ONLY_NO_GROUP_COMPUTATION'
        result['group_triviality_asserted'] = False
        result['checker_sha256'] = sha(__file__)
        atomic_json(out / 'LINK_EXTRACTION.json', result)
        print(json.dumps(result), flush=True)
        return result
    group = link.group()
    result['initial_group_generators'] = group.countGenerators()
    result['initial_group_relations'] = group.countRelations()
    progress('FUNDAMENTAL_GROUP_COMPUTED')
    if group.countGenerators():
        result['further_presentation_simplification_used'] = True
        group.simplify()
    else:
        result['further_presentation_simplification_used'] = False
    result['final_group_generators'] = group.countGenerators()
    result['final_group_relations'] = group.countRelations()
    trivial = group.countGenerators() == 0
    result['status'] = ('PASS_LINK_IS_3_BALL_BY_TRIVIAL_GROUP' if boundary
                        else 'PASS_LINK_IS_3_SPHERE_BY_TRIVIAL_GROUP') if trivial \
                       else 'INCONCLUSIVE_GROUP_PRESENTATION'
    result['recognised'] = trivial
    result['theorem'] = (
        'A finite compact connected orientable 3-manifold with one S2 boundary '
        'and trivial fundamental group is a 3-ball, by capping the boundary, '
        'Poincare, and the PL 3-dimensional Schoenflies theorem.' if boundary else
        'A finite closed connected orientable 3-manifold with trivial '
        'fundamental group is a 3-sphere, by Poincare and dimension-three '
        'topological/PL compatibility.')
    result['trust_boundary'] = ('Regina 7.4.1 three-dimensional validity, '
                                'fundamental-group construction and exact '
                                'presentation simplification; explicit raw '
                                'link assembly and the stated 3-manifold theorem.')
    result['seconds'] = time.monotonic() - started
    result['checker_sha256'] = sha(__file__)
    atomic_json(out / 'LINK_GROUP_ADMISSION.json', result)
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--vertex', type=int, default=0)
    parser.add_argument('--finite', action='store_true')
    parser.add_argument('--extract-only', action='store_true')
    args = parser.parse_args()
    probe(args.source, args.out, args.vertex, args.finite, args.extract_only)

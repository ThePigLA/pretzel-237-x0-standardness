#!/usr/bin/env python3
"""Admit a raw finite/closed PL 4-manifold using algebraic 3D link tests.

No 4D skeleton, geometric link simplification, isSphere(), or isBall() is
used.  Valid finite closed 3-manifold links with zero-generator fundamental
groups are spheres by Poincare and dimension-three PL compatibility.
Boundary links must have exactly one explicitly verified S2 boundary; the
same group condition then gives a ball by capping and PL Schoenflies.

Failure to simplify a presentation to zero generators is inconclusive.
The primal method reconstructs the cellular presentation and independently
replays every signed elimination.  The optional dual method instead trusts
Regina 7.4.1 group construction and presentation transformations.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import regina

from audit_raw_closed4 import raw_pairings,face_stabilizers,explicit_vertex_links,require


METHOD='dual_2_skeleton_trivial_group_poincare_pl_schoenflies'
PRIMAL_METHOD='primal_triangle_presentation_signed_elimination_with_independent_replay'
CLOSED_THEOREM=(
    'A finite closed connected orientable 3-manifold with trivial fundamental '
    'group is homeomorphic to S3 by the three-dimensional Poincare theorem. '
    'Uniqueness of PL structures in dimension three identifies its given PL '
    'triangulation with the standard PL 3-sphere.')
BALL_THEOREM=(
    'A finite compact connected orientable PL 3-manifold with exactly one '
    'boundary component, verified to be PL S2, and trivial fundamental group '
    'is a PL 3-ball: cap its boundary by a PL 3-ball; van Kampen preserves '
    'trivial fundamental group, Poincare and dimension-three PL compatibility '
    'give PL S3, and the locally flat PL 2-sphere bounds PL 3-balls by '
    'three-dimensional PL Schoenflies.')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def link_hypotheses(link,allow_boundary=False):
    """Verify exactly the manifold/boundary hypotheses for the group criterion."""
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina distribution 7.4.1')
    require(link.size()>0,'empty 3D link')
    require(link.isValid(),'invalid 3D link')
    require(not link.isIdeal(),'3D link has ideal vertices; a finite manifold is required')
    require(link.isConnected(),'disconnected 3D link')
    require(link.isOrientable(),'nonorientable 3D link')
    has_boundary=link.hasBoundaryFacets()
    if has_boundary:
        require(allow_boundary,'closed 4D input has a link with real boundary')
        require(link.countBoundaryComponents()==1,'ball criterion requires exactly one boundary component')
        bc=link.boundaryComponent(0)
        require(bc.isReal(),'boundary component is not real')
        surface=bc.build()
        require(surface.size()>0 and surface.isValid() and surface.isConnected()
                and surface.isClosed() and surface.isOrientable()
                and surface.eulerChar()==2,
                'ball criterion requires a valid connected closed orientable Euler-characteristic-2 boundary surface')
        boundary={'components':1,'real':True,'triangles':surface.size(),
                  'surface_valid':True,'surface_connected':True,'surface_closed':True,
                  'surface_orientable':True,'surface_euler_characteristic':2,
                  'surface_is_PL_S2_by_surface_classification':True}
        kind='ball';theorem=BALL_THEOREM
    else:
        require(link.isClosed() and link.countBoundaryComponents()==0,
                'sphere criterion requires a closed finite 3-manifold')
        boundary={'components':0};kind='sphere';theorem=CLOSED_THEOREM
    facts={'recognition_method':METHOD,'tetrahedra':link.size(),
           'valid_finite_connected_orientable_3_manifold':True,
           'closed':not has_boundary,'boundary':boundary,
           'link_tight_sha256':hashlib.sha256(link.tightEncoding().encode()).hexdigest(),
           'recognised_kind':kind,'theorem':theorem,
           'geometric_simplification_used':False,
           'native_3d_sphere_ball_recognition_used':False,
           'cached_flags_modified':False}
    return facts


def recognise_link_by_group(link,allow_boundary=False,progress=None):
    """A sufficient sphere/ball criterion using Regina's dual presentation."""
    started=time.monotonic();facts=link_hypotheses(link,allow_boundary)
    kind=facts['recognised_kind']
    if progress:progress('LINK_HYPOTHESES_VERIFIED_BEGIN_GROUP',facts)
    group_start=time.monotonic()
    # Copy the returned presentation; do not mutate a cached group object.
    group=regina.GroupPresentation(link.group())
    require(group.isValid(),'Regina returned an invalid group presentation')
    facts.update(group_construction='Regina Triangulation3.group(), computed from the dual 2-skeleton',
                 returned_group_generators=group.countGenerators(),
                 returned_group_relations=group.countRelations(),
                 returned_group_relator_length=group.relatorLength(),
                 group_construction_seconds=time.monotonic()-group_start)
    if progress:progress('GROUP_CONSTRUCTED',facts)
    simplify_start=time.monotonic();further=group.countGenerators()!=0
    change=group.simplify() if further else None
    require(group.isValid(),'simplified group presentation is invalid')
    recognised=group.countGenerators()==0
    facts.update(further_group_simplification_used=further,
                 group_simplification_algorithm='Regina GroupPresentation.simplify(): small-cancellation and Nielsen moves; returned changes are declared isomorphisms',
                 simplification_returned_declared_isomorphism=change is not None,
                 group_simplification_seconds=time.monotonic()-simplify_start,
                 final_group_generators=group.countGenerators(),
                 final_group_relations=group.countRelations(),
                 final_group_relator_length=group.relatorLength(),
                 recognised=recognised,is_sphere=recognised and kind=='sphere',
                 is_ball=recognised and kind=='ball',
                 status=('PASS_LINK_IS_3_BALL_BY_TRIVIAL_GROUP' if kind=='ball'
                         else 'PASS_LINK_IS_3_SPHERE_BY_TRIVIAL_GROUP') if recognised
                        else 'INCONCLUSIVE_GROUP_PRESENTATION',
                 seconds=time.monotonic()-started,
                 trust_boundary='Regina 7.4.1 3D validity and finite-boundary classification, dual-2-skeleton fundamental-group construction, and exact presentation operations; surface classification, Poincare, Moise PL compatibility, and PL 3D Schoenflies.')
    return facts


def recognise_link_by_primal(link,proof_dir,allow_boundary=False,seconds=120,progress=None):
    """The same theorem, supplied with an independently replayed primal proof."""
    from primal_triangle_group import solve
    from verify_primal_triangle_group import verify
    started=time.monotonic();facts=link_hypotheses(link,allow_boundary)
    facts['recognition_method']=PRIMAL_METHOD
    kind=facts['recognised_kind'];proof_dir=Path(proof_dir)
    proof_dir.mkdir(parents=True,exist_ok=False)
    source=proof_dir/'explicit_link.rga'
    require(regina.PacketOfTriangulation3(link).save(str(source)),'could not save explicit primal-proof link')
    link_sha=sha(source)
    require(hashlib.sha256(regina.open(str(source)).tightEncoding().encode()).hexdigest()==facts['link_tight_sha256'],
            'saved primal input differs from the explicitly assembled link')
    if progress:progress('LINK_HYPOTHESES_VERIFIED_BEGIN_PRIMAL_ELIMINATION',facts)
    produced=solve(source,proof_dir/'elimination',seconds=seconds)
    if progress:progress('PRIMAL_PROOF_PRODUCED_BEGIN_INDEPENDENT_REPLAY',facts)
    checked=verify(source,proof_dir/'elimination',proof_dir/'INDEPENDENT_REPLAY.json')
    require(sha(source)==link_sha,'explicit link changed during primal proof')
    recognised=(produced['status']=='PASS_TRIVIAL_PRIMAL_FUNDAMENTAL_GROUP'
                and checked['status']=='PASS_INDEPENDENT_TRIVIAL_GROUP_PROOF'
                and checked['remaining_generator_classes']==0)
    if recognised:
        require(checked['recognised_link_kind_by_poincare_and_3d_schoenflies']==kind,
                'independent link-kind conclusion differs')
    facts.update(recognised=recognised,is_sphere=recognised and kind=='sphere',
                 is_ball=recognised and kind=='ball',
                 status=('PASS_LINK_IS_3_BALL_BY_TRIVIAL_GROUP' if kind=='ball'
                         else 'PASS_LINK_IS_3_SPHERE_BY_TRIVIAL_GROUP') if recognised
                        else 'INCONCLUSIVE_PRIMAL_GROUP_ELIMINATION',
                 link_rga_sha256=link_sha,group_proof_directory=str(proof_dir),
                 group_construction='Primal spanning-tree edge generators and oriented triangle relators',
                 group_simplification_algorithm='Only justified one-generator kills and signed merges of two distinct roots; torsion relations are retained',
                 primal_initial_generators=produced['generators'],
                 primal_eliminations=checked['eliminations'],
                 final_group_generators=checked['remaining_generator_classes'],
                 independent_triangle_relators_regenerated=checked['triangle_relators_regenerated_via_tetrahedron_embeddings'],
                 independent_elimination_replay_passed=checked['status'].startswith('PASS_'),
                 primal_result_sha256=sha(proof_dir/'elimination/RESULT.json'),
                 primal_elimination_log_sha256=sha(proof_dir/'elimination/eliminations.jsonl'),
                 independent_group_replay_sha256=sha(proof_dir/'INDEPENDENT_REPLAY.json'),
                 regina_group_computation_used=False,
                 seconds=time.monotonic()-started,
                 trust_boundary='Regina 7.4.1 3D validity and face/edge mappings; independently regenerated cellular group presentation and exact signed-elimination replay; surface classification, Poincare, Moise PL compatibility and PL 3D Schoenflies.')
    return facts


def audit(path,out,save_links=False,closed=True,method='dual',group_seconds=120):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina distribution 7.4.1')
    require(method in ('dual','primal'),'unknown group method')
    path,out=Path(path),Path(out);out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();source_sha=sha(path)
    packet=regina.open(str(path));require(isinstance(packet,regina.Triangulation4),'not a 4D input')
    require(sha(path)==source_sha,'source changed while loading')
    tri=regina.Triangulation4(packet)
    pairs,orientation=raw_pairings(tri,closed=closed)
    counts=face_stabilizers(tri.size(),pairs)
    print(json.dumps({'stage':'FACET_AND_FACE_CHECKS_PASSED','pentachora':tri.size(),
                      'seconds':time.monotonic()-started}),flush=True)
    vertices=[]
    for index,embeddings,link in explicit_vertex_links(tri,pairs):
        def progress(stage,facts):
            value={'stage':stage,'vertex':index,'source_sha256':source_sha,
                   'completed_vertex_links':len(vertices),'seconds':time.monotonic()-started,
                   'link_facts':facts}
            atomic_json(out/'PROGRESS.json',value)
            print(json.dumps({'stage':stage,'vertex':index,'tetrahedra':link.size(),
                              'seconds':time.monotonic()-started}),flush=True)
        try:
            if method=='primal':
                row=recognise_link_by_primal(link,out/f'link_{index:04d}_primal',
                    allow_boundary=not closed,seconds=group_seconds,progress=progress)
            else:row=recognise_link_by_group(link,allow_boundary=not closed,progress=progress)
        except ValueError as error:
            atomic_json(out/'REJECTION.json',{'status':'FAIL_LINK_HYPOTHESES',
                'source_sha256':source_sha,'vertex':index,'reason':str(error),
                'four_manifold_admission_asserted':False,'standardness_asserted':False})
            raise
        row.update(vertex=index,original_tetrahedra=link.size(),first_local_corner=list(embeddings[0]),
                   native_3d_simplification_used=False)
        if save_links and method=='dual':
            link_path=out/f'link_{index:04d}.rga'
            require(regina.PacketOfTriangulation3(link).save(str(link_path)),'could not save explicit link')
            row['link_rga_sha256']=sha(link_path)
        vertices.append(row);atomic_json(out/'VERTEX_LINKS.json',vertices)
        print(json.dumps({'vertex':index,'recognised':row['recognised'],
                          'kind':row['recognised_kind'],'generators':row['final_group_generators'],
                          'seconds':row['seconds']}),flush=True)
        if not row['recognised']:
            atomic_json(out/'REJECTION.json',{'status':'INCONCLUSIVE_LINK_GROUP_PRESENTATION',
                'source_sha256':source_sha,'vertex':index,
                'reason':'The presentation did not simplify to zero generators. This does not prove the link is not a sphere or ball.',
                'four_manifold_admission_asserted':False,'standardness_asserted':False})
            raise ValueError(f'algebraic criterion did not certify vertex link {index}')
    require(sha(path)==source_sha,'source changed during verification')
    vector=[len(vertices),counts[1],counts[2],counts[3],tri.size()]
    result={'status':'PASS_CLOSED_ORIENTABLE_PL_4_MANIFOLD' if closed else 'PASS_FINITE_ORIENTABLE_PL_4_MANIFOLD',
            'closed_required':closed,'recognition_method':PRIMAL_METHOD if method=='primal' else METHOD,
            'real_boundary_facets':sum(s.adjacentSimplex(f) is None for s in tri.simplices() for f in range(5)),
            'source_file':str(path.resolve()),'source_sha256':source_sha,
            'source_tight_sha256':hashlib.sha256(tri.tightEncoding().encode()).hexdigest(),
            'f_vector':vector,'euler_characteristic':sum((-1)**i*n for i,n in enumerate(vector)),
            'all_face_stabilizers_trivial':True,'facet_pairings_reciprocal':True,
            'simplex_orientations':orientation,'vertex_links':vertices,
            'four_dimensional_skeleton_used':False,'native_4d_recognition_used':False,
            'native_3d_recognition_used':False,'native_3d_group_operations_used':method=='dual',
            'independent_group_elimination_replay_required':method=='primal',
            'geometric_simplification_used':False,'cached_flags_modified':False,
            'simply_connectedness_asserted':False,'standardness_asserted':False,
            'group_algorithm':('Primal spanning-tree triangle presentation, signed length-one/length-two elimination, and independent reconstruction/replay; every original generator must be proved identity.' if method=='primal' else 'Regina 7.4.1 dual-2-skeleton group construction, followed when needed by deterministic small-cancellation/Nielsen presentation simplification; zero generators is the only accepting group criterion.'),
            'theorem_basis':{'closed_links':CLOSED_THEOREM,'boundary_links':BALL_THEOREM},
            'trust_boundary':('Raw 4D facet and face-stabilizer checks and explicit vertex-link assembly; Regina 7.4.1 3D validity/boundary and edge-mapping APIs; independently replayed exact primal group elimination; the stated 3-manifold theorems.' if method=='primal' else 'Raw 4D facet and face-stabilizer checks and explicit vertex-link assembly; Regina 7.4.1 3D validity/boundary classification and exact group algorithms; the stated 3-manifold theorems. No independent Tietze trace is supplied.'),
            'regina_distribution':'7.4.1','seconds':time.monotonic()-started,'checker_sha256':sha(__file__),
            'raw_assembly_checker_sha256':sha(Path(__file__).with_name('audit_raw_closed4.py'))}
    output=out/('CLOSED_MANIFOLD_ADMISSION.json' if closed else 'FINITE_MANIFOLD_ADMISSION.json')
    atomic_json(output,result)
    print(json.dumps({k:result[k] for k in ('status','recognition_method','f_vector','euler_characteristic','seconds')}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--finite',action='store_true');p.add_argument('--save-links',action='store_true')
    p.add_argument('--method',choices=('dual','primal'),default='dual')
    p.add_argument('--group-seconds',type=float,default=120,
                   help='per-link primal search bound, including extraction; independent replay is separate')
    a=p.parse_args();audit(a.source,a.out,a.save_links,closed=not a.finite,
                          method=a.method,group_seconds=a.group_seconds)

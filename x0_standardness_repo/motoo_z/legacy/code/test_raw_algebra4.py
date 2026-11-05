#!/usr/bin/env python3
"""Small positive/negative controls for the algebraic 3D link criterion."""
import argparse
import hashlib
import json
from pathlib import Path
import regina
from audit_raw_algebra4 import recognise_link_by_group,audit,atomic_json,require


def controls(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    tests=[('sphere',regina.Example3.sphere(),False,True),
           ('ball',regina.Example3.ball(),True,True),
           ('S2_times_S1',regina.Example3.s2xs1(),False,False),
           ('Poincare_homology_sphere',regina.Example3.poincare(),False,False),
           ('solid_torus',regina.Example3.handlebody(1),True,False),
           ('ideal_figure_eight',regina.Example3.figureEight(),False,False)]
    # S3 with interiors of two disjoint tetrahedra removed: simply connected,
    # but two S2 boundary components. This must not be called a 3-ball.
    twice_punctured=regina.Example3.simplicialSphere();twice_punctured.subdivide()
    first={twice_punctured.simplex(0).vertex(v).index() for v in range(4)}
    other=next(s.index() for s in twice_punctured.simplices()
               if first.isdisjoint(s.vertex(v).index() for v in range(4)))
    twice_punctured.removeSimplexAt(other);twice_punctured.removeSimplexAt(0)
    require(twice_punctured.isValid() and not twice_punctured.isIdeal()
            and twice_punctured.countBoundaryComponents()==2,'two-boundary control construction failed')
    require(twice_punctured.group().countGenerators()==0,'two-boundary control did not have certified trivial group')
    tests.append(('two_S2_boundaries_trivial_group',twice_punctured,True,False))
    disconnected=regina.Triangulation3();disconnected.newSimplices(2)
    tests.append(('two_disjoint_balls',disconnected,True,False))
    results=[]
    for name,tri,allow_boundary,wanted in tests:
        try:
            facts=recognise_link_by_group(tri,allow_boundary)
            accepted=facts['recognised'];error=None
        except ValueError as exc:
            facts=None;accepted=False;error=str(exc)
        require(accepted==wanted,'unexpected link-control result: '+name)
        row={'name':name,'tetrahedra':tri.size(),'expected_accept':wanted,
             'accepted':accepted,'rejection_reason':error,'facts':facts}
        if name=='Poincare_homology_sphere':row['first_homology']=str(tri.homology())
        if name=='two_S2_boundaries_trivial_group':row['certified_group_generators']=0
        results.append(row)
    # Full 4D negative: cone over the Poincare sphere has a non-spherical
    # apex link despite trivial first homology. It must not get finite admission.
    cone=regina.Example4.singleCone(regina.Example3.poincare())
    source=out/'cone_over_poincare.rga'
    require(regina.PacketOfTriangulation4(cone).save(str(source)),'could not save negative cone')
    rejected=False
    try:audit(source,out/'cone_negative_audit',closed=False)
    except ValueError:rejected=True
    require(rejected and not (out/'cone_negative_audit/FINITE_MANIFOLD_ADMISSION.json').exists(),
            'negative cone received finite-manifold admission')
    result={'status':'PASS_ALGEBRAIC_LINK_POSITIVE_AND_NEGATIVE_CONTROLS',
            'link_controls':results,'cone_over_poincare_rejected':True,
            'negative_results_do_not_assert_nontriviality_for_arbitrary_presentations':True,
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    atomic_json(out/'CONTROL_RESULTS.json',result)
    print(json.dumps({'status':result['status'],'link_controls':len(results),
                      'cone_over_poincare_rejected':True}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True,type=Path)
    a=p.parse_args();controls(a.out)

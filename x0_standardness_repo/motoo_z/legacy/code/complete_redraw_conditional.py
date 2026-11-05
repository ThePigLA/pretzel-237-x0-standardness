#!/usr/bin/env python3
"""Prepare/complete the conditional dashed-back-arc reading of redraw Figure2.

The sole source-convention assumption is stated in CONDITIONAL_SCOPE.json.
No source isotopy or closed-manifold standardness is asserted. Native Regina
ideal truncation and initial finite-body simplify are explicit trust boundaries.
Subsequent boundary Pachner moves, layering and gluing use existing replay tools.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import threading
import time

import regina

from complete_upper_handles import boundary_data, glue_pieces, layer_path, require, save3, save4, save_json
from boundary_pachner_path import replay as replay_boundary, search


PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / 'results/redraw_20261102/katie_figure2_dashed_v1/raw_katie.rga'
CAP = PROJECT / 'results/boundary/standard_cap_six_boundary_v1/cap_with_matched_boundary.rga'
OUT = PROJECT / 'results/redraw_20261102/figure2_completion_conditional_v1'
GOAL = 'gLALQccbefefdpcghfc'
NATIVE_SECONDS = 60


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load4(path):
    p = regina.open(str(path))
    require(isinstance(p, regina.Triangulation4), 'expected a four-dimensional triangulation')
    return regina.Triangulation4(p)


def bounded_native_simplify(tri):
    """Cancel via Regina's tracker so completed native improvements survive."""
    tracker=regina.ProgressTrackerObjective(tri.size())
    result={}
    def work():
        try:
            result['changed']=bool(tri.simplify(tracker))
        except BaseException as exc:
            result['error']=repr(exc)
    worker=threading.Thread(target=work)
    now=time.monotonic();last=now
    worker.start()
    while worker.is_alive():
        worker.join(timeout=1)
        elapsed=time.monotonic()-now
        if time.monotonic()-last>=10:
            print(json.dumps({'phase':'native_simplify_progress','elapsed':elapsed,
                              'best_pentachora':tracker.objective(),
                              'stage':tracker.description(),'cancel_requested':tracker.isCancelled()}),flush=True)
            last=time.monotonic()
        if elapsed>=NATIVE_SECONDS and not tracker.isCancelled():
            tracker.cancel()
    require('error' not in result,result.get('error','native simplifier failed'))
    result.update({'elapsed_seconds':time.monotonic()-now,
                   'cancel_requested':tracker.isCancelled(),
                   'tracked_best_pentachora':tracker.objective(),
                   'requested_seconds':NATIVE_SECONDS})
    return result


def prepare():
    require(not OUT.exists(), 'output already exists; preserve preparation rather than overwrite')
    OUT.mkdir(parents=True)
    shutil.copy2(RAW, OUT / 'raw_katie_bound_source.rga')
    scope = {
        'status': 'CONDITIONAL_SOURCE_INTERPRETATION_NOT_CONFIRMED',
        'figure': 2,
        'source_pdf': 'inputs/redraw_20261102/figure2.pdf',
        'source_pdf_sha256': sha(PROJECT/'inputs/redraw_20261102/figure2.pdf'),
        'assumption': 'At assembled-F2-v1 crossing68, PDF(825.6894441053033,442.6894441053033), the red strand passes OVER the dashed green back arc. Literal paint order instead puts green over red.',
        'source_confirmation_received': False,
        'isotopy_to_previous_Figure2_asserted': False,
        'standardness_asserted': False,
        'source_katie_raw': str(RAW.relative_to(PROJECT)),
        'source_katie_raw_sha256': sha(RAW),
        'bound_raw_copy_sha256': sha(OUT/'raw_katie_bound_source.rga'),
        'source_katie_input_sha256': sha(RAW.parent/'input.json'),
        'source_katie_receipt_sha256': sha(RAW.parent/'RUN_RECEIPT.json'),
        'reference_upper_cap': str(CAP.relative_to(PROJECT)),
        'reference_upper_cap_sha256': sha(CAP),
        'regina_python_distribution': importlib.metadata.version('regina'),
        'native_trust_boundaries': ['Katie conversion', 'Regina ideal truncation', 'Regina finite-body simplify'],
    }
    save_json(OUT/'CONDITIONAL_SCOPE.json', scope)
    lower = load4(RAW)
    require(lower.isValid() and lower.isConnected() and lower.isOrientable() and lower.isIdeal(), 'raw source is not one valid orientable ideal body')
    require(lower.countBoundaryComponents()==1, 'raw source must have one boundary')
    before={'pentachora':lower.size(),'vertices':lower.countVertices(),'ideal_boundary_components':lower.countBoundaryComponents()}
    print(json.dumps({'phase':'raw_validated',**before}),flush=True)
    now=time.monotonic()
    require(lower.truncateIdeal(),'ideal truncation made no change')
    save4(lower,OUT/'truncated_raw.rga')
    truncated_size=lower.size()
    print(json.dumps({'phase':'truncated_saved','pentachora':truncated_size,'elapsed':time.monotonic()-now}),flush=True)
    save_json(OUT/'TRUNCATION_PREPARATION.json',{
        'raw':before,'source_raw_sha256':sha(RAW),
        'truncated_sha256':sha(OUT/'truncated_raw.rga'),
        'truncated_pentachora':truncated_size,
        'operation':'Regina Triangulation4.truncateIdeal()',
        'large_truncated_vertex_classification':'deferred until after native simplification',
        'conditional_source_interpretation':True,
    })
    simplify_prepared(lower=lower,before=before)


def simplify_prepared(lower=None,before=None):
    require((OUT/'truncated_raw.rga').exists(),'saved ideal truncation required')
    if lower is None:
        lower=load4(OUT/'truncated_raw.rga')
    if before is None:
        before={'pentachora':1488,'vertices':12,'ideal_boundary_components':1}
    truncated_size=lower.size()
    now=time.monotonic()
    print(json.dumps({'phase':'finite_native_simplify_started','pentachora':lower.size()}),flush=True)
    native=bounded_native_simplify(lower)
    elapsed=time.monotonic()-now
    save4(lower,OUT/'lower_reference.rga')
    print(json.dumps({'phase':'finite_native_simplify_saved','pentachora':lower.size(),'elapsed':elapsed}),flush=True)
    require(lower.isValid() and not lower.isIdeal() and lower.countBoundaryComponents()==1 and lower.boundaryComponent(0).isReal(), 'simplified truncation is not finite with one real boundary')
    boundary,locations=boundary_data(lower)
    save3(boundary,OUT/'lower_boundary.rga')
    save_json(OUT/'lower_boundary_locations.json',locations)
    report={
        'status':'PASS_FINITE_PREPARATION_CONDITIONAL',
        'raw':before,'truncated_pentachora':truncated_size,
        'operation':'Regina Triangulation4.simplify()',
        'source_raw_sha256':sha(RAW),
        'source_truncated_sha256':sha(OUT/'truncated_raw.rga'),
        'result_reference_sha256':sha(OUT/'lower_reference.rga'),
        'native_regina_simplification_used':True,
        'individual_move_log_available':False,
        'native_simplify_seconds':elapsed,
        'native_progress_tracker':native,
        'regina_python_distribution':importlib.metadata.version('regina'),
        'pentachora':lower.size(),'vertices':lower.countVertices(),
        'isosig':lower.isoSig(),'boundary_tetrahedra':boundary.size(),
        'boundary_isosig':boundary.isoSig(),
        'conditional_source_interpretation':True,'standardness_asserted':False,
    }
    save_json(OUT/'LOWER_REFERENCE.json',report)
    print(json.dumps({k:report[k] for k in ['status','pentachora','boundary_tetrahedra','boundary_isosig','native_simplify_seconds']}),flush=True)


def ideal_precondition():
    global OUT
    parent=OUT
    require((parent/'CONDITIONAL_SCOPE.json').exists(),'bound parent scope required')
    OUT=parent/'ideal_preconditioned'
    require(not (OUT/'simplified_ideal.rga').exists(),'completed ideal-preconditioning branch already exists')
    OUT.mkdir(exist_ok=True)
    scope=json.loads((parent/'CONDITIONAL_SCOPE.json').read_text())
    scope['additional_native_trust_boundary']='Regina ideal-body simplify before truncation'
    scope['parent_scope_sha256']=sha(parent/'CONDITIONAL_SCOPE.json')
    save_json(OUT/'CONDITIONAL_SCOPE.json',scope)
    lower=load4(RAW)
    now=time.monotonic()
    print(json.dumps({'phase':'ideal_native_simplify_started','pentachora':lower.size()}),flush=True)
    native=bounded_native_simplify(lower)
    elapsed=time.monotonic()-now
    save4(lower,OUT/'simplified_ideal.rga')
    save_json(OUT/'IDEAL_PRECONDITIONING.json',{
        'source_raw_sha256':sha(RAW),'simplified_ideal_sha256':sha(OUT/'simplified_ideal.rga'),
        'operation':'Regina Triangulation4.simplify() on ideal lower body',
        'native_simplify_seconds':elapsed,'pentachora':lower.size(),
        'native_progress_tracker':native,
        'individual_move_log_available':False,'conditional_source_interpretation':True,
        'source_confirmation_received':False,'standardness_asserted':False,
    })
    print(json.dumps({'phase':'ideal_native_simplify_saved','pentachora':lower.size(),'elapsed':elapsed}),flush=True)
    require(lower.truncateIdeal(),'preconditioned ideal truncation made no change')
    save4(lower,OUT/'truncated_raw.rga')
    print(json.dumps({'phase':'preconditioned_truncation_saved','pentachora':lower.size()}),flush=True)
    simplify_prepared(lower=lower)


def find_path(seconds,seed):
    require((OUT/'LOWER_REFERENCE.json').exists(),'prepare first')
    start=json.loads((OUT/'LOWER_REFERENCE.json').read_text())['boundary_isosig']
    return search(start,OUT/f'boundary_search_seed{seed}',seconds=seconds,seed=seed,goals={GOAL})


def finish(path_file):
    path_file=Path(path_file)
    path=json.loads(path_file.read_text())
    require(path['end']==GOAL,'path does not end at the reference upper-cap boundary')
    replay_boundary(path)
    shutil.copy2(path_file,OUT/'selected_boundary_path.json')
    layer_path(OUT/'lower_reference.rga',OUT/'selected_boundary_path.json',OUT/'lower_layered')
    report=glue_pieces(OUT/'lower_layered/cap_with_matched_boundary.rga',CAP,OUT/'closed')
    report['conditional_source_interpretation']=True
    report['source_confirmation_received']=False
    report['scope_sha256']=sha(OUT/'CONDITIONAL_SCOPE.json')
    report['isotopy_to_previous_Figure2_asserted']=False
    report['standardness_asserted']=False
    save_json(OUT/'CONDITIONAL_COMPLETION.json',report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--prepare',action='store_true')
    choice.add_argument('--simplify-prepared',action='store_true')
    choice.add_argument('--ideal-precondition',action='store_true')
    choice.add_argument('--search',action='store_true')
    choice.add_argument('--finish',type=Path)
    parser.add_argument('--seconds',type=float,default=120)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--out',type=Path,default=OUT)
    parser.add_argument('--native-seconds',type=float,default=NATIVE_SECONDS)
    a=parser.parse_args()
    OUT=a.out.resolve()
    NATIVE_SECONDS=a.native_seconds
    if a.prepare:prepare()
    elif a.simplify_prepared:simplify_prepared()
    elif a.ideal_precondition:ideal_precondition()
    elif a.search:find_path(a.seconds,a.seed)
    else:finish(a.finish)

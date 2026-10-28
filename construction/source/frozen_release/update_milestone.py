#!/usr/bin/env python3
"""Combine separately verified construction, invariant and standardization gates."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import argparse,json,re
from torus_adaptor import load,digest
ROOT=Path(__file__).resolve().parent
def fh(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def req(x,m):
    if not x:raise ValueError(m)
def main(search_finished=False):
    g=load(ROOT/'control/GEOMETRIC_CONSTRUCTION_ADMISSION.json')
    q=load(ROOT/'group_audit/independent_replay.json')
    c=load(ROOT/'unmarked_contractions_v1/independent_replay.json')
    m=load(ROOT/'control/marked_control.json.gz')
    req(g['result']=='PASS','geometric gate')
    for row in g['integrity']['input_artifacts']:
        req(fh(row['path'])==row['file_sha256'],'geometric gate artifact changed: '+row['path'])
    req(fh('control/marked_control.json.gz')==g['immutable_marked_master']['file_sha256']==q['source_file_sha256'],'master binding')
    req(q['passed'] and q['fundamental_group']=='trivial' and q['remaining_generators']==q['remaining_relators']==0,'group gate')
    req(q['source_facets_sha256']==m['facets_sha256']==c['source_facets_sha256'],'group/prefix/master facet binding')
    req(c['passed'] and c['endpoint_native_file_sha256']==fh('unmarked_contractions_v1/endpoint_regina.json.gz'),'contraction prefix gate')
    standard=[];partial=[]
    for p in sorted((ROOT/'control').glob('standardization_seed*_replay.json')):
        if not re.fullmatch(r'standardization_seed[0-9]+_replay\.json',p.name):continue
        r=load(p);seed=int(p.name.split('seed')[1].split('_')[0]);cp=f'control/standardization_seed{seed}.json';cert=load(ROOT/cp)
        req(r['certificate_sha256']==fh(cp),'standardization certificate changed after replay')
        req(r['source_sha256']==cert['source_sha256']==c['endpoint_native_file_sha256'],'standardization source binding')
        req(cert['source_facets_sha256']==m['facets_sha256'],'standardization master binding')
        req(r['replay_passed'],'native replay failed')
        req(r['blocks']==len(cert['blocks']) and r['blocks']>0,'native replay block count')
        req(r['guarded_kicks']==sum(len(b['kicks']) for b in cert['blocks']),'native replay kick count')
        req(r['endpoint_pentachora']==cert['final_size']==cert['blocks'][-1]['pentachora'],'native replay endpoint size')
        req(r.get('endpoint_numbered_sha256')==cert['blocks'][-1]['hash'] and r.get('endpoint_tight_encoding_match') is True,'native replay numbered endpoint')
        req(r.get('source_facets_sha256')==m['facets_sha256'] and r.get('puncture')==cert['puncture'] and r.get('source_puncture_vertices')==cert['source_puncture_vertices'],'native replay explicit puncture binding')
        item={'seed':seed,'replay':str(p.relative_to(ROOT)),'replay_sha256':fh(p.relative_to(ROOT)),
              'certificate':cp,'certificate_sha256':fh(cp),'endpoint_pentachora':r['endpoint_pentachora'],
              'endpoint_isoSig':r['endpoint_isoSig'],'blocks':r['blocks'],'guarded_kicks':r['guarded_kicks']}
        if r['standardization_passed']:
            req(r['ball'] and r['endpoint_pentachora']==1 and r['endpoint_isoSig']=='baa','reference ball endpoint')
            standard.append(item)
        else:partial.append(item)
    status='STANDARDIZED_X2' if standard else ('CONTROL_NOT_STANDARDIZED_WITHIN_DECLARED_BATCH' if search_finished else 'CONSTRUCTION_AND_INVARIANTS_PASSED_STANDARDIZATION_RUNNING')
    result={'schema':'marked-control-milestone-status-v1','updated_utc':datetime.now(timezone.utc).isoformat(),
        'status':status,'parameters':{'knot':'P(-2,3,7)','roll':1,'twist':2},
        'homotopy_type':'S4','homotopy_type_basis':'The source-bound trivial fundamental group and integral homology of S4 imply homotopy equivalence to S4 by Hurewicz and the homological Whitehead theorem. This is not a PL standardization.',
        'geometric_construction':{'passed':True,'receipt':'control/GEOMETRIC_CONSTRUCTION_ADMISSION.json','receipt_sha256':fh('control/GEOMETRIC_CONSTRUCTION_ADMISSION.json'),
            'trust_basis':'Exact finite geometry and signed markings, pinned native PL transport/recognition, and the reviewed relative PL normal-neighborhood and boundary-extension theorem. A finite global normal-bundle matrix atlas is not claimed.'},
        'fundamental_group':{'value':'trivial','independent_source_replay':True,'receipt':'group_audit/independent_replay.json','receipt_sha256':fh('group_audit/independent_replay.json')},
        'integral_homology':q['integral_homology'],'homology_derivation':q['homology_derivation'],
        'marked_master':g['immutable_marked_master'],
        'certified_unmarked_closed_prefix':{'pentachora':c['endpoint_pentachora'],'contractions':c['contractions'],'receipt':'unmarked_contractions_v1/independent_replay.json','receipt_sha256':fh('unmarked_contractions_v1/independent_replay.json')},
        'standardization':{'passed':bool(standard),'successful_replays':standard,'partial_replays':partial,'search_finished_for_this_batch':search_finished,
            'acceptance_rule':'Immutable source -> independently replayed closed PL contraction prefix -> embedded simplex puncture -> all guarded native macros replayed with exact hashes -> isomorphism to Example4.ball(), isoSig baa.'},
        'closed_X0_constructed':False,'exoticness_inference':None,
        'scope':'Homotopy-sphere invariants and faithful geometric identification are separate from a standardization certificate. Failed or bounded heuristic search implies no exoticness claim.'}
    p=ROOT/'control/MILESTONE_STATUS.json';p.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':status,'master_sha256':g['immutable_marked_master']['file_sha256'],'standardization_successes':len(standard),'partial_replays':len(partial)}))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--search-finished',action='store_true');main(a.parse_args().search_finished)

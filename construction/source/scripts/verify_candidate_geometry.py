#!/usr/bin/env python3
"""Bind independently checked X0 assembly to the unchanged admitted source.

All frozen bytes are verified. The m=0 specialization uses the u-only roll
transfer and identity spherical filling; no m=2 twist localization is used
as a substitute for checking the actual zero-twist attachment.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, time

ROOT=Path(__file__).resolve().parents[1]
FROZEN=ROOT/'frozen_release'
PINNED_MANIFEST='b05c8a1c7af0c0dcf3092fc1d8492fd813aa7e831e487ef44c9116be4a869fe9'
PINNED_MASTER='af21a7dd107b765d5fef0a257dcb4df6d6df8689c594f584e4318d15c49aa709'


def require(ok, why):
    if not ok:raise ValueError(why)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def load(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)


def main(output=None):
    started=time.monotonic()
    require(sha(FROZEN/'MANIFEST_SHA256.json')==PINNED_MANIFEST,'Wrong frozen source release manifest')
    manifest=load(FROZEN/'MANIFEST_SHA256.json')
    for name,row in manifest['files'].items():
        p=FROZEN/name
        require(p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'Frozen source file changed: '+name)
    require(sha(FROZEN/'control/marked_control.json.gz')==PINNED_MASTER,'Wrong X2 original master')
    gate=load(ROOT/'x2_acceptance.json')
    require(gate.get('gate_open_for_X0') is True and gate.get('original_master_file_sha256')==PINNED_MASTER,'Required X2 control gate not open')
    old_gate=load(FROZEN/'control/GEOMETRIC_CONSTRUCTION_ADMISSION.json')
    require(old_gate['result']=='PASS' and old_gate['immutable_marked_master']['file_sha256']==PINNED_MASTER,'Original source construction not admitted')
    for row in old_gate['integrity']['input_artifacts']:
        require(sha(FROZEN/row['path'])==row['file_sha256'],'Prior geometric proof source changed')

    m0=load(ROOT/'candidate/marked_candidate.json.gz')
    m2=load(FROZEN/'control/marked_control.json.gz')
    a=load(ROOT/'candidate/assembly_verification.json')
    require(a['passed'] and a['source_model_canonical_sha256']==digest(m0) and a['source_facets_sha256']==m0['facets_sha256'],'X0 assembly replay not current')
    for name,h in a['file_sha256'].items():
        path=(ROOT if name.startswith(('candidate/','scripts/')) else FROZEN)/name
        require(sha(path)==h,'X0 assembly bound input changed: '+name)
    require(m0['twist']==0 and m0['roll']==1 and m0['grid']==[3,4,4],'Wrong X0 parameters')
    require(m0['vertices']==m2['vertices'] and m0['deck']==m2['deck'],'Shared source vertex/deck correspondence changed')
    for part in ('N','E','central','bridge'):
        require(m0['parts'][part]==m2['parts'][part],'Changed noncap piece '+part)
    for key in ('adaptor','central','bridge'):
        require(m0['input_bindings'][key]==m2['input_bindings'][key],'Changed shared input '+key)
    for key in ('bulk_boundary','source_port_map','cap_port_map','orientation_anchor'):
        require(m0[key]==m2[key],'Changed source/boundary marking '+key)
    c0=load(FROZEN/'cap/data/cap_m0.json.gz')
    c2=load(FROZEN/'cap/data/cap_m2.json.gz')
    require(m0['input_bindings']['cap']==digest(c0),'X0 cap data not bound')
    require(c0['twist']==0 and c0['beta']['cap_full_turns']==0 and c0['beta']['downstairs_winding']==0 and c0['beta']['unwrapped_displacement']==[0,4],'Actual cap is not zero winding')
    require(c0['boundary']==c2['boundary'] and c0['vertices']==c2['vertices'] and c0['deck']==c2['deck'],'Caps do not have the literal shared attachment mesh')
    require(c0['facets']!=c2['facets'] and c0['boundary_to_source_refinement']!=c2['boundary_to_source_refinement'],'m0 attachment was merely relabelled m2')
    fresh=load(ROOT/'candidate/cap_recheck.json')
    require(fresh['result']=='PASS' and fresh['cases']['0']['result']=='PASS','Current actual cap verification missing')
    require(fresh['cases']['0']['unwrapped_displacement']==[0,4] and fresh['cases']['0']['downstairs_winding']==0,'Current cap verification did not establish zero lift')
    require(all(row['rejected'] for row in fresh['negative_controls']),'Cap corruption controls did not reject')
    zero_path=ROOT/'review/zero_cap_geometry.json'
    zero=load(zero_path)
    require(zero['result']=='PASS' and zero['checker_sha256']==sha(ROOT/'review/check_zero_cap_geometry.py'),'Independent zero-cap checker binding changed')
    require(zero['candidate_comparison']['result']=='PASS' and zero['candidate_comparison']['candidate_file_sha256']==sha(ROOT/'candidate/marked_candidate.json.gz'),'Independent zero-cap review bound to another candidate')
    zero_files={'bridge':'cap/data/source_cap_bridge.json.gz','central_certificate':'central_roll/certificate.json.gz','central_receipt':'central_roll/independent_verification.json','m0':'cap/data/cap_m0.json.gz','m2':'cap/data/cap_m2.json.gz','shared':'cap/data/shared_boundary.json.gz'}
    require(set(zero['zero_cap']['input_file_sha256'])==set(zero_files),'Independent zero-cap input list incomplete')
    for name,path in zero_files.items():
        require(zero['zero_cap']['input_file_sha256'][name]==sha(FROZEN/path),'Independent zero-cap input changed: '+path)
    require(zero['zero_cap']['identity_spherical_tetrahedra']==1792 and zero['zero_cap']['F0']=='identity' and zero['zero_cap']['Phi0']=='identity','Full spherical zero-twist identity not established')
    central=load(FROZEN/'central_roll/independent_verification.json')
    require(central['passed'] and central['model_sha256']==m0['input_bindings']['central'] and central['straightening']['full_integer_winding']==1,'Actual positive roll not certified')
    f0=set(map(tuple,m0['facets']));f2=set(map(tuple,m2['facets']))
    s0={tuple(f):s for f,s in zip(m0['facets'],m0['orientation'])}
    s2={tuple(f):s for f,s in zip(m2['facets'],m2['orientation'])}
    require(all(s0[f]==s2[f] for f in f0&f2),'Shared source orientation changed')
    require(m0['ramification_triangles']==m2['ramification_triangles'] and m0['twin_triangles']==m2['twin_triangles'],'Marked spheres do not use the shared source labels')
    native={}
    for which in ('cover','quotient'):
        row=load(ROOT/'candidate'/f'{which}_native_admission.json')
        require(row['passed'] and row['source_sha256']==sha(ROOT/'candidate'/f'{which}_regina.json.gz'),'Current native manifold admission missing')
        native[which]=row
    result={
        'schema':'marked-X0-geometric-construction-admission-v1','result':'PASS',
        'classification':'GEOMETRICALLY_IDENTIFIED_MARKED_X0_M0_N1',
        'parameters':{'knot':'oriented P(-2,3,7)','twist':0,'roll':1,'upstairs_source_winding':[1,0],'downstairs_source_winding':[1,0],'cap_unwrapped_displacement':[0,4]},
        'source_control_gate':{'path':'x2_acceptance.json','file_sha256':sha(ROOT/'x2_acceptance.json'),'original_master_file_sha256':PINNED_MASTER,'gate_open_for_X0':True},
        'source_release':{'manifest_file_sha256':PINNED_MANIFEST,'payload_files_checked':len(manifest['files']),'original_geometry_admission_file_sha256':sha(FROZEN/'control/GEOMETRIC_CONSTRUCTION_ADMISSION.json'),'prior_geometry_input_files_checked':len(old_gate['integrity']['input_artifacts'])},
        'immutable_marked_master':{'path':'candidate/marked_candidate.json.gz','file_sha256':sha(ROOT/'candidate/marked_candidate.json.gz'),'canonical_sha256':digest(m0),'facets_sha256':m0['facets_sha256'],'pentachora':len(f0),'vertices':len(m0['vertices']),'parts':{k:len(v) for k,v in m0['parts'].items()}},
        'shared_source_comparison':{'identical_noncap_pentachora':sum(len(m0['parts'][k]) for k in ('N','E','central','bridge')),'common_total_pentachora':len(f0&f2),'removed_m2_cap_pentachora':len(f2-f0),'inserted_m0_cap_pentachora':len(f0-f2),'identical_vertices':True,'identical_deck':True,'identical_marked_surface_labels':True,'identical_source_and_cap_port_maps':True,'orientation_on_all_shared_cells':True,'source_port_tetrahedra':a['source_port_tetrahedra'],'cap_port_tetrahedra':a['cap_port_tetrahedra']},
        'actual_assembly':{'receipt':'candidate/assembly_verification.json','file_sha256':sha(ROOT/'candidate/assembly_verification.json'),'all_cells_regenerated':True,'all_covering_face_squares_commute':a['all_covering_face_squares_commute'],'full_fixed_subcomplex_checked':a['full_fixed_subcomplex_checked'],'normal_branch_circles_checked':a['normal_branch_circles_checked'],'ramification_f_vector':a['ramification_f_vector'],'twin_f_vector':a['twin_f_vector']},
        'zero_cap_direct_proof':{'receipt':'review/zero_cap_geometry.json','receipt_file_sha256':sha(zero_path),'checker_file_sha256':zero['checker_sha256'],'identity_spherical_tetrahedra':1792,'F0':'identity','Phi0':'identity','actual_cap_recheck_file_sha256':sha(ROOT/'candidate/cap_recheck.json')},
        'geometric_identification':[
            'The exact oriented pretzel source, preferred meridian/longitude, marked puncture and deck cover are inherited byte-for-byte from the admitted X2 construction.',
            'The finite u-only source-profile transfer and central-roll certificate identify the unchanged bulk with the relative mapping torus of positive roll n=1 and zero meridional twist. The map fixes the puncture and preserves its marking.',
            'The unchanged source-cap bridge identifies the literal 96-tetrahedron source product boundary with the shared 1792-tetrahedron cap boundary by checked relative refinement collars.',
            'The current independent cap replay verifies the actual m0 rational barycentric attachment, zero lifted angle displacement [0,4], all facet maps and deck compatibility. This is the identity spherical filling class; F0 and Phi0 can be the identity.',
            'All-strata local branch models follow from the admitted source pair, products, equivariant relative ball replacements, unchanged bridge collars and the independently checked m0 cap cone pair, at the actual independently replayed assembly ports.'
        ],
        'trust_basis':[
            'The inherited signed-diagram identification retains the originally declared native Spherogram/SnapPy peripheral transport semantics, with the exact oriented endpoint and both peripheral classes checked.',
            'The X0 cover and quotient local PL-manifold admission uses Regina sphere-link recognition; X2 reviewer link certificates are not misrepresented as X0 link certificates.',
            'The elementary product, refinement collar and cone-pair geometry and the standard relative PL gluing principles retain their mathematical meanings. The extra m2 torus-to-puncture twist-localization theorem is not required for this m0 specialization.'
        ],
        'excluded_conclusions':{'PL_standardization_of_X0':'NOT_ATTEMPTED','global_fundamental_group':'SEPARATE_DIRECT_SOURCE_REPLAY','global_homology':'SEPARATE_DIRECT_SOURCE_REPLAY','exoticness':'NO_INFERENCE'},
        'wall_seconds':time.monotonic()-started}
    output=output or ROOT/'diagnostics/current_geometric_admission.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    main(parser.parse_args().output)

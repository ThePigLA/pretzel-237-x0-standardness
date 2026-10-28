#!/usr/bin/env python3
"""Bind current geometric certificates into the marked m=2,n=1 admission.

Read-only on every source, mesh and component certificate.  This checks
current byte and canonical-object hashes and the inter-component identities;
it writes only the two GEOMETRIC_CONSTRUCTION_ADMISSION files under control.
It neither runs nor admits group, homology or standardization calculations.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import gzip,hashlib,json

ROOT=Path(__file__).resolve().parent
CACHE={};FILES={};LINKS=[];CHECKS=[]
def req(ok,msg):
    if not ok:raise ValueError(msg)
    CHECKS.append(msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def digest(x):return sha(json.dumps(x,sort_keys=True,separators=(',',':')).encode())
def filehash(path):
    path=str(path)
    if path not in FILES:
        b=(ROOT/path).read_bytes();FILES[path]={'path':path,'bytes':len(b),'file_sha256':sha(b)}
    return FILES[path]['file_sha256']
def load(path):
    path=str(path)
    if path not in CACHE:
        b=(ROOT/path).read_bytes();filehash(path)
        CACHE[path]=json.loads(gzip.decompress(b) if path.endswith('.gz') else b)
        FILES[path]['canonical_sha256']=digest(CACHE[path])
    return CACHE[path]
def fb(expected,path,why):
    actual=filehash(path);req(expected==actual,why)
    LINKS.append({'kind':'file_sha256','relation':why,'target':str(path),'sha256':actual})
def cb(expected,path,why):
    actual=digest(load(path));req(expected==actual,why)
    LINKS.append({'kind':'canonical_sha256','relation':why,'target':str(path),'sha256':actual})
def mapped_hashes(row,why,prefix=''):
    for path,h in row.items():fb(h,prefix+path,why+': '+path)

def native_master_export_check(master,cover,quotient):
    """Bind both native 3D tables to actual master facets and deck maps."""
    fs=master['facets'];D=master['deck'];half={}
    req(cover['simplex_count']==len(fs),'Master native cover count')
    for i,k,j,p in cover['gluings']:
        req(sorted(p)==list(range(4)),'Master cover face permutation')
        l=p[k];req((i,k) not in half and (j,l) not in half,'Master cover repeated half facet')
        req(all(fs[i][a]==fs[j][p[a]] for a in range(4) if a!=k),'Master cover face map preserves actual source points')
        half[i,k]=(j,p);ip=[p.index(a) for a in range(4)];half[j,l]=(i,ip)
    req(len(half)==4*len(fs),'Master cover table covers every half facet')
    reps=master['quotient_representatives'];td=master['tetrahedron_deck'];qhalf={}
    req(quotient['simplex_count']==len(reps) and 2*len(reps)==len(fs),'Master quotient count')
    orbit={r:i for i,r in enumerate(reps)};orbit.update({td[r]:i for i,r in enumerate(reps)})
    for i,k,j,p in quotient['gluings']:
        req(sorted(p)==list(range(4)),'Master quotient face permutation');l=p[k]
        req((i,k) not in qhalf and (j,l) not in qhalf,'Master quotient repeated half facet')
        for a,b,face,perm in [(i,j,k,p),(j,i,l,[p.index(z) for z in range(4)])]:
            r=reps[a];dest,cperm=half[r,face];req(orbit[dest]==b,'Master quotient adjacent cell orbit')
            pts=[fs[dest][cperm[z]] for z in range(4)]
            if dest!=reps[b]:pts=[D[x] for x in pts]
            expect=[fs[reps[b]].index(x) for x in pts]
            req(expect==perm,'Master quotient face square commutes with actual deck')
        qhalf[i,k]=(j,p);qhalf[j,l]=(i,[p.index(z) for z in range(4)])
    req(len(qhalf)==4*len(reps),'Master quotient table covers every half facet')
    return {'cover_tetrahedra':len(fs),'quotient_tetrahedra':len(reps),
            'cover_half_facets_checked':len(half),'quotient_half_facets_checked':len(qhalf)}

def admit():
    # Every existing independent audit's complete byte-hash list is checked.
    source_audit=load('diagram/source_audit_certificate.json')
    assembly=load('diagram/control_assembly_verification.json')
    review=load('cap/geometric_review.json');manifest=load('cap/MANIFEST.json')
    req(source_audit['passed'] and assembly['passed'],'Signed source and actual assembly independent gates passed')
    mapped_hashes(source_audit['file_sha256'],'Signed source audit binding')
    mapped_hashes(assembly['file_sha256'],'Actual assembly audit binding')
    for row in review['bindings'].values():fb(row['file_sha256'],row['path'],'Independent localization review binding')
    for row in manifest['files']:
        p='cap/'+row['path'];fb(row['sha256'],p,'Cap package current manifest binding')
        req(FILES[p]['bytes']==row['bytes'],'Cap manifest byte length: '+p)
    req(review['finite_inputs_result']=='PASS' and review['negative_control']['rejected'],'Mixed profile and proof independent review passed')

    # Signed diagram, exact endpoint and native transport ancestry.
    diagram=load('diagram/diagram_verification.json');iso=load('diagram/marked_isomorphism_certificate.json')
    req(diagram['result']=='PASS' and diagram['components']==1 and diagram['crossings']==12,'Actual signed pretzel diagram gate')
    fb(diagram['input_sha256'],'diagram/pretzel_237_diagram.json','Diagram input file binding')
    fb(diagram['checker_sha256'],'diagram/verify_diagram.py','Diagram checker binding')
    req(iso['verified'] and iso['tetrahedra']==3,'Exact marked endpoint isomorphism gate')
    fb(iso['source_file_sha256'],iso['source_file'],'Exact endpoint source file binding')
    fb(iso['target_file_sha256'],iso['target_file'],'Exact endpoint signed-diagram file binding')
    selected=iso['all_isomorphisms'][iso['selected_identity_peripheral_isomorphism']]
    req(selected['orientation_sign']==1 and selected['cusp_matrix_columns_images_of_source_M_L']==[[1,0],[0,1]],'Endpoint preserves ambient orientation and both signed peripheral classes')
    coned=load('peripheral/marked_coned_exterior.json.gz');transport=load('peripheral/marked_exterior_transport.json')
    cb(transport['source_sha256'],'peripheral/marked_coned_exterior.json.gz','Native transport actual marked exterior source')
    cb(coned['cut_sha256'],'peripheral/cut_source.json.gz','Marked exterior cut binding')
    cb(coned['chart_sha256'],'peripheral/peripheral_chart.json.gz','Marked exterior raw chart binding')
    load('peripheral/quotient_exterior_native_path.json')

    # Geometric master, actual descendant and exact orientation/cut ancestry.
    master=load('coarsened/master_source.json.gz');rich=load('seifert/data_oriented/source.json.gz')
    source=load('coarsened/source.json.gz');mo=load('seifert/data/master_orientation.json')
    req(master['facets']==master['closed_cover_facets']==rich['facets']==rich['closed_cover_facets'],'Oriented geometric source and immutable master have identical facets')
    req(master['deck']==rich['deck'] and master['vertex_ranks']==rich['vertex_ranks'],'Oriented geometric source and master have identical marked deck/order')
    req(digest(master['facets'])==master['facets_sha256']==rich['facets_sha256']==mo['facets_sha256'],'Geometric master facet digest')
    req(rich['tetrahedron_orientations']==mo['tetrahedron_orientations'],'Product-derived master orientations')
    cb(source['parent_source_sha256'],'coarsened/master_source.json.gz','Compressed source parent object')
    fb(source['parent_file_sha256'],'coarsened/master_source.json.gz','Compressed source parent bytes')
    cr=load('coarsened/coarsening_verified.json');cc=load('coarsened/coarsening_certificate.json.gz')
    req(cr['passed'] and cr['deck_pair_operations']==len(cc['steps'])==588,'Source relative ball replay count')
    cb(cr['source_sha256'],'coarsened/master_source.json.gz','Source coarsening replay master')
    cb(cr['endpoint_sha256'],'coarsened/source.json.gz','Source coarsening replay endpoint')
    req(cc['initial_sha256']==digest(cc['initial']) and cc['endpoint_sha256']==digest(cc['endpoint']),'Coarsening event-list endpoint digests')
    ori=load('coarsened/source_orientation.json');cut=load('peripheral/cut_source.json.gz');co=load('peripheral/cut_orientation.json')
    opaths={'master':'coarsened/master_source.json.gz','mo':'seifert/data/master_orientation.json','cert':'coarsened/coarsening_certificate.json.gz','source':'coarsened/source.json.gz'}
    for key,path in opaths.items():fb(ori['input_file_sha256'][key],path,'Orientation transport binding '+key)
    req(ori['passed'] and ori['source_facets_sha256']==digest(source['facets']),'Compressed source facet orientation binding')
    cb(cut['source_sha256'],'coarsened/source.json.gz','Actual clipped source binding')
    cb(co['source_sha256'],'coarsened/source.json.gz','Cut determinant source binding')
    cb(co['cut_sha256'],'peripheral/cut_source.json.gz','Cut determinant output binding')
    fb(co['source_orientation_file_sha256'],'coarsened/source_orientation.json','Cut orientation chain binding')
    req(co['passed'] and co['parent_tetrahedra_checked']==len(source['facets'])==628 and co['child_tetrahedra_checked']==len(cut['N'])+len(cut['E'])==1586,'Exact source cut parent/child coverage')
    native3=load('seifert/data/native_admission.json')
    c3=load('seifert/data/cover_regina.json.gz');q3=load('seifert/data/quotient_regina.json.gz')
    master_exports=native_master_export_check(master,c3,q3)
    for which in ('cover','quotient'):
        row=native3['cases'][which];fb(row['input_sha256'],f'seifert/data/{which}_regina.json.gz','Native master admission input '+which)
        req(all(row[k] for k in ('valid','closed','orientable','connected')),'Native master manifold predicates '+which)
    req(native3['cases']['quotient']['S3_recognition'] is True,'Native master quotient three-sphere recognition')

    # Raw and corrected peripheral chart, adaptor, central and mixed profile.
    raw=load('peripheral/peripheral_chart.json.gz');chart=load('peripheral/diagram_chart.json.gz')
    rawgate=load('peripheral/chart_independent_verification.json');iface=load('peripheral/puncture_product_interface.json')
    fb(rawgate['input_file_sha256']['cut'],'peripheral/cut_source.json.gz','Raw chart actual cut bytes')
    fb(rawgate['input_file_sha256']['chart'],'peripheral/peripheral_chart.json.gz','Raw chart data bytes')
    cb(chart['raw_chart_sha256'],'peripheral/peripheral_chart.json.gz','Corrected chart raw ancestry')
    cb(chart['native_transport_sha256'],'peripheral/marked_exterior_transport.json','Corrected chart signed native transport')
    cb(chart['cut_sha256'],'peripheral/cut_source.json.gz','Corrected chart actual cut')
    cb(iface['diagram_chart_sha256'],'peripheral/diagram_chart.json.gz','Puncture product chart binding')
    cb(iface['source_cut_sha256'],'peripheral/cut_source.json.gz','Puncture product source binding')
    req(chart['diagram_sign_transport']['ambient_orientation_multiplier_vs_Seifert_master']==1 and chart['diagram_sign_transport']['boundary_N_vs_du_dv']==-1,'Physical source orientation calibration')
    req(source_audit['diagram_chart_verification']['passed'] and source_audit['core']['N_cocycle_degree']==1,'Source core direction and exact corrected torus chart')
    canonical={'cap_bridge':'cap/data/source_cap_bridge.json.gz','cut':'peripheral/cut_source.json.gz','diagram_chart':'peripheral/diagram_chart.json.gz','puncture_interface':'peripheral/puncture_product_interface.json','raw_chart':'peripheral/peripheral_chart.json.gz','torus_adaptor':'adaptor/adaptor.json.gz','torus_model':'adaptor/model.json.gz'}
    for key,path in canonical.items():cb(source_audit['canonical_digest_links'][key],path,'Signed source canonical link '+key)
    adaptor=load('adaptor/model.json.gz');ac=load('adaptor/adaptor.json.gz');adgate=load('adaptor/independent_verification.json');adco=load('adaptor/coarsening_certificate.json.gz')
    req(adgate['passed'] and adgate['coarsening']['passed'],'Independent adaptor geometry and relative replay')
    cb(adaptor['source_cut_sha256'],'peripheral/cut_source.json.gz','Adaptor model source')
    cb(adaptor['diagram_chart_sha256'],'peripheral/diagram_chart.json.gz','Adaptor model signed chart')
    cb(adaptor['adaptor_sha256'],'adaptor/adaptor.json.gz','Adaptor model geometry')
    cb(adaptor['coarsening_certificate_sha256'],'adaptor/coarsening_certificate.json.gz','Adaptor model coarsening')
    cb(adgate['model_sha256'],'adaptor/model.json.gz','Independent adaptor endpoint')
    cb(ac['cut_sha256'],'peripheral/cut_source.json.gz','Adaptor chart source')
    cb(ac['chart_sha256'],'peripheral/diagram_chart.json.gz','Adaptor chart signed marking')
    cm=load('central_roll/model.json.gz');cent=load('central_roll/certificate.json.gz');cgate=load('central_roll/independent_verification.json')
    cb(cm['central_certificate_sha256'],'central_roll/certificate.json.gz','Central model complete geometry')
    cb(cgate['model_sha256'],'central_roll/model.json.gz','Independent central model')
    cb(cgate['certificate_sha256'],'central_roll/certificate.json.gz','Independent central certificate')
    cpaths={'source_cut_sha256':'peripheral/cut_source.json.gz','diagram_chart_sha256':'peripheral/diagram_chart.json.gz','adaptor_model_sha256':'adaptor/model.json.gz','adaptor_coarsening_certificate_sha256':'adaptor/coarsening_certificate.json.gz'}
    req(cm['source_bindings']==cent['source_bindings'],'Central certificate/model source identities')
    for key,path in cpaths.items():cb(cent['source_bindings'][key],path,'Central source binding '+key)
    req(cgate['passed'] and cgate['straightening']['full_integer_winding']==1 and len(cm['tetrahedra'])==cgate['tetrahedra']==696,'Central positive roll full winding')
    mixed=load('central_roll/mixed_profile_transfer.json');local=load('central_roll/localization_inputs.json')
    for key,path in [('central_certificate_sha256','central_roll/certificate.json.gz'),('adaptor_model_sha256','adaptor/model.json.gz'),('source_cut_sha256','peripheral/cut_source.json.gz')]:cb(mixed[key],path,'Mixed profile binding '+key)
    for key,path in [('source_sha256','coarsened/source.json.gz'),('cut_sha256','peripheral/cut_source.json.gz'),('diagram_chart_sha256','peripheral/diagram_chart.json.gz'),('mixed_profile_transfer_sha256','central_roll/mixed_profile_transfer.json')]:cb(local[key],path,'Localization hypothesis binding '+key)
    req(mixed['passed'] and local['passed'] and mixed['source_winding_upstairs']==[1,1] and mixed['source_winding_downstairs']==[1,2],'Mixed source integer winding and local-flatness input gates')
    req(mixed['enlarged_N']==local['enlarged_N'] and mixed['enlarged_N']['puncture_disjoint_from_outer_grid'],'Corrected enlarged-neighborhood geometry')

    # Finite cap, actual bridge and source-linked local F2/Phi.
    cap=load('cap/data/cap_m2.json.gz');shared=load('cap/data/shared_boundary.json.gz');bridge=load('cap/data/source_cap_bridge.json.gz')
    transfer=load('cap/data/finite_F2_Phi_transfer.json.gz');capgate=load('cap/verification.json');bgate=load('cap/bridge_verification.json');tgate=load('cap/transfer_verification.json')
    req(capgate['result']==bgate['result']==tgate['result']=='PASS','Finite cap, bridge and local transfer gates')
    cb(cap['shared_boundary_sha256'],'cap/data/shared_boundary.json.gz','Cap shared boundary object')
    cb(bridge['cap_shared_canonical_sha256'],'cap/data/shared_boundary.json.gz','Bridge cap boundary object')
    fb(bridge['source_cut_binding']['file_sha256'],'peripheral/cut_source.json.gz','Bridge actual source cut bytes')
    req(bridge['source_cut_binding']['puncture_canonical_sha256']==digest(cut['puncture']),'Bridge actual puncture subrecord')
    fb(bridge['source_interface_file_sha256'],'peripheral/puncture_product_interface.json','Bridge actual puncture product interface')
    req(bridge['source_interface']==iface==transfer['bindings']['source_interface'],'Bridge and local transfer share actual source markings')
    for key,name in [('shared_boundary_file_sha256','shared_boundary.json.gz'),('cap_m2_file_sha256','cap_m2.json.gz'),('source_bridge_file_sha256','source_cap_bridge.json.gz')]:fb(transfer['bindings'][key],'cap/data/'+name,'Finite F2/Phi input '+key)
    req(cap['twist']==2 and cap['beta']['cap_full_turns']==-1 and cap['beta']['unwrapped_displacement']==capgate['cases']['2']['unwrapped_displacement']==tgate['cap_displacement']==[-8,4],'Actual beta2 full integer cap path')
    req(transfer['lifts']['constant_lattice_adjustment']==transfer['lifts']['bulk_D2_unwrapped_shift']==1,'Retained full positive bulk lift')
    req(len(transfer['gamma']['cells'])==tgate['gamma_tetrahedra']==tgate['gamma_radius_identities']==1792 and len(transfer['Phi']['cells'])==tgate['Phi_affine_4_cells']==3584,'Local family and pullback coverage counts')
    req(all(tgate[k] for k in ('R0_R1_identity','deck_equivariant','initial_time_and_bulk_identity','puncture_restriction_is_actual_beta2','unwrapped_D2_Phi1_equals_F2')),'Local collar inverse, boundary, deck and seam predicates')

    # Actual marked four-dimensional model and both current native exports.
    m=load('control/marked_control.json.gz');cover=load('control/cover_regina.json.gz');quotient=load('control/quotient_regina.json.gz')
    ar=load('control/assembly_receipt.json');fs=m['facets'];master_before=filehash('control/marked_control.json.gz')
    req((m['twist'],m['roll'],m['grid'])==(2,1,[3,4,4]),'Current marked model parameters')
    paths={'adaptor':'adaptor/model.json.gz','central':'central_roll/model.json.gz','bridge':'cap/data/source_cap_bridge.json.gz','cap':'cap/data/cap_m2.json.gz'}
    req(m['input_bindings']==assembly['input_bindings'],'Current model and independent assembly input identities')
    for key,path in paths.items():cb(m['input_bindings'][key],path,'Actual marked assembly input '+key)
    cb(assembly['source_model_canonical_sha256'],'control/marked_control.json.gz','Independent assembly exact marked model')
    facet_hash=digest(fs);req(facet_hash==m['facets_sha256']==cover['source_facets_sha256']==quotient['source_facets_sha256']==assembly['source_facets_sha256']==ar['facet_sha256'],'All current four-dimensional artifacts share the same master facets')
    req(len(fs)==cover['simplex_count']==ar['pentachora']==103744 and quotient['simplex_count']==ar['quotient_pentachora']==51872,'Current cover and quotient cell counts')
    counts={k:len(v) for k,v in m['parts'].items()};req(counts==assembly['parts_regenerated']==ar['parts'],'Actual assembled part counts')
    req(assembly['source_port_tetrahedra']==bgate['ports']['source']==96 and assembly['cap_port_tetrahedra']==bgate['ports']['cap']==1792,'Actual source and cap ports close reviewed branch models')
    req(all(assembly[k] for k in ('all_covering_face_squares_commute','deck_orientation_preserving','full_fixed_subcomplex_checked','no_simplex_inversions')) and assembly['omitted_fixed_vertices_or_edges']==0,'Actual covering, deck and fixed-subcomplex gates')
    req(assembly['normal_branch_circles_checked']==208 and len(assembly['source_seam_carriers'])==10,'Actual normal-circle and source-carrier membership checks')
    cb(m['orientation_anchor']['cut_orientation_sha256'],'peripheral/cut_orientation.json','Actual four-dimensional source orientation anchor')
    native={}
    for name,n in [('cover',103744),('quotient',51872)]:
        row=load(f'control/{name}_native_admission.json');fb(row['source_sha256'],f'control/{name}_regina.json.gz','Current native '+name+' admission')
        req(all(row[k] for k in ('passed','valid','closed','orientable','connected')) and not row['ideal'] and row['boundary_components']==0 and row['simplex_count']==n,'Native '+name+' closed nonideal manifold predicates')
        native[name]={k:row[k] for k in ('regina','simplex_count','valid','closed','orientable','connected','ideal','boundary_components','source_sha256','vertex_link_scope')}

    finite=[
        {'stage':'signed diagram and actual source','status':'PASS','basis':'explicit Montesinos cells, complete master export/deck squares, source-relative ball replay, exact marked cusp endpoint and orientation calibration'},
        {'stage':'source cut and peripheral data','status':'PASS','basis':'exact barycentric child partitions and orientation transport; rational torus chart, signed meridian/longitude, deck translation and double-cover projection'},
        {'stage':'torus adaptors and central roll','status':'PASS','basis':'exact overlay/prism charts and inverses, actual ordered ports, relative equivariant coarsening, full roll winding and ten actual 4D seam carriers'},
        {'stage':'mixed profile and local branch inputs','status':'PASS','basis':'six stored simultaneous u/v affine maps/inverses; true enlarged N-plus; complete fixed-circle normal and vertex links; source-positive puncture'},
        {'stage':'actual cap bridge and beta2','status':'PASS','basis':'complete finite refinement collars and cap attachment; full unwrapped displacement [-8,4]; actual 96/1792 ports'},
        {'stage':'local F2/Phi','status':'PASS','basis':'1792 joint-PL disk-family cells and 3584 affine pullback 4-polytopes, inverses, radial identities, seams and deck; canonical simplicial expansion specified but not materialized'},
        {'stage':'actual assembly and quotient','status':'PASS','basis':'complete regeneration of every part; all face pairings and covering squares; source orientation; full fixed subcomplex and no simplex inversions'}]
    trust=[
        {'kind':'native_PL_move_semantics','statement':'Spherogram/SnapPy signed-diagram exterior construction, finite-vertex removal and simplification transport use pinned native PL semantics. Exact source/import signs and the final +I peripheral isomorphism are independently checked. A complete elementary move-by-move peripheral ancestry is not exported.'},
        {'kind':'native_manifold_predicates','statement':'The current cover and quotient native exports satisfy Regina validity, closedness, orientability, connectedness and nonideal predicates, bound by exact file hashes. These predicates are not a standardization certificate.'},
        {'kind':'relative_PL_normal_neighborhood_theorem','statement':'After the finite simultaneous profile transfer, enlarged N-plus admits a framed pair chart agreeing with its outer marked grid and the interior puncture normal chart. Relative normal-neighborhood uniqueness, the solid-torus meridian extension and its branched lift are theorem-level steps; a global finite source-cell atlas is not exported.',
         'proof':'central_roll/LOCALIZATION_LEMMA.txt','proof_file_sha256':filehash('central_roll/LOCALIZATION_LEMMA.txt'),
         'independent_review':'cap/geometric_review.json','review_file_sha256':filehash('cap/geometric_review.json')},
        {'kind':'all_strata_branch_geometry','statement':'The standard local double-branch pair is supplied at every stratum by source normal links, products, equivariant relative ball replacements, bridge collars, and the cap cone pair. This reviewed geometric argument is joined to the actual checked ports and full assembly, rather than inferred from fixed-triangle normal circles alone.',
         'review':'cap/REVIEW_BRANCH_AND_LOCALIZATION.txt','review_file_sha256':filehash('cap/REVIEW_BRANCH_AND_LOCALIZATION.txt')}]

    # Bind the directly relevant replay/geometry code even where an older
    # component receipt did not carry its own checker-byte field.
    for path in ('seifert/build_source.py','seifert/native_admit.py',
                 'central_roll.py','verify_central_roll.py','verify_mixed_profile_transfer.py',
                 'verify_twist_localization_inputs.py','torus_adaptor.py',
                 'lib/exact_flat_mesh.py','lib/verify_adaptor_mesh.py','lib/compact3.py',
                 'lib/verify_compact3.py','lib/branch_dpll.py','geometric_construction_admission.py'):
        filehash(path)
    # Detect concurrent changes before issuing the gate. No source file is written.
    for path,row in list(FILES.items()):req(sha((ROOT/path).read_bytes())==row['file_sha256'],'Input unchanged during admission: '+path)
    req(filehash('control/marked_control.json.gz')==master_before,'Immutable marked master preserved')
    # De-duplicate semantic check labels; repeated native half-facet checks are
    # represented by their exact counts instead of tens of thousands of lines.
    result={'schema':'marked-X2-geometric-construction-admission-v1','result':'PASS',
        'classification':'GEOMETRICALLY_IDENTIFIED_MARKED_X2_M2_N1',
        'meaning':'The actual closed marked control is identified with the intended signed J=P(-2,3,7), m=2, n=1 construction, under the stated native semantics and relative PL theorem trust basis.',
        'parameters':{'knot':'oriented P(-2,3,7)','twist':2,'roll':1,'upstairs_source_winding':[1,1],'downstairs_source_winding':[1,2],'cap_unwrapped_displacement':[-8,4]},
        'immutable_marked_master':{'path':'control/marked_control.json.gz','file_sha256':master_before,'canonical_sha256':digest(m),'facets_sha256':facet_hash,'pentachora':len(fs),'quotient_pentachora':quotient['simplex_count'],'parts':counts,'unchanged':True},
        'finite_component_gates':finite,'master_native_export_direct_check':master_exports,
        'native_closed_admission':native,'geometric_trust_basis':trust,
        'orientation':{'source_ambient_vs_Seifert_master':1,'boundary_N_vs_du_dv':-1,'boundary_E_vs_du_dv':1,'four_dimensional':'positive source orientation followed by positive time'},
        'actual_interfaces':{'source_port_tetrahedra':96,'cap_port_tetrahedra':1792,'normal_branch_triangles_checked':208,'source_roll_seam_carriers':10,'ramification_f_vector':assembly['ramification_f_vector']},
        'excluded_conclusions':{'global_fundamental_group':'NOT_ADMITTED_BY_THIS_GATE','global_homology':'NOT_ADMITTED_BY_THIS_GATE','PL_standardization_of_X2':'NOT_ADMITTED_BY_THIS_GATE','X0_construction':'REMAINS_GUARDED_BY_REQUIRED_X2_STANDARDIZATION','fully_finite_global_N_plus_conjugacy':'NOT_EXPORTED'},
        'integrity':{'mode':'current byte hashes, canonical-object hashes, cross-certificate identities and direct master export/deck checks; existing finite geometry receipts are aggregated without altering them',
            'checked_hash_links':len(LINKS),'bound_artifacts':len(FILES),'semantic_check_count':len(CHECKS),
            'hash_links':LINKS,'input_artifacts':sorted(FILES.values(),key=lambda r:r['path'])},
        'created_utc':datetime.now(timezone.utc).isoformat()}
    out=ROOT/'control/GEOMETRIC_CONSTRUCTION_ADMISSION.json';out.write_text(json.dumps(result,indent=2)+'\n')
    txt=f'''GEOMETRIC CONSTRUCTION ADMISSION: MARKED X2, m=2, n=1
=====================================================

Result: PASS -- GEOMETRICALLY_IDENTIFIED_MARKED_X2_M2_N1.

The actual immutable marked model is geometrically identified with the
intended oriented J=P(-2,3,7), positive twist m=2 and roll n=1 construction.
This admission includes the finite geometry, the bound native predicates,
and the explicit theorem dependencies below. It is not an X2 standardization.

Actual model
------------
Cover: 103,744 pentachora. Orbit quotient: 51,872 pentachora.
N product: 25,216; E product: 37,344; central roll: 11,136;
source-cap bridge: 24,224; cap: 5,824.
The source and cap ports have respectively 96 and 1,792 tetrahedra.
All part cells, face pairings, deck/covering squares, fixed faces, and source
orientation are bound to the independent assembly replay. All ten stored
source-roll seam carriers occur in actual master pentachora.

Master file SHA-256: {master_before}
Master facets SHA-256: {facet_hash}
The marked master and all component data were read without modification.

What passed
-----------
The signed source-to-diagram identification, exact source cut, rational
peripheral charts, torus adaptors, relative equivariant coarsenings, central
roll and integer source winding are current and hash-bound. The simultaneous
u/v profile correction is checked, and the actual enlarged N-plus separates
the puncture from the outer grid. The cap retains the full path (-8,4).
The finite local F2/Phi has 1,792 disk-family cells and 3,584 affine pullback
4-polytopes, with inverse, boundary, deck and integer seam certificates.
Both current native exports pass validity, closedness, connectedness,
orientability and nonideal checks.

The geometric branch review covers every stratum using the actual source
pair, products, relative deck-paired ball replacements, refinement collars,
and the cap cone pair. The actual assembly ports satisfy that review's
gluing hypotheses. This conclusion does not use homology to infer branch
models or to recognize the closed cover as a sphere.

Explicit trust basis
--------------------
Native diagram/exterior construction and intermediate legal PL moves retain
their native semantics. Exact endpoint cusp signs and the +I peripheral
identification are checked independently; a complete elementary peripheral
move transcript is not exported.

The global framed normal-disk chart on enlarged N-plus is supplied by the
relative PL regular-neighborhood and marked solid-torus extension theorems.
LOCALIZATION_LEMMA.txt explains its branched lift, the source twist-to-cap
isotopy, and the joint-PL Alexander comparison. Independent review accepted
that corrected argument. Software does not claim to have enumerated the
entire global N-plus conjugacy as a source-cell matrix atlas.

The machine-readable admission records {len(LINKS)} checked hash links,
{len(FILES)} bound artifacts, complete current hashes and the component gates.

Separate unfinished gates
-------------------------
This admission does not mark global fundamental group, global homology or
PL standardization of X2 as passed. X0 remains behind the required
source-bound X2 standardization gate. Native validity and geometric
identification do not replace that standardization certificate.
'''
    (ROOT/'control/GEOMETRIC_CONSTRUCTION_ADMISSION.txt').write_text(txt)
    print(json.dumps({'result':result['result'],'classification':result['classification'],'hash_links':len(LINKS),'artifacts':len(FILES),'master_file_sha256':master_before,'admission_file_sha256':sha(out.read_bytes())},indent=2),flush=True)

if __name__=='__main__':admit()

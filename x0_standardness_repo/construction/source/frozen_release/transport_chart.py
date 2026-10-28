"""Apply the explicitly recorded source-to-diagram peripheral sign transition."""
from copy import deepcopy
from fractions import Fraction
from prepare_source import read,save,ROOT
from group_cells import inv,digest

def build():
    p=ROOT/'peripheral';raw=read(p/'peripheral_chart.json.gz');native=read(p/'marked_exterior_transport.json')
    if native['exact_import_normalization']['meridian_sign']!=-1:raise ValueError('unexpected native normalization')
    if not any(i['cusp_maps']==[[[1,0],[0,1]]] for i in native['exact_combinatorial_isomorphisms']):raise ValueError('identity cusp map absent')
    x=deepcopy(raw);x['schema']='diagram-aligned-pretzel-peripheral-chart-v1'
    x['raw_chart_sha256']=digest(raw);x['native_transport_sha256']=digest(native)
    x['E_cocycle']=[-z for z in raw['E_cocycle']]
    x['E_H1_certificate']['sign_normalization']*=-1
    x['downstairs_periods']=[[u,-v] for u,v in raw['downstairs_periods']]
    x['upstairs_periods']=[[u,-v] for u,v in raw['upstairs_periods']]
    x['downstairs_vertex_coordinates']=[[a,u,str(-Fraction(v))] for a,u,v in raw['downstairs_vertex_coordinates']]
    x['upstairs_vertex_coordinates']=[[a,u,str(-Fraction(v))] for a,u,v in raw['upstairs_vertex_coordinates']]
    for key in ('torus_edge_loop','raw_exterior_word','reduced_exterior_word'):
        x['peripheral_loops']['meridian'][key]=list(inv(raw['peripheral_loops']['meridian'][key]))
    x['oriented_base_meridian_lift']=list(reversed(raw['oriented_base_meridian_lift']))
    x['coordinate_meaning']['v_downstairs']='positive meridian transported to the signed diagram by the recorded native normalization and exact cusp identity'
    x['diagram_sign_transport']={'coordinate_matrix':[[1,0],[0,-1]],'coordinate_order':['longitude','meridian'],
                                'ambient_orientation_multiplier_vs_Seifert_master':1,
                                'imported_ordered_simplex_multiplier_vs_master':-1,
                                'snappy_right_handed_vs_standard_ordered_simplex':-1,
                                'boundary_N_vs_du_dv':-1,
                                'boundary_E_vs_du_dv':1,
                                'ambient_orientation_convention':'positive core longitude followed by positive normal disk; ambient du wedge dr wedge dv, outward boundary N is -du wedge dv',
                                'reason':'source core has N-cocycle degree +1; reversed puncture equator is positive normal disk in the master orientation; imported ordered-simplex sign -1 and SnapPy right-handed convention -1 cancel; exact marked endpoint isomorphism is even and carries M,L with their signs'}
    save(p/'diagram_chart.json.gz',x)
    cut=read(p/'cut_source.json.gz');iface=read(p/'puncture_product_interface.json');iface['source_cut_sha256']=digest(cut)
    iface['diagram_chart_sha256']=digest(x);save(p/'puncture_product_interface.json',iface)
    return x

if __name__=='__main__':
    x=build();print('DIAGRAM CHART',digest(x))

"""Source-bound cusp markings before native finite-vertex removal."""
from collections import defaultdict
from fractions import Fraction as Q
import json,time
from prepare_source import read,save,export3,ROOT
from group_cells import digest

def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def inverse(p):return [p.index(i) for i in range(len(p))]

def source_ideal(cut,chart):
    D=cut['deck'];base=export3(cut['E'],D);reps=list(map(tuple,base['tetrahedron_representatives']));n=len(reps)
    data=[[[-1]*4,[None]*4] for _ in reps];labels=[list(t) for t in reps]
    def join(i,f,j,p):
        g=p[f]
        if data[i][0][f]>=0 or data[j][0][g]>=0:raise ValueError('already paired')
        data[i][0][f]=j;data[i][1][f]=p;data[j][0][g]=i;data[j][1][g]=inverse(p)
    for i,f,j,p in base['gluings']:join(i,f,j,p)
    boundary=[(i,f) for i in range(n) for f in range(4) if data[i][0][f]<0]
    sides=defaultdict(list);cones=[]
    for i,f in boundary:
        t=reps[i];face=[v for k,v in enumerate(t) if k!=f]
        j=len(data);data.append([[-1]*4,[None]*4]);labels.append(face+[-1]);cones.append(j)
        join(j,3,i,[t.index(v) for v in face]+[f])
        for g in range(3):
            e=tuple(sorted(v for k,v in enumerate(face) if k!=g));r=min(e,tuple(sorted(D[v] for v in e)))
            sides[r].append((j,g,e))
    for r,hh in sides.items():
        if len(hh)!=2:raise ValueError('boundary not closed torus')
        (i,f,e),(j,g,h)=hh;direct=e==h;p=[]
        for k,v in enumerate(labels[i]):
            p.append(g if k==f else 3 if v==-1 else labels[j].index(v if direct else D[v]))
        join(i,f,j,p)
    if any(z<0 for t in data for z in t[0]):raise ValueError('unpaired ideal face')
    sk=chart['torus_downstairs'];per=chart['downstairs_periods'];xy={v:(Q(u),Q(w)) for v,u,w in chart['downstairs_vertex_coordinates']}
    ei={tuple(e):i for i,e in enumerate(sk['edge_representatives'])}
    def period(a,b):
        e=tuple(sorted((a,b)));r=min(e,tuple(sorted(D[v] for v in e)));sg=1 if (a,b)==r or (D[a],D[b])==r else -1
        return tuple(sg*z for z in per[ei[r]])
    def position(v):return xy[min(v,D[v])]
    def area(t):
        a,b,c=t;pa=position(a);pb=tuple(x+y for x,y in zip(position(b),period(a,b)));pc=tuple(x+y for x,y in zip(position(c),period(a,c)))
        return (pb[0]-pa[0])*(pc[1]-pa[1])-(pb[1]-pa[1])*(pc[0]-pa[0])
    i,f=boundary[0];ar=area(tuple(v for k,v in enumerate(reps[i]) if k!=f))
    if not ar:raise ValueError('degenerate source chart')
    orientation={i:(1 if ar>0 else -1)*(-1)**f};todo=[i]
    for i in todo:
        for f,j in enumerate(data[i][0]):
            s=-parity(data[i][1][f])*orientation[i]
            if j in orientation:
                if s!=orientation[j]:raise ValueError('nonorientable exterior')
            else:orientation[j]=s;todo.append(j)
    if len(orientation)!=len(data):raise ValueError('disconnected ideal exterior')
    for i,f in boundary:
        if area(tuple(v for k,v in enumerate(reps[i]) if k!=f))*orientation[i]*(-1)**f<=0:raise ValueError('chart orientation inconsistent')
    P=[list(range(4)) if orientation[i]>0 else [1,0,2,3] for i in range(len(data))]
    out=[];newlabels=[];cusp=[];curves=[]
    for i,(neighbors,glue) in enumerate(data):
        newlabels.append([labels[i][k] for k in P[i]]);ne=[];gg=[]
        for f in range(4):
            oldf=P[i][f];j=neighbors[oldf];pj=inverse(P[j])
            ne.append(j);gg.append([pj[glue[oldf][P[i][k]]] for k in range(4)])
        if any(parity(p)!=-1 for p in gg):raise ValueError('not consistently oriented')
        out.append([ne,gg]);cusp.append([0 if v==-1 else -1 for v in newlabels[-1]])
        rows=[[0]*16 for _ in range(4)]
        if i in cones:
            v=newlabels[-1].index(-1);cc=[k for k in range(4) if k!=v]
            if (-1)**(v+1)<0:cc[0],cc[1]=cc[1],cc[0]
            for k,f in enumerate(cc):
                a,b=(newlabels[-1][cc[(k+j)%3]] for j in (1,2));u,w=period(a,b)
                rows[0][4*v+f]=u;rows[2][4*v+f]=-w
        curves.extend(rows)
    return {'schema':'source-marked-coned-exterior-v1','data':out,'cusp_indices':cusp,'curves':curves,
            'tetrahedron_source_vertices':newlabels,'tetrahedron_orientations_before_relabeling':[orientation[i] for i in range(len(data))],
            'positive_relabelings':P,'compact_tetrahedra':n,'cone_tetrahedra':cones,
            'cut_sha256':digest(cut),'chart_sha256':digest(chart),
            'orientation':'outward compact-exterior boundary chart du wedge dv',
            'flow_convention':'positive flow ENTERS the cusp triangle; M=du, L=-dv in the displayed cusp-side orientation; importer read-back audited'}

def run():
    import snappy,spherogram
    start=time.monotonic();out=ROOT/'peripheral';cut=read(out/'cut_source.json.gz');chart=read(out/'peripheral_chart.json.gz')
    x=source_ideal(cut,chart);save(out/'marked_coned_exterior.json.gz',x)
    T=snappy.Triangulation('empty');T._from_tetrahedra_gluing_data(x['data'],num_or_cusps=1,cusp_indices=x['cusp_indices'],peripheral_curves=x['curves'],remove_finite_vertices=False)
    got=T._get_cusp_indices_and_peripheral_curve_data()
    info={'schema':'source-marked-exterior-transport-v1','source_sha256':digest(x),'input_tetrahedra':T.num_tetrahedra(),
          'curve_arrays_unchanged_by_import':got[1]==x['curves'],'cusp_indices_unchanged_by_import':got[0]==x['cusp_indices'],
          'kernel':'SnapPy 3.3.2 / SnapPea; explicit source arrays and trusted native move transport'}
    flipped=[[-z for z in row] if i%4==0 else row for i,row in enumerate(x['curves'])]
    if not info['cusp_indices_unchanged_by_import']:raise ValueError('native importer changed source cusp indices')
    if got[1]==x['curves']:normalization=1
    elif got[1]==flipped:normalization=-1
    else:raise ValueError('unexplained change of source peripheral arrays')
    info['exact_import_normalization']={'meridian_sign':normalization,'longitude_sign':1,
        'description':'literal equality of every source entry after this diagonal sign change; kernel enforces its peripheral intersection convention',
        'native_meridian_periods_in_source_u_v':[0,normalization],'native_longitude_periods_in_source_u_v':[1,0]}
    save(out/'native_import_peripheral_arrays.json.gz',{'cusp_indices':got[0],'curves':got[1],'raw_source_sha256':digest(x),'normalization':normalization})
    (out/'marked_coned_exterior.tri').write_text(T._to_string());print('IMPORTED',info,flush=True)
    T._remove_finite_vertices();T.simplify();print('REMOVED',T.num_tetrahedra(),T.triangulation_isosig(),flush=True)
    (out/'marked_ideal_exterior.tri').write_text(T._to_string())
    info['output_tetrahedra']=T.num_tetrahedra();info['output_decorated_isoSig']=T.triangulation_isosig()
    G=T.fundamental_group(False,False);info['presentation']={'generators':G.generators(),'relators':G.relators(),'peripheral_curves':G.peripheral_curves()}
    import sympy as S,math
    gs=G.generators();mat=S.Matrix([[w.count(g)-w.count(g.upper()) for g in gs] for w in G.relators()]);basis=mat.nullspace()
    if len(basis)!=1:raise ValueError('exterior abelian rank')
    den=math.lcm(*(int(z.q) for z in basis[0]));weights=[int(z*den) for z in basis[0]];gg=math.gcd(*weights);weights=[z//gg for z in weights]
    hom=lambda w:sum((w.count(g)-w.count(g.upper()))*z for g,z in zip(gs,weights))
    mu,lam=G.peripheral_curves()[0]
    if abs(hom(mu))!=1 or hom(lam)!=0:raise ValueError('marked zero longitude or primitive meridian lost')
    info['H1']=str(T.homology());info['peripheral_abelian_check']={'generator_weights':weights,'meridian':hom(mu),'longitude':hom(lam)}
    J=(spherogram.RationalTangle(-1,2)+spherogram.RationalTangle(1,3)+spherogram.RationalTangle(1,7)).numerator_closure().exterior(with_hyperbolic_structure=False)
    (out/'diagram_exterior.tri').write_text(J._to_string());isos=T.isomorphisms_to(J);history=[]
    for seed in range(1,9):
        if isos:break
        snappy.set_rand_seed(seed);T.randomize();T.simplify()
        history.append({'seed':seed,'size':T.num_tetrahedra(),'decorated_isoSig':T.triangulation_isosig()});isos=T.isomorphisms_to(J)
    info['matching_macro_history']=history
    info['exact_combinatorial_isomorphisms']=[{'cusp_images':list(i.cusp_images()),'cusp_maps':[[[int(z) for z in row] for row in M] for M in i.cusp_maps()],
                                             'extends_to_link':i.extends_to_link()} for i in isos]
    (out/'matched_marked_exterior.tri').write_text(T._to_string())
    info['matched_decorated_isoSig']=T.triangulation_isosig();info['diagram_decorated_isoSig']=J.triangulation_isosig()
    F=T.copy();F.dehn_fill((1,0));FG=F.fundamental_group();info['meridian_filling_presentation']={'generators':FG.generators(),'relators':FG.relators()}
    info['seconds']=time.monotonic()-start;save(out/'marked_exterior_transport.json',info);print(json.dumps(info,indent=2),flush=True)
    return info

if __name__=='__main__':run()

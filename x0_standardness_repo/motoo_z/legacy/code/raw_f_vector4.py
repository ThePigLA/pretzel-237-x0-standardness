#!/usr/bin/env python3
"""Count raw local-face orbits of a closed4-dimensional triangulation.

Only raw pentachoron facet gluings are read. Each dimension is processed
with its own compact union-find array, without a Regina4D skeleton or
vertex recognition. The alternating sum is the combinatorial Euler
count; its topological interpretation assumes trivial face stabilizers.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,time
from array import array
from pathlib import Path
import regina

def require(ok,msg):
    if not ok:raise ValueError(msg)

def raw_f_vector(source,out,expected=2):
    source,out=Path(source),Path(out);then=time.monotonic()
    tri=regina.Triangulation4(regina.open(str(source)))
    n=tri.size();counts=[];work=[]
    for dimension in range(4):
        faces=list(itertools.combinations(range(5),dimension+1));length=len(faces)
        number={sum(1<<v for v in row):i for i,row in enumerate(faces)}
        in_facet=[[i for i,row in enumerate(faces) if f not in row] for f in range(5)]
        parent=array('I',range(n*length));rank=bytearray(n*length);components=n*length
        def find(x):
            while parent[x]!=x:
                parent[x]=parent[parent[x]];x=parent[x]
            return x
        pairs=unions=0
        for simplex in tri.simplices():
            i=simplex.index()
            for facet in range(5):
                adjacent=simplex.adjacentSimplex(facet)
                require(adjacent is not None,'source has an unpaired facet')
                perm=simplex.adjacentGluing(facet);j=adjacent.index()
                if (i,facet)>(j,perm[facet]):continue
                pairs+=1
                for local in in_facet[facet]:
                    target=number[sum(1<<perm[v] for v in faces[local])]
                    a=find(length*i+local);b=find(length*j+target)
                    if a==b:continue
                    if rank[a]<rank[b]:a,b=b,a
                    parent[b]=a
                    if rank[a]==rank[b]:rank[a]+=1
                    components-=1;unions+=1
        require(pairs*2==5*n,'facet-pair count is inconsistent with a closed triangulation')
        counts.append(components);work.append({'dimension':dimension,'local_faces':n*length,
                                               'facet_pairs_processed':pairs,'successful_unions':unions})
        print(json.dumps({'dimension':dimension,'face_orbits':components,'seconds':time.monotonic()-then}),flush=True)
    counts.append(n);euler=sum((-1)**i*f for i,f in enumerate(counts))
    require(euler==expected,f'combinatorial Euler count is{euler}, expected{expected}')
    report={'status':'PASS_RAW_CLOSED_F_VECTOR_EXPECTED_EULER','source':str(source),
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'source_tight_sha256':hashlib.sha256(tri.tightEncoding().encode()).hexdigest(),
            'f_vector':counts,'combinatorial_euler_characteristic':euler,'expected_euler_characteristic':expected,
            'method':'union-find of local simplex faces under raw facet pairing permutations',
            'work_by_dimension':work,'four_dimensional_skeleton_used':False,
            'vertex_recognition_performed':False,'face_stabilizers_checked_by_this_tool':False,
            'euler_interpretation':'alternating raw face-orbit count; equals topological Euler characteristic when face stabilizers are trivial',
            'standardness_asserted':False,'seconds':time.monotonic()-then,
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    tmp=out.with_suffix(out.suffix+'.tmp');tmp.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');tmp.replace(out)
    print(json.dumps(report),flush=True);return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('--out',required=True)
    p.add_argument('--expected-euler',type=int,default=2);a=p.parse_args();raw_f_vector(a.source,a.out,a.expected_euler)

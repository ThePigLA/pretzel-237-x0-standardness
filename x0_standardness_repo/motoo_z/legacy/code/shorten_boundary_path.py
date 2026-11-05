#!/usr/bin/env python3
"""Erase isomorphism loops, then reindex a pure labelled Pachner path.

The output has no relabelling operations: when a loop is skipped, each
subsequent original move is transported through a checked combinatorial
isomorphism to the live triangulation, and is again applied with its
Pachner guard enabled.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import regina
from boundary_pachner_fast import digest,encoding,face,mapped_face_index,replay,require,save

def concatenate(paths):
    data=json.loads(Path(paths[0]).read_text())
    data['moves']=list(data['moves'])
    for path in paths[1:]:
        other=json.loads(Path(path).read_text())
        require(data['end_tight']==other['start_tight'],'labelled path segments do not meet exactly')
        data['moves'].extend(other['moves']);data['end_tight']=other['end_tight']
        data['end_tetrahedra']=other['end_tetrahedra']
    return data

def shorten(paths,output,limit=200):
    now=time.monotonic();data=concatenate(paths)
    original=regina.Triangulation3.tightDecoding(data['start_tight'])
    keys={};states={};last={}
    for i in range(len(data['moves'])+1):
        if original.size()<=limit:
            key=original.isoSig();keys[i]=key;states[i]=encoding(original);last[key]=i
        if i==len(data['moves']):break
        m=data['moves'][i]
        require(digest(encoding(original))==m['before_sha256'],f'source digest differs at{i}')
        require(original.pachner(face(original,m['dimension'],m['face'])),f'illegal original move{i}')
        require(digest(encoding(original))==m['after_sha256'],f'target digest differs at{i}')
        if (i+1)%1000==0:
            print(json.dumps({'analysed_moves':i+1,'tetrahedra':original.size(),
                              'seconds':time.monotonic()-now}),flush=True)
    require(encoding(original)==data['end_tight'],'original path endpoint differs')
    target=regina.Triangulation3(original)
    original=regina.Triangulation3.tightDecoding(data['start_tight'])
    live=regina.Triangulation3(original);moves=[];i=0;skips=[]
    while i<len(data['moves']):
        key=keys.get(i)
        if key is not None and last[key]>i:
            j=last[key]
            jumped=regina.Triangulation3.tightDecoding(states[j])
            require(live.isIsomorphicTo(jumped) is not None,'loop endpoints are not isomorphic')
            skips.append({'from_original_step':i,'to_original_step':j,'isosig':key})
            original=jumped;i=j;continue
        m=data['moves'][i];d=m['dimension'];index=m['face']
        if encoding(original)==encoding(live):live_index=index
        else:
            iso=original.isIsomorphicTo(live)
            require(iso is not None,f'transported path lost isomorphism at{i}')
            live_index=mapped_face_index(original,live,iso,d,index)
        before=digest(encoding(live));n=live.size()
        require(live.pachner(face(live,d,live_index)),f'transported guard failed at{i}')
        require(original.pachner(face(original,d,index)),f'original guard failed on second replay at{i}')
        moves.append({'dimension':d,'face':live_index,'before_sha256':before,
                      'after_sha256':digest(encoding(live)),'source_tetrahedra':n,'target_tetrahedra':live.size()})
        i+=1
    require(live.isIsomorphicTo(target) is not None,'shortening changed final isomorphism type')
    result={'schema':'motoo.boundary.labelled-pachner.v1','status':'LOOPS_ERASED_AND_PURE_PACHNER_REPLAYED',
            'source_paths':list(map(str,paths)),'start_tight':data['start_tight'],'end_tight':encoding(live),
            'start_tetrahedra':data['start_tetrahedra'],'end_tetrahedra':live.size(),
            'moves':moves,'original_moves':len(data['moves']),'erased_moves':len(data['moves'])-len(moves),
            'loop_scan_max_tetrahedra':limit,'skipped_isomorphism_loops':skips,
            'all_moves_guarded':True,'native_simplify_used':False,
            'explicit_relabelling_moves_in_output':0,'seconds':time.monotonic()-now}
    _,result['replay']=replay(result);save(output,result)
    endpoint=Path(output).with_suffix('.endpoint.rga')
    require(regina.PacketOfTriangulation3(live).save(str(endpoint)),'endpoint save failed')
    print(json.dumps({'status':result['status'],'original_moves':len(data['moves']),'retained_moves':len(moves),
                      'end_isosig':live.isoSig(),'seconds':result['seconds']}),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('paths',nargs='+')
    p.add_argument('--out',required=True);p.add_argument('--limit',type=int,default=200)
    a=p.parse_args();shorten(a.paths,a.out,a.limit)

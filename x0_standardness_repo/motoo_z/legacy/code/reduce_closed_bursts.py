#!/usr/bin/env python3
"""Guarded closed-4D excursions with growth bursts and forced vertex births.

No native 4D simplify routine is used. Immutable committed log chunks and
numbered RGA checkpoints make an atomic CHECKPOINT.json sufficient to recover
the last completed chunk. Pending rows stay in bounded memory; no log writer
remains open between moves. A crash loses at most one uncommitted chunk.
The only standardness target is the explicit double of a four-simplex.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import signal
import time
import regina


def require(ok,message):
    if not ok:raise ValueError(message)


def sha_file(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def state_sha(tri):return hashlib.sha256(tri.tightEncoding().encode()).hexdigest()


def atomic_json(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def save_tri(path,tri):
    path=Path(path);tmp=path.with_suffix('.tmp.rga')
    require(regina.PacketOfTriangulation4(tri).save(str(tmp)),'RGA save failed')
    tmp.replace(path)


def load_tri(path):
    packet=regina.open(str(path));require(isinstance(packet,regina.Triangulation4),'not a4D triangulation')
    return regina.Triangulation4(packet)


def standard_double():
    tri=regina.Triangulation4();a,b=tri.newSimplices(2)
    for facet in range(5):a.join(facet,b,regina.Perm5())
    return tri


def standard_target(tri):
    return tri.size()==2 and tri.isIsomorphicTo(standard_double()) is not None


def admit(tri):
    require(tri.size()>0 and tri.isValid() and tri.isConnected() and
            tri.isOrientable() and not tri.isIdeal() and tri.isClosed(),
            'requires a finite closed connected orientable PL4-manifold')
    result={'pentachora':tri.size(),'vertices':tri.countVertices(),'closed':True,
            'valid':True,'orientable':True,'connected':True,
            'recognition':'Regina7.4.1 source manifold admission'}
    runtime=os.environ.get('MOTOO_RECOGNITION_RUNTIME')
    if runtime:
        result['recognition']='Regina7.4.1 kernel with separately audited positive primal-group 3D recognition shortcut'
        result['recognition_runtime_file']=runtime
    return result


def perform(tri,operation):
    require(isinstance(operation,list) and len(operation)==3,'malformed operation')
    kind,dim,index=operation
    require(type(dim) is int and dim in range(5),'bad dimension')
    require(type(index) is int and 0<=index<(tri.size() if dim==4 else tri.countFaces(dim)),'bad face index')
    face=tri.simplex(index) if dim==4 else tri.face(dim,index)
    n=tri.size();vertices=tri.countVertices()
    if kind=='P':
        require(tri.hasPachner(face),'Pachner guard rejected')
        require(tri.pachner(face) and tri.size()==n+2*dim-4,'checked Pachner failed')
        if dim==4:require(tri.countVertices()==vertices+1,'1-5 did not create one vertex')
        if dim==0:require(tri.countVertices()==vertices-1,'5-1 did not remove one vertex')
    elif kind=='Z':
        require(dim in (0,1,2) and tri.has20(face),'2-0 guard rejected')
        require(tri.move20(face) and tri.size()==n-2,'checked2-0 failed')
    elif kind=='C':
        require(dim==1 and face.vertex(0)!=face.vertex(1) and tri.hasCollapseEdge(face),'collapse guard rejected')
        require(tri.collapseEdge(face) and tri.size()<n and tri.countVertices()==vertices-1,'checked collapse failed')
    elif kind=='F':
        require(dim==1 and tri.has44(face),'4-4 guard rejected')
        require(tri.move44(face) and tri.size()==n,'checked4-4 failed')
    else:raise ValueError('unsupported move code')


def candidates(tri,kind,dim):
    guard={'P':tri.hasPachner,'Z':tri.has20,'C':tri.hasCollapseEdge,'F':tri.has44}[kind]
    answer=[]
    if kind=='C' and tri.countVertices()<2:return answer
    for face in (tri.simplices() if dim==4 else tri.faces(dim)):
        if kind=='P' and dim<4 and face.degree()!=5-dim:continue
        if kind=='Z' and face.degree()!=2:continue
        if kind=='C' and face.vertex(0)==face.vertex(1):continue
        if kind=='F' and face.degree()!=4:continue
        if guard(face):answer.append([kind,dim,face.index()])
    return answer


class ChunkJournal:
    def __init__(self,out,tri,chunk_size,started):
        self.out=Path(out);self.chunk_size=chunk_size;self.started=started
        (self.out/'moves').mkdir();(self.out/'states').mkdir()
        self.n=0;self.current=state_sha(tri);self.best=regina.Triangulation4(tri)
        self.best_n=0;self.best_hash=self.current;self.max_size=tri.size()
        self.counts=Counter();self.chunks=[];self.pending_rows=[];self.pending_first=0
        self.commit(tri,'RUNNING',0)

    def record(self,tri,operation,phase,cycle):
        before=self.current;before_size=tri.size()
        perform(tri,operation)
        self.current=state_sha(tri)
        if not self.pending_rows:self.pending_first=self.n
        row={'index':self.n,'operation':operation,'before_sha256':before,
             'after_sha256':self.current,'before_pentachora':before_size,
             'pentachora':tri.size(),'phase':phase,'cycle':cycle}
        self.pending_rows.append(json.dumps(row,separators=(',',':'))+'\n')
        self.n+=1;self.counts[str(operation[:2])]+=1
        self.max_size=max(self.max_size,tri.size())
        if tri.size()<self.best.size():
            self.best=regina.Triangulation4(tri);self.best_n=self.n;self.best_hash=self.current
        if len(self.pending_rows)>=self.chunk_size:self.commit(tri,'RUNNING',cycle)

    def commit(self,tri,status,cycle):
        if self.pending_rows:
            final=self.out/'moves'/('%06d.jsonl'%len(self.chunks))
            temporary=final.with_suffix('.tmp')
            with temporary.open('w') as stream:
                stream.write(''.join(self.pending_rows));stream.flush();os.fsync(stream.fileno())
            temporary.replace(final)
            self.chunks.append({'file':str(final.relative_to(self.out)),
                                'sha256':sha_file(final),'first_index':self.pending_first,
                                'moves':len(self.pending_rows)})
            self.pending_rows=[]
        current_file=self.out/'states'/('state_%09d.rga'%self.n)
        if not current_file.exists():save_tri(current_file,tri)
        best_file=self.out/'states'/('best_%09d.rga'%self.best_n)
        if not best_file.exists():save_tri(best_file,self.best)
        report={'schema':'motoo.closed-bursts.checkpoint.v1','status':status,
                'moves':self.n,'current_pentachora':tri.size(),'current_tight_sha256':self.current,
                'current_file':str(current_file.relative_to(self.out)),'current_rga_sha256':sha_file(current_file),
                'best_moves':self.best_n,'best_pentachora':self.best.size(),'best_tight_sha256':self.best_hash,
                'best_file':str(best_file.relative_to(self.out)),'best_rga_sha256':sha_file(best_file),
                'chunks':self.chunks,'counts':dict(self.counts),'max_pentachora':self.max_size,
                'cycle':cycle,'elapsed_seconds':time.monotonic()-self.started,
                'standard_endpoint':standard_target(tri),'replay_verified':False,
                'native_4d_simplify_used':False}
        atomic_json(self.out/'CHECKPOINT.json',report)
        print(json.dumps({k:report[k] for k in ['status','moves','current_pentachora','best_pentachora','cycle','elapsed_seconds','standard_endpoint']}),flush=True)
        return report


def search(source,out,seconds,seed,growth,burst24,mix,birth_every,birth_vertices,chunk_size,initial_burst=False,tri=None):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina7.4.1')
    require(seconds>0 and growth>=4 and burst24>0 and mix>=0 and birth_every>0 and birth_vertices>=2 and chunk_size>0,'invalid search budget/policy')
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    if tri is None:tri=load_tri(source)
    else:require(tri.tightEncoding()==load_tri(source).tightEncoding(),'preloaded source differs')
    print(json.dumps({'stage':'source_admission','pentachora':tri.size()}),flush=True)
    source_info=admit(tri);started=time.monotonic();rng=random.Random(seed)
    policy={'growth':growth,'burst24':burst24,'mix':mix,'birth_every':birth_every,
            'birth_target_vertices':birth_vertices,'chunk_size':chunk_size,'initial_burst':initial_burst,
            'birth_cycles_preserve_vertices_from_edge_collapse':True}
    atomic_json(out/'SOURCE.json',{'schema':'motoo.closed-bursts.source.v1','source_file':str(source.resolve()),
                'source_sha256':sha_file(source),'start_tight_sha256':state_sha(tri),'source_admission':source_info,
                'regina':'7.4.1','seed':seed,'policy':policy,'code_sha256':sha_file(__file__),
                'target':'explicit double of Delta4, all five gluings identity','native_4d_simplify_used':False})
    journal=ChunkJournal(out,tri,chunk_size,started);stopping=False
    def stop(signum,frame):
        nonlocal stopping
        stopping=True
    handlers={s:signal.signal(s,stop) for s in (signal.SIGINT,signal.SIGTERM)}
    cycle=0;reason='WALL_BUDGET'
    def live():return not stopping and time.monotonic()-started<seconds and not standard_target(tri)
    def choose(kind,dim,phase):
        options=candidates(tri,kind,dim)
        if not options:return False
        journal.record(tri,rng.choice(options),phase,cycle);return True
    def descent(allow_collapse):
        changed=False
        while live():
            found=False
            order=([('C',1)] if allow_collapse else [])+[('Z',2),('Z',1),('P',1),('P',0),('Z',0)]
            for kind,dim in order:
                if choose(kind,dim,'descent'):
                    found=True;changed=True;break
            if not found:break
        return changed
    try:
        if not initial_burst:descent(True)
        while live():
            cycle+=1;made_move=False
            ceiling=journal.best.size()+growth
            birth_cycle=(cycle-1)%birth_every==0
            if birth_cycle:
                births=max(1,birth_vertices-tri.countVertices())
                for _ in range(births):
                    if not live() or tri.size()+4>ceiling:break
                    if choose('P',4,'forced_vertex_birth'):made_move=True
            for _ in range(burst24):
                if not live() or tri.size()+2>ceiling:break
                if not choose('P',3,'growth_burst'):break
                made_move=True
            for _ in range(mix):
                if not live():break
                options=candidates(tri,'P',2)+candidates(tri,'F',1)
                if not options:break
                journal.record(tri,rng.choice(options),'mix',cycle);made_move=True
            made_move=descent(not birth_cycle) or made_move
            if not made_move:
                reason='NO_ALLOWED_MOVE_WITHIN_POLICY';break
        if standard_target(tri):reason='STANDARD_SPHERE_REACHED'
        elif stopping:reason='STOP_REQUESTED_AFTER_COMPLETE_MOVE'
        report=journal.commit(tri,reason,cycle)
    finally:
        for sig,handler in handlers.items():signal.signal(sig,handler)
    return report


def replay(directory,out,source_override=None):
    require(importlib.metadata.version('regina')=='7.4.1','requires Regina7.4.1')
    directory=Path(directory);meta=json.loads((directory/'SOURCE.json').read_text())
    check=json.loads((directory/'CHECKPOINT.json').read_text())
    source=Path(source_override) if source_override is not None else Path(meta['source_file'])
    require(sha_file(source)==meta['source_sha256'],'source checksum mismatch')
    tri=load_tri(source);initial=admit(tri);current=state_sha(tri)
    require(current==meta['start_tight_sha256'],'source state mismatch')
    started=time.monotonic();n=0;counts=Counter();best_seen=tri.size()
    for chunk in check['chunks']:
        path=directory/chunk['file']
        require(sha_file(path)==chunk['sha256'] and chunk['first_index']==n,'chunk checksum/order mismatch')
        local_count=0
        for line in path.read_text().splitlines():
            row=json.loads(line)
            require(row['index']==n and row['before_sha256']==current,'broken move chain')
            perform(tri,row['operation']);current=state_sha(tri)
            require(current==row['after_sha256'] and tri.size()==row['pentachora'],'replay state mismatch')
            n+=1;local_count+=1;counts[str(row['operation'][:2])]+=1;best_seen=min(best_seen,tri.size())
        require(local_count==chunk['moves'],'wrong chunk length')
    require(n==check['moves'] and current==check['current_tight_sha256'],'wrong replay endpoint')
    saved=directory/check['current_file']
    require(sha_file(saved)==check['current_rga_sha256'],'saved RGA checksum mismatch')
    require(load_tri(saved).tightEncoding()==tri.tightEncoding(),'saved RGA state differs')
    require(dict(counts)==check['counts'],'move counts differ')
    require(best_seen==check['best_pentachora'],'minimum size differs')
    final=admit(tri);standard=standard_target(tri)
    require(standard==check['standard_endpoint'],'claimed standard endpoint differs')
    result={'status':'PASS_STANDARDIZATION' if standard else 'PASS_PARTIAL_PATH','moves':n,'counts':dict(counts),
            'source_sha256':meta['source_sha256'],'source_admission':initial,'end_admission':final,
            'end_tight_sha256':current,'standard_endpoint':standard,'seconds':time.monotonic()-started,
            'trust':'Regina7.4.1 source admission and checked local moves; explicit target isomorphism',
            'native_4d_simplify_used':False}
    atomic_json(out,result);print(json.dumps(result),flush=True);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    s=sub.add_parser('search');s.add_argument('source',type=Path);s.add_argument('--out',type=Path,required=True)
    s.add_argument('--seconds',type=float,default=600);s.add_argument('--seed',type=int,default=0)
    s.add_argument('--growth',type=int,default=40);s.add_argument('--burst24',type=int,default=12)
    s.add_argument('--mix',type=int,default=80);s.add_argument('--birth-every',type=int,default=3)
    s.add_argument('--birth-vertices',type=int,default=5);s.add_argument('--chunk-size',type=int,default=200)
    s.add_argument('--initial-burst',action='store_true')
    r=sub.add_parser('replay');r.add_argument('directory',type=Path);r.add_argument('--out',type=Path,required=True)
    r.add_argument('--source',type=Path,help='relocated source file; certificate SHA-256 is still enforced')
    a=parser.parse_args()
    if a.command=='search':search(a.source,a.out,a.seconds,a.seed,a.growth,a.burst24,a.mix,a.birth_every,a.birth_vertices,a.chunk_size,a.initial_burst)
    else:replay(a.directory,a.out,a.source)


if __name__=='__main__':main()

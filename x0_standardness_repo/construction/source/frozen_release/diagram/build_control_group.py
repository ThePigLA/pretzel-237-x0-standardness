#!/usr/bin/env python3
"""Bounded exact edge/triangle presentation of the immutable closed control."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip, heapq, json, resource, sys, time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lib'))
from group_cells import presentation, cyc, inv

def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def filehash(p):return sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    b=Path(p).read_bytes();return json.loads(gzip.decompress(b) if str(p).endswith('.gz') else b)
def save(p,x):
    b=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(p).write_bytes(gzip.compress(b,mtime=0) if str(p).endswith('.gz') else b)

def bounded_eliminate(p,deadline,max_word=200000,max_total_letters=20000000):
    """One-occurrence Tietze moves; preflight each complete substitution step."""
    rel={};occ=defaultdict(set);heap=[];next_id=0;letters=0
    active=set(range(1,p['generators']+1));events=[];peak_word=0;peak_letters=0
    def put(w):
        nonlocal next_id,letters,peak_word,peak_letters
        if not w:return
        next_id+=1;rel[next_id]=w;letters+=len(w)
        peak_word=max(peak_word,len(w));peak_letters=max(peak_letters,letters)
        for a in set(map(abs,w)):occ[a].add(next_id)
        heapq.heappush(heap,(len(w),next_id))
    def erase(i):
        nonlocal letters
        w=rel.pop(i);letters-=len(w)
        for a in set(map(abs,w)):occ[a].remove(i)
    for w in p['relators']:put(tuple(w))
    stop='one_occurrence_queue_exhausted';notice=time.monotonic()
    while heap:
        if time.monotonic()>=deadline:stop='producer_time_cap';break
        _,i=heapq.heappop(heap)
        if i not in rel:continue
        w=rel[i];counts=Counter(map(abs,w));single=[a for a in counts if counts[a]==1]
        if not single:continue
        g=min(single,key=lambda a:(len(occ[a]),a));j=next(k for k,a in enumerate(w) if abs(a)==g)
        tail=w[j+1:]+w[:j];replacement=inv(tail) if w[j]>0 else tail
        inverse=inv(replacement);affected=sorted(occ[g]);prepared={};bad=None
        remaining_letters=letters-sum(len(rel[z]) for z in affected)
        for z in affected:
            if z==i:continue
            expanded=[]
            for a in rel[z]:
                expanded.extend(replacement if a==g else inverse if a==-g else (a,))
                if len(expanded)>max_word*4:bad='word_expansion_cap';break
            if bad:break
            reduced=cyc(expanded)
            if len(reduced)>max_word:bad='word_cap';break
            prepared[z]=reduced;remaining_letters+=len(reduced)
            if remaining_letters>max_total_letters:bad='total_word_cap';break
            if time.monotonic()>=deadline:bad='producer_time_cap';break
        if bad:stop=bad;break
        # All substitutions have been checked before this atomic Tietze step.
        events.append({'generator':g,'relator':list(w),'replacement':list(replacement)})
        for z in affected:erase(z)
        for z in affected:
            if z!=i:put(prepared[z])
        active.remove(g)
        if time.monotonic()-notice>=10:
            print(f'{len(events)} exact eliminations; {len(active)} generators; {len(rel)} relators',flush=True)
            notice=time.monotonic()
    return {'eliminations':events,'generators':sorted(active),
            'relators':[list(w) for w in sorted(set(rel.values()))],
            'stop_reason':stop,'peak_word_length':peak_word,'peak_live_letters':peak_letters}

def main():
    started=time.monotonic();cpu=time.process_time()
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    folder=ROOT/'group_audit';folder.mkdir(exist_ok=False)
    master_path=ROOT/'control/marked_control.json.gz';master=load(master_path)
    if digest(master['facets'])!=master['facets_sha256']:raise ValueError('Master facet digest')
    vertices=sorted({v for f in master['facets'] for v in f})
    edges=sorted({e for f in master['facets'] for e in combinations(f,2)})
    triangles=sorted({t for f in master['facets'] for t in combinations(f,3)})
    ei={e:i+1 for i,e in enumerate(edges)}
    triangle_boundaries=[[ei[a,b],ei[b,c],-ei[a,c]] for a,b,c in triangles]
    sk={'vertices':vertices,'edge_endpoints':edges,'triangle_boundaries':triangle_boundaries}
    p,_,_=presentation(sk)
    print(f'Full source skeleton: {len(vertices)} vertices, {len(edges)} edges, {len(triangles)} triangles; '
          f'{p["generators"]} generators, {len(p["relators"])} nonempty distinct relators',flush=True)
    skeleton={'schema':'full-immutable-control-two-skeleton-v1','vertices':vertices,'edges':edges,'triangles':triangles,
              'source_facets_sha256':master['facets_sha256']}
    binding={'schema':'immutable-control-group-source-binding-v1','source_file':'control/marked_control.json.gz',
             'source_file_sha256':filehash(master_path),'source_facets_sha256':master['facets_sha256'],
             'native_admission_file_sha256':filehash(ROOT/'control/cover_native_admission.json'),
             'native_file_sha256':filehash(ROOT/'control/cover_regina.json.gz'),
             'skeleton_sha256':digest(skeleton),'initial_presentation_sha256':digest(p)}
    save(folder/'source_binding.json',binding)
    save(folder/'skeleton.json.gz',skeleton);save(folder/'initial_presentation.json.gz',p)
    reduced=bounded_eliminate(p,started+80)
    certificate={'schema':'exact-one-occurrence-Tietze-control-group-v1',
                 'source_binding_sha256':digest(binding),'initial_presentation_sha256':digest(p),
                 'permitted_operations':['free word reduction','cyclic relator reduction',
                     'duplicate/empty relator deletion','one-occurrence generator and defining-relator elimination'],
                 **reduced}
    save(folder/'tietze_certificate.json.gz',certificate)
    result={'schema':'immutable-control-group-producer-v1','source_facets_sha256':master['facets_sha256'],
            'source_skeleton_f_vector':[len(vertices),len(edges),len(triangles)],
            'initial_generators':p['generators'],'initial_relators':len(p['relators']),
            'eliminations':len(reduced['eliminations']),'remaining_generators':len(reduced['generators']),
            'remaining_relators':len(reduced['relators']),'stop_reason':reduced['stop_reason'],
            'peak_word_length':reduced['peak_word_length'],'peak_live_letters':reduced['peak_live_letters'],
            'wall_seconds':time.monotonic()-started,'cpu_seconds':time.process_time()-cpu,
            'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            'budget':{'producer_seconds':90,'elimination_deadline_from_start_seconds':80,'address_space_bytes':8*1024**3},
            'scope':'Producer checkpoint; independent replay required for the group conclusion.'}
    save(folder/'producer_summary.json',result);print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()

"""Source-bound fresh replay of reviewer X0 elementary certificate, steps 0--2.

The mathematical checkers are the supplied, unchanged reviewer implementations.
This analyst wrapper verifies all input and move schemas they assume, compares
the original release bytes, reconstructs complete face-pairing exports, executes
the actual checkers in bounded unoptimized Python subprocesses, and binds the
fresh endpoint to the exact native source. No Regina is imported or required.
"""
from pathlib import Path
from collections import Counter, defaultdict
from itertools import combinations
import argparse
import gzip
import hashlib
import json
import os
import resource
import signal
import subprocess
import sys
import time


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def validate_facets(facets):
    require(type(facets) is list and facets, 'empty or non-list facets')
    for f in facets:
        require(type(f) is list and len(f) == 5, 'facet arity')
        require(all(type(v) is int and v >= 0 for v in f), 'facet vertex type')
        require(len(set(f)) == 5, 'folded facet')
    require(len({tuple(sorted(f)) for f in facets}) == len(facets), 'duplicate abstract facet')


def closed_export(facets, data):
    """Require each abstract ridge's unique incidence pair and actual label transport."""
    n = len(facets)
    require(data['simplex_count'] == n, 'export simplex count')
    incidence = defaultdict(list)
    for i, f in enumerate(facets):
        for k in range(5):
            incidence[tuple(sorted(f[:k] + f[k+1:]))].append((i, k))
    require(all(len(o) == 2 for o in incidence.values()), 'nonclosed or branched ridge')
    expected = {frozenset(o) for o in incidence.values()}
    seen, halves, adj = set(), set(), defaultdict(list)
    for row in data['gluings']:
        require(type(row) is list and len(row) == 4, 'gluing row schema')
        i, k, j, p = row
        require(all(type(x) is int for x in (i,k,j)), 'gluing index type')
        require(0 <= i < n and 0 <= j < n and 0 <= k < 5 and i != j, 'gluing range')
        require(type(p) is list and len(p) == 5 and all(type(x) is int for x in p)
                and sorted(p) == list(range(5)), 'gluing permutation')
        pair = frozenset(((i,k),(j,p[k])))
        require(pair in expected and pair not in seen and not (pair & halves), 'incorrect/repeated face pairing')
        require(all(facets[i][x] == facets[j][p[x]] for x in range(5) if x != k), 'incorrect vertex transport')
        seen.add(pair); halves.update(pair)
        parity = -1 if sum(p[x]>p[y] for x in range(5) for y in range(x+1,5)) % 2 else 1
        adj[i].append((j,-parity)); adj[j].append((i,-parity))
    require(seen == expected and len(halves) == 5*n, 'incomplete export')
    orient = {0:1}; todo = [0]
    while todo:
        i = todo.pop()
        for j, transport in adj[i]:
            wanted = orient[i]*transport
            if j in orient:
                require(orient[j] == wanted, 'orientation conflict')
            else:
                orient[j] = wanted; todo.append(j)
    require(len(orient) == n, 'disconnected facet dual graph')
    return {'simplex_count':n,'tetrahedron_pairings':len(seen),'half_facets_checked':len(halves),
            'all_half_facets_used_once':True,'every_gluing_is_exact_abstract_face_identification':True,
            'connected':True,'orientable':True,'closed':True}


def run_checked(label, script, args, output):
    logfile = output/(label+'.log'); metrics = output/(label+'_resource.json')
    wrapper = r'''
import json,resource,runpy,sys,time
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
resource.setrlimit(resource.RLIMIT_CPU,(900,901))
out,script,*args=sys.argv[1:]
if sys.flags.optimize:raise RuntimeError('assertions disabled')
sys.argv=[script,*args]
start=time.monotonic()
try:
    runpy.run_path(script,run_name='__main__')
finally:
    r=resource.getrusage(resource.RUSAGE_SELF)
    json.dump({'wall_seconds':time.monotonic()-start,'cpu_seconds':r.ru_utime+r.ru_stime,
      'peak_rss_bytes':r.ru_maxrss*(1 if sys.platform=='darwin' else 1024),
      'rlimit_as_bytes':8*1024**3,'rlimit_cpu_seconds':900,'assertions_enabled':not sys.flags.optimize},open(out,'w'),indent=2)
'''
    env = dict(os.environ)
    env.pop('PYTHONOPTIMIZE',None); env.pop('PYTHONPATH',None)
    env['PYTHONDONTWRITEBYTECODE']='1'
    start=time.monotonic()
    with logfile.open('w') as stream:
        child=subprocess.Popen([sys.executable,'-E','-B','-u','-c',wrapper,str(metrics),str(script),*map(str,args)],
            stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,env=env)
        try:
            while True:
                try:
                    code=child.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    elapsed=time.monotonic()-start
                    print(f'{label}: running, {elapsed:.0f}s elapsed',flush=True)
                    if elapsed>1200:raise TimeoutError(label+' wall limit')
        except BaseException:
            os.killpg(child.pid,signal.SIGKILL);child.wait();raise
    require(code==0,label+' checker failed; see '+str(logfile))
    r=json.loads(metrics.read_text())
    r.update({'returncode':code,'wall_ceiling_seconds':1200,'script_sha256':sha(script),
              'log_file':logfile.name,'log_sha256':sha(logfile),'resource_file':metrics.name,
              'resource_sha256':sha(metrics),'parent_wall_seconds':time.monotonic()-start})
    print(label+': '+logfile.read_text().strip(),flush=True)
    return r


def main():
    base=Path(__file__).resolve().parents[2]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reviewer-root',type=Path,default=base/'reviewer_x0_census/x0_census_reduction')
    p.add_argument('--original-root',type=Path,default=base/'x0_preparation_v1')
    p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'fresh_replay')
    p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args(); review=a.reviewer_root.resolve(); original=a.original_root.resolve(); out=a.output_dir.resolve()
    require(not sys.flags.optimize,'Run without -O; reviewer checks use assertions')
    require(not out.exists(),'Output directory must be new')
    require(review not in out.parents and original not in out.parents,'Output cannot be in frozen inputs')
    out.mkdir(parents=True)
    start=time.monotonic(); cpu=time.process_time()
    manifest_path=original/'MANIFEST_SHA256.json'
    require(sha(manifest_path)=='77a3af20d9d26308d113b4cb6cc48961c874e8675227516d01ffa4b4836e5d0b','original release manifest changed')
    manifest=json.loads(manifest_path.read_text())['files']
    pairs={
      'inputs/marked_control.json.gz':'candidate/marked_candidate.json.gz',
      'inputs/cover_regina.json.gz':'candidate/cover_regina.json.gz',
      'inputs/unmarked_contractions_v1/source_facets.json.gz':'unmarked_contractions_v1/source_facets.json.gz',
      'inputs/unmarked_contractions_v1/certificate.json.gz':'unmarked_contractions_v1/certificate.json.gz',
      'inputs/unmarked_contractions_v1/endpoint.json.gz':'unmarked_contractions_v1/endpoint.json.gz',
      'inputs/unmarked_contractions_v1/endpoint_regina.json.gz':'unmarked_contractions_v1/endpoint_regina.json.gz',
      'step0_links/x0_vertex_links.json.gz':'portable_links/x0_vertex_links.json.gz',
    }
    hashes={}
    for name,old in pairs.items():
        h=sha(review/name);require(h==sha(original/old)==manifest[old]['sha256'],'original input mismatch '+name)
        require((original/old).stat().st_size==manifest[old]['bytes'],'original manifest size '+old)
        hashes[name]=h
    all_review_hashes={str(f.relative_to(review)):sha(f) for f in review.rglob('*') if f.is_file()}
    master=read(review/'inputs/marked_control.json.gz'); endpoint=read(review/'inputs/unmarked_contractions_v1/endpoint.json.gz')
    source=read(review/'inputs/unmarked_contractions_v1/source_facets.json.gz')
    contraction=read(review/'inputs/unmarked_contractions_v1/certificate.json.gz')
    links=read(review/'step0_links/x0_vertex_links.json.gz')
    for d in (master,source,endpoint):validate_facets(d['facets'])
    require(len(master['facets'])==103744 and len(endpoint['facets'])==44620,'endpoint sizes')
    require(source['facets']==master['facets'],'trace source facet order differs')
    require(canonical_digest(master['facets'])==master['facets_sha256']=='f40b43c09373dff94cde542a3076a8c2f482b1632726b63eed2c5610ddd1e92b','master facet hash')
    require(canonical_digest(endpoint['facets'])==endpoint['facets_sha256'],'endpoint facet hash')
    vertices={v for f in master['facets'] for v in f}
    require(vertices==set(range(4310)) and len(master['vertices'])==4310,'master vertices')
    require(set(links['links'])=={str(v) for v in vertices},'link key coverage')
    require(links['complete'] is True and links['vertices']==4310,'incomplete links certificate')
    require(links['input_bindings']['x0_master_file_sha256']==hashes['inputs/marked_control.json.gz'],'links source binding')
    require(sha(review/'step0_links/check_links.py')==links['input_bindings']['checker_file_sha256'],'links checker binding')
    link_arities=Counter()
    for vertex,moves in links['links'].items():
        require(type(moves) is list,'link log schema')
        for m in moves:
            require(type(m) is list and len(m)==2,'link move schema')
            s,t=m
            require(type(s) is list and type(t) is list and s and t,'empty link flip face')
            require(all(type(v) is int and v>=0 for v in s+t),'link vertex schema')
            require(len(set(s))==len(s) and len(set(t))==len(t),'repeated link vertex')
            require(not set(s)&set(t) and len(s)+len(t)==5,'invalid link flip arity')
            link_arities[f'{len(t)}-{len(s)}']+=1
    require(sum(link_arities.values())==links['flips_total']==529460,'link flip total')
    require(len(contraction['steps'])==3308,'wrong prefix contraction count')
    for s in contraction['steps']:
        require(all(type(s[k]) is int for k in ('u','v','before','after')),'contraction schema')
        require(s['u']>=0 and s['v']>=0 and s['u']!=s['v'],'invalid contraction endpoints')
    exports={'master':closed_export(master['facets'],read(review/'inputs/cover_regina.json.gz')),
             'prefix_endpoint':closed_export(endpoint['facets'],read(review/'inputs/unmarked_contractions_v1/endpoint_regina.json.gz'))}
    # Verify the reviewer-described cap-only difference using the original frozen X2 master.
    x2path=original/'frozen_release/control/marked_control.json.gz'
    require(sha(x2path)==manifest['frozen_release/control/marked_control.json.gz']['sha256'],'X2 master changed')
    x2=read(x2path); f0={tuple(f) for f in master['facets']};f2={tuple(f) for f in x2['facets']}
    require(len(f0&f2)==100736 and len(f0-f2)==len(f2-f0)==3008,'cap difference count')
    for key in ('vertices','deck','ramification_triangles','twin_triangles'):
        require(master[key]==x2[key],'marked '+key+' changed')
    for key in ('N','E','central','bridge'):
        require(master['parts'][key]==x2['parts'][key],'noncap '+key+' changed')
    require(f0-f2<={tuple(f) for f in master['parts']['cap']},'new noncap facet')
    require(f2-f0<={tuple(f) for f in x2['parts']['cap']},'removed noncap facet')
    cap_comparison={'shared_facets':100736,'removed_X2_cap_facets':3008,'new_X0_cap_facets':3008,
      'vertices_deck_branch_twin_identical':True,'all_noncap_parts_identical':True,'x2_master_sha256':sha(x2path)}
    # Stream the complete move vocabulary; do not silently accept unknown tags or empty faces.
    move_path=review/'step2_flips/moves_44620_to_1332.jsonl.gz'
    counts=Counter();types=Counter();lines=0
    with gzip.open(move_path,'rt') as stream:
        for n,line in enumerate(stream,1):
            require(bool(line.strip()),f'blank move line {n}')
            m=json.loads(line)
            require(type(m) is list and len(m)==3 and m[0] in ('c','f'),f'move vocabulary line {n}')
            if m[0]=='c':
                require(type(m[1]) is int and type(m[2]) is int and m[1]>=0 and m[2]>=0 and m[1]!=m[2],f'self/invalid contraction line {n}')
            else:
                s,t=m[1:]
                require(type(s) is list and type(t) is list and s and t,f'empty flip line {n}')
                require(all(type(v) is int and v>=0 for v in s+t),f'flip label line {n}')
                require(len(set(s))==len(s) and len(set(t))==len(t) and not set(s)&set(t) and len(s)+len(t)==6,f'flip schema line {n}')
                types[f'{len(t)}-{len(s)}']+=1
            counts[m[0]]+=1;lines+=1
    require(lines==629874,'move count differs')
    result={'status':'PREFLIGHT_PASS','schema':'X0-census-elementary-audit-v1','python':sys.version,
      'original_manifest_sha256':sha(manifest_path),'input_file_sha256':hashes,'reviewer_file_sha256':all_review_hashes,
      'audit_script_sha256':sha(Path(__file__)),'master_facets_sha256':master['facets_sha256'],
      'prefix_endpoint_facets_sha256':endpoint['facets_sha256'],'link_moves':sum(link_arities.values()),
      'link_bistellar_types':dict(sorted(link_arities.items())),'prefix_contractions':3308,
      'step2_moves':lines,'step2_contractions':counts['c'],'step2_bistellar_flips':counts['f'],
      'step2_bistellar_types':dict(sorted(types.items())),'closed_export_checks':exports,
      'cap_only_comparison':cap_comparison,'all_move_schemas_checked':True,'all_contraction_endpoints_distinct':True,
      'every_vertex_link_key_present_exactly_once':True,'native_engine_used':False,
      'preflight_cpu_seconds':time.process_time()-cpu,'preflight_wall_seconds':time.monotonic()-start}
    (out/'preflight.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PREFLIGHT PASS: seven original copied inputs, complete exports, cap-only comparison, all move schemas',flush=True)
    if a.preflight_only:return
    result['step0']=run_checked('step0_links',review/'step0_links/check_links.py',
        ['verify',review/'inputs/marked_control.json.gz',review/'step0_links/x0_vertex_links.json.gz'],out)
    result['step1']=run_checked('step1_contractions',review/'step1_contractions/check_contractions.py',[review/'inputs'],out)
    replayed=out/'replayed_1332.json.gz'
    result['step2']=run_checked('step2_moves',review/'step2_flips/replay_moves.py',
        [review/'inputs/unmarked_contractions_v1/endpoint.json.gz',move_path,replayed],out)
    actual=read(replayed); target=read(review/'step3_macro/x0_C1332_regina.json.gz')
    validate_facets(actual['facets'])
    require((actual['simplex_count'],actual['facets'],actual['gluings'])==
            (target['simplex_count'],target['facets'],target['gluings']),'actual replay differs from macro source')
    require(actual['simplex_count']==1332,'wrong final facet count')
    require(len({v for f in actual['facets'] for v in f})==45,'wrong final vertex count')
    endpoint_export=closed_export(actual['facets'],actual)
    fv=[len({s for f in actual['facets'] for s in combinations(f,k)}) for k in range(1,6)]
    require(sum((-1)**i*n for i,n in enumerate(fv))==2,'Euler characteristic not 2')
    require(sha(review/'inputs/marked_control.json.gz')==hashes['inputs/marked_control.json.gz'],'source changed')
    for name,h in all_review_hashes.items():require(sha(review/name)==h,'reviewer file changed during replay: '+name)
    result.update({'status':'PASS','source_files_unchanged':True,'endpoint_facets':1332,'endpoint_vertices':45,
      'endpoint_f_vector':fv,'endpoint_export_check':endpoint_export,
      'replay_equals_macro_source_exactly_including_ordered_facets_and_gluings':True,
      'endpoint_facets_sha256_json_default':hashlib.sha256(json.dumps(actual['facets']).encode()).hexdigest(),
      'endpoint_facets_sha256_json_canonical':canonical_digest(actual['facets']),
      'replayed_export_sha256':sha(replayed),'macro_source_sha256':sha(review/'step3_macro/x0_C1332_regina.json.gz'),
      'total_wall_seconds':time.monotonic()-start,'wrapper_cpu_seconds':time.process_time()-cpu,
      'proof_scope':'The exact X0 marked master is a closed combinatorial four-manifold; the 3308 full-link prefix contractions and 629874 subsequent full-link contractions/bistellar flips preserve PL type and reach exactly the closed 1332-facet native source. No claim of PL standardness is made.'})
    (out/'audit_steps012.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','step2_moves','step2_contractions','step2_bistellar_types','endpoint_f_vector','total_wall_seconds')},indent=2),flush=True)


if __name__=='__main__':
    signal.signal(signal.SIGTERM,lambda *_:(_ for _ in ()).throw(KeyboardInterrupt()))
    main()

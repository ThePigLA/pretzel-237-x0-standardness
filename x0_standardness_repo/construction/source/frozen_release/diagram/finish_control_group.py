#!/usr/bin/env python3
"""Exact logged relator shortening on the small residual master presentation."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import gzip,json,resource,sys,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lib'))
from group_cells import cyc,inv
def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()))
def save(p,x):
    b=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
    p.write_bytes(gzip.compress(b,mtime=0) if str(p).endswith('.gz') else b)

def shortening(rels,deadline):
    # Every rule is U=V^-1 from a retained defining relator UV.  Strictly more
    # than half the defining relation must match, so total length decreases.
    ordered=sorted(rels,key=lambda w:(len(w),w))
    for bound in (8,max(map(len,ordered),default=0)):
        rules={}
        for r in ordered:
            if len(r)>bound:continue
            for sign in (1,-1):
                q=r if sign==1 else inv(r)
                for j in range(len(q)):
                    rotated=q[j:]+q[:j]
                    for n in range(len(q)//2+1,len(q)+1):
                        rules.setdefault((n,rotated[:n]),[]).append((r,sign,j,inv(rotated[n:])))
            if time.monotonic()>=deadline:return None,'time_cap'
        sizes=sorted({k[0] for k in rules},reverse=True)
        for w in ordered:
            doubled=w+w
            for n in sizes:
                if n>len(w):continue
                for p in range(len(w)):
                    options=rules.get((n,doubled[p:p+n]),())
                    for r,sign,j,replacement in options:
                        if r==w:continue # Retain the rule while changing a different relator.
                        rotated=w[p:]+w[:p]
                        result=cyc(replacement+rotated[n:])
                        if len(result)>=len(w):raise ValueError('Nonshortening rule')
                        return {'kind':'relator_shortening','rule_relator':list(r),'rule_sign':sign,
                                'rule_start':j,'target':list(w),'match_start':p,'match_length':n,
                                'result':list(result)},None
            if time.monotonic()>=deadline:return None,'time_cap'
        if bound==max(map(len,ordered),default=0):break
    return None,'no_shortening_rule'

def main():
    started=time.monotonic();cpu=time.process_time()
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    folder=ROOT/'group_audit';prefix_path=folder/'tietze_certificate.json.gz';prefix=load(prefix_path)
    active=set(prefix['generators']);rels=set(map(tuple,prefix['relators']));events=[];stop='empty_presentation'
    deadline=started+25
    while active:
        if time.monotonic()>=deadline:stop='time_cap';break
        selected=None
        for w in sorted(rels,key=lambda x:(len(x),x)):
            counts=Counter(map(abs,w));single=[g for g,n in counts.items() if n==1]
            if single:
                g=min(single);j=next(i for i,a in enumerate(w) if abs(a)==g)
                tail=w[j+1:]+w[:j];replacement=inv(tail) if w[j]>0 else tail
                selected={'kind':'eliminate','generator':g,'relator':list(w),'replacement':list(replacement)};break
        if selected is not None:
            g=selected['generator'];w=tuple(selected['relator']);s=tuple(selected['replacement']);si=inv(s)
            rels.remove(w)
            rels={cyc(a for z in r for a in (s if z==g else si if z==-g else (z,))) for r in rels}
            rels.discard(());active.remove(g);events.append(selected);continue
        event,reason=shortening(rels,deadline)
        if event is None:stop=reason;break
        rels.remove(tuple(event['target']))
        if event['result']:rels.add(tuple(event['result']))
        events.append(event)
    result={'schema':'exact-residual-Tietze-control-group-v1','prefix_file_sha256':sha256(prefix_path.read_bytes()).hexdigest(),
            'initial_generators':prefix['generators'],'initial_relators_sha256':digest(prefix['relators']),
            'events':events,'generators':sorted(active),'relators':[list(w) for w in sorted(rels)],'stop_reason':stop,
            'scope':'Only retained-relator consequences and exact one-occurrence eliminations; no added quotient relation.',
            'wall_seconds':time.monotonic()-started,'cpu_seconds':time.process_time()-cpu,
            'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    save(folder/'residual_tietze_certificate.json.gz',result)
    summary={k:v for k,v in result.items() if k not in ('events','relators')}
    summary['events']=len(events);summary['remaining_relators']=len(rels)
    save(folder/'residual_producer_summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()

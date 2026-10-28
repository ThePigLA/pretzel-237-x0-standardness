#!/usr/bin/env python3
"""Combine every native worker ledger before allocating further computation."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json

ROOT=Path(__file__).resolve().parent

def main(reserved_replay):
    names=['native_budget_ledger.json',*[f'native_seed{i}_budget.json' for i in (1,2,3)],'native_final_replay_budget.json']
    runs=[];hashes={};active=[]
    for name in names:
        path=ROOT/'control'/name
        if not path.exists():continue
        hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        for i,r in enumerate(json.loads(path.read_text())['runs']):
            runs.append(dict(r,ledger=name,index=i))
            if r['status']=='RUNNING':active.append({'ledger':name,'index':i})
    if active:raise RuntimeError('Wait for active ledger runs before a new allocation: '+repr(active))
    diagnostic=45.0
    completed=sum(r.get('cpu_seconds',0.0) for r in runs)+diagnostic
    seeds=sorted({r['seed'] for r in runs if r.get('mode')=='search'})
    result={'schema':'aggregate-native-budget-v2','updated_utc':datetime.now(timezone.utc).isoformat(),
        'cpu_cap_seconds':14400,'max_distinct_seeds':16,'distinct_seeds':seeds,
        'completed_cpu_seconds':completed,'diagnostic_cpu_charge_seconds':diagnostic,
        'diagnostic_charge_basis':'Conservative full 45-second bound for the tight-decoding skeleton diagnostic; no result is claimed from it.',
        'reserved_final_replay_cpu_seconds':reserved_replay,
        'total_used_and_reserved_seconds':completed+reserved_replay,
        'within_cap':completed+reserved_replay<=14400 and len(seeds)<=16,
        'ledger_file_hashes':hashes,'runs':runs,
        'scope':'Native control admission, attempted native invariants, certified closed PL prefix, standardization searches, diagnostics, and final native replay. Includes superseded and interrupted runs. Source construction and other exact combinatorial audits are separate.',
        'search_scope':'Four-seed bounded batch. The 16-seed and four-CPU-hour values are ceilings, not claims that either ceiling was exhausted.'}
    if not result['within_cap']:raise RuntimeError('Aggregate reservation exceeds original cap')
    (ROOT/'control/native_budget_aggregate.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('completed_cpu_seconds','reserved_final_replay_cpu_seconds','total_used_and_reserved_seconds','within_cap','distinct_seeds')}))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--reserve-replay',type=float,default=0);args=a.parse_args()
    if args.reserve_replay<0:raise ValueError('negative reservation')
    main(args.reserve_replay)

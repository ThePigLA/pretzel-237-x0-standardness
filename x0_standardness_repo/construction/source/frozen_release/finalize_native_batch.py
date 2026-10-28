#!/usr/bin/env python3
"""Summarize a finished bounded batch without promoting partial replay to S4."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json

ROOT=Path(__file__).resolve().parent
def read(name):return json.loads((ROOT/name).read_text())
def sha(name):return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

def main():
    budget=read('control/native_budget_aggregate.json')
    milestone=read('control/MILESTONE_STATUS.json')
    if budget['reserved_final_replay_cpu_seconds']!=0 or not budget['within_cap']:raise ValueError('budget is not finalized')
    if any(r['status']=='RUNNING' for r in budget['runs']):raise ValueError('native workers still running')
    if not milestone['standardization']['search_finished_for_this_batch']:raise ValueError('batch still open')
    rows=[]
    for seed in budget['distinct_seeds']:
        cp=f'control/standardization_seed{seed}.json';c=read(cp)
        rp=f'control/standardization_seed{seed}_replay.json';r=read(rp) if (ROOT/rp).exists() else None
        runs=[x for x in budget['runs'] if x.get('mode')=='search' and x.get('seed')==seed]
        if r and (r['certificate_sha256']!=sha(cp) or not r['replay_passed']):raise ValueError('stale or failed final replay')
        rows.append({'seed':seed,'source_puncture_index':c['puncture'],'source_puncture_vertices':c['source_puncture_vertices'],
            'initial_local_pentachora':c['blocks'][0]['pentachora'],'saved_endpoint_pentachora':c['final_size'],
            'completed_macro_blocks':len(c['blocks']),'guarded_kicks':sum(len(b['kicks']) for b in c['blocks']),
            'last_search_worker_status':runs[-1]['status'],'charged_search_cpu_seconds':sum(x['cpu_seconds'] for x in runs),
            'certificate':cp,'certificate_sha256':sha(cp),'final_replay':rp if r else None,
            'final_replay_sha256':sha(rp) if r else None,'complete_path_replayed':bool(r),'standardized':bool(r and r['standardization_passed'])})
    best=min(rows,key=lambda x:x['saved_endpoint_pentachora'])
    if not best['complete_path_replayed']:raise ValueError('best endpoint has not received complete replay')
    replay=read(best['final_replay']);cert=read(best['certificate'])
    result={'schema':'finished-bounded-native-batch-v1','updated_utc':datetime.now(timezone.utc).isoformat(),
        'status':milestone['status'],'standardized':milestone['standardization']['passed'],
        'seeds':rows,'best_replayed_seed':best['seed'],'best_replayed_endpoint_pentachora':replay['endpoint_pentachora'],
        'best_replayed_endpoint_numbered_sha256':replay['endpoint_numbered_sha256'],
        'source_master_file_sha256':milestone['marked_master']['file_sha256'],
        'milestone_sha256':sha('control/MILESTONE_STATUS.json'),'aggregate_budget_sha256':sha('control/native_budget_aggregate.json'),
        'charged_cpu_seconds':budget['completed_cpu_seconds'],'cpu_ceiling_seconds':budget['cpu_cap_seconds'],
        'scope':'A complete replay to a non-reference endpoint certifies the partial PL simplification only. Interrupted macro objectives are not saved endpoints. The declared four-seed batch is finished; the original maxima of sixteen seeds and four CPU hours were not exhausted.',
        'closed_X0_constructed':False,'exoticness_inference':None}
    (ROOT/'control/NATIVE_BATCH_RESULT.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    table='\n'.join(f"| {x['seed']} | {x['initial_local_pentachora']:,} | {x['saved_endpoint_pentachora']:,} | {x['completed_macro_blocks']} | {'Passed' if x['complete_path_replayed'] else 'No separate final replay'} |" for x in rows)
    replay_charge=next(x['cpu_seconds'] for x in reversed(budget['runs']) if x.get('mode')=='replay')
    text=f'''# Finished native batch: exact outcome and replay

**Result:** `{milestone['status']}`.
The best complete replay ends at **{replay['endpoint_pentachora']:,} pentachora**.
Its `replay_passed` value is true and its `standardization_passed` value is
{str(replay['standardization_passed']).lower()}. The original marked control,
geometric identification, direct trivial-group certificate, and integral
homology certificate are unchanged.

## Accepted source-to-endpoint path

The path begins with the immutable 103,744-pentachoron marked master.
The independently replayed 3,307 full-link contractions produce a closed,
unmarked 44,728-pentachoron abstract simplicial complex of the same oriented
PL type. All markings remain on the master.

Seed {best['seed']} removes the interior of its embedded pentachoron at
zero-based index {cert['puncture']}, with vertices
`{cert['source_puncture_vertices']}`. The punctured input has 44,727 pentachora
and one standard three-sphere boundary. This is a puncture of the recorded
source, not a substitute model.

| Step | Pentachora after the step | Verification |
| --- | ---: | --- |
| Immutable closed master | 103,744 | Marked assembly and native PL-manifold admission |
| 3,307 closed full-link contractions | 44,728 | Independent standard-library replay |
| Recorded embedded-simplex puncture | 44,727 | Source facet labels and native gluing index |
| Initial local native simplification | {cert['blocks'][0]['pentachora']:,} | Exact numbered hash reproduced |
| {sum(len(b['kicks']) for b in cert['blocks'])} guarded 2-4 moves | {cert['blocks'][0]['pentachora']+2*sum(len(b['kicks']) for b in cert['blocks']):,} | Every native legality check replayed |
| Full native simplification | {replay['endpoint_pentachora']:,} | Exact numbered hash and final tight encoding reproduced |

The last numbered gluing hash is
`{replay['endpoint_numbered_sha256']}`.
The final replay checks the exact source file, master facets, contraction
prefix, puncture labels, Regina version, each move, each per-block macro
policy, each endpoint hash/size, and final tight encoding. Its receipt is
`control/standardization_seed{best['seed']}_replay.json`.

Native simplification is a version-pinned software trust boundary. These
receipts are macro replays, not an exported elementary Pachner transcript.
The final native replay took {replay['wall_seconds']:.3f} wall seconds and
was charged {replay_charge:.3f} CPU seconds by the external supervisor.

## Four-seed results

| Seed | Replayed initial local prefix | Last saved complete endpoint | Saved macro blocks | Separate complete-path replay |
| ---: | ---: | ---: | ---: | --- |
{table}

Seeds 0 and 3 finished their current worker normally. Seed 1 saved its
2,072-cell endpoint before a later trial reached its CPU cap. That later
candidate was discarded. Seed 2's full macro reached a reported objective
of 2,395, but it was interrupted before returning; the accepted checkpoint
remains 9,708. In-progress objective counts are not promoted to certificates.
Every initial local prefix was replayed during restoration. Only the final
receipt listed above certifies the entire best saved path.

## Budget and operational lessons

The charged native-control total is **{budget['completed_cpu_seconds']:.3f}
CPU seconds ({budget['completed_cpu_seconds']/3600:.4f} CPU hours)**. It includes
native admissions, abandoned invariant calls, the closed contraction prefix
and its replay, initial trials, stopped resumes, the conservative 45-second
diagnostic charge, the final four workers, and final replay. The complete
ledger is `control/native_budget_aggregate.json`.

The original ceilings remain 250,000 pentachora, 8 GiB, at most 16 seeds,
and 14,400 cumulative CPU seconds for this native-control campaign. This
batch used four distinct seeds. Each of the four concurrent final search
workers had a 2 GiB address-space cap. The final replay had a separate
2 GiB cap and 1,200 CPU-second allocation. The CPU and seed ceilings were
not exhausted. Source construction and the other exact combinatorial audits
are outside this native-budget ledger.

The first two historical RSS samples were invalidated because their monitor
used a namespace-local process ID against the host process table. Their
waited-child CPU accounting remains valid. Later workers used corrected
process identification and kernel limits; the package makes no claim that
the unavailable historical RSS samples were measured successfully.

Bare tight decoding loses Regina's cached manifold information and can make
the next skeleton query redo expensive three-sphere recognition. The saved
checkpoint therefore remains a storage encoding; a resumed search restores
its working state by replaying the preceding legal history from the original
input. That restoration also verifies the prefix before extending it.

## Remaining acceptance gate

The criterion remains a full source-bound replay to one pentachoron,
an actual isomorphism with `Example4.ball()`, and signature `baa`.
No saved path in this batch satisfies it. The bounded outcome gives no
exoticness inference. The stored cover is a homotopy four-sphere by the
separate group/homology proof, while the requested PL standardization
certificate remains outstanding. No closed X0 has been built.

For the release's authoritative combined status, read
`control/MILESTONE_STATUS.json`. Earlier source-stage reports and final
review snapshots retain their original scope and timestamps; their pending
status statements are superseded by that combined status and this result.

To repeat the best native replay, use the pinned environment in a working
copy, aggregate all earlier native charges, and run `native_control.py replay
--seed {best['seed']} --work-source` through `supervise_native.py`. Running a
replayer refreshes its receipt, so first verify the untouched release with
`python verify_package.py`.
'''
    (ROOT/'NATIVE_BATCH_RESULT.md').write_text(text)
    print(json.dumps({'status':result['status'],'best_replayed_seed':best['seed'],'endpoint_pentachora':replay['endpoint_pentachora'],'charged_cpu_seconds':budget['completed_cpu_seconds']}))

if __name__=='__main__':main()

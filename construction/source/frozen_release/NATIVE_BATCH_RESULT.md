# Finished native batch: exact outcome and replay

**Result:** `CONTROL_NOT_STANDARDIZED_WITHIN_DECLARED_BATCH`.
The best complete replay ends at **2,072 pentachora**.
Its `replay_passed` value is true and its `standardization_passed` value is
false. The original marked control,
geometric identification, direct trivial-group certificate, and integral
homology certificate are unchanged.

## Accepted source-to-endpoint path

The path begins with the immutable 103,744-pentachoron marked master.
The independently replayed 3,307 full-link contractions produce a closed,
unmarked 44,728-pentachoron abstract simplicial complex of the same oriented
PL type. All markings remain on the master.

Seed 1 removes the interior of its embedded pentachoron at
zero-based index 8805, with vertices
`[66, 142, 241, 242, 699]`. The punctured input has 44,727 pentachora
and one standard three-sphere boundary. This is a puncture of the recorded
source, not a substitute model.

| Step | Pentachora after the step | Verification |
| --- | ---: | --- |
| Immutable closed master | 103,744 | Marked assembly and native PL-manifold admission |
| 3,307 closed full-link contractions | 44,728 | Independent standard-library replay |
| Recorded embedded-simplex puncture | 44,727 | Source facet labels and native gluing index |
| Initial local native simplification | 9,164 | Exact numbered hash reproduced |
| 9 guarded 2-4 moves | 9,182 | Every native legality check replayed |
| Full native simplification | 2,072 | Exact numbered hash and final tight encoding reproduced |

The last numbered gluing hash is
`4f908feff00b9893668de1afeffa712d2311b995588bfdb23df9983d10f93316`.
The final replay checks the exact source file, master facets, contraction
prefix, puncture labels, Regina version, each move, each per-block macro
policy, each endpoint hash/size, and final tight encoding. Its receipt is
`control/standardization_seed1_replay.json`.

Native simplification is a version-pinned software trust boundary. These
receipts are macro replays, not an exported elementary Pachner transcript.
The final native replay took 364.967 wall seconds and
was charged 365.187 CPU seconds by the external supervisor.

## Four-seed results

| Seed | Replayed initial local prefix | Last saved complete endpoint | Saved macro blocks | Separate complete-path replay |
| ---: | ---: | ---: | ---: | --- |
| 0 | 9,090 | 2,252 | 2 | No separate final replay |
| 1 | 9,164 | 2,072 | 2 | Passed |
| 2 | 9,708 | 9,708 | 1 | No separate final replay |
| 3 | 9,057 | 2,113 | 2 | No separate final replay |

Seeds 0 and 3 finished their current worker normally. Seed 1 saved its
2,072-cell endpoint before a later trial reached its CPU cap. That later
candidate was discarded. Seed 2's full macro reached a reported objective
of 2,395, but it was interrupted before returning; the accepted checkpoint
remains 9,708. In-progress objective counts are not promoted to certificates.
Every initial local prefix was replayed during restoration. Only the final
receipt listed above certifies the entire best saved path.

## Budget and operational lessons

The charged native-control total is **5802.218
CPU seconds (1.6117 CPU hours)**. It includes
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
--seed 1 --work-source` through `supervise_native.py`. Running a
replayer refreshes its receipt, so first verify the untouched release with
`python verify_package.py`.

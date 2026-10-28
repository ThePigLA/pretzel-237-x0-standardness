# Independent review of final claims

**Review snapshot:** 2026-10-26, 19:35 UTC, before the complete final seed-1
replay receipt was available.  This review reads the current README, milestone,
geometric proof and reviews, exact group and contraction receipts, and saved
native checkpoints.  It does not rerun the heavy geometry or native searches.
Fourteen current byte-hash links among the milestone, group certificate,
contraction prefix, and native work source were checked and matched.

**Decision:** the construction, direct fundamental-group, integral-homology,
and unmarked PL-contraction claims pass with their stated trust boundaries.
The stored cover may also be called a **homotopy four-sphere** by the argument
below.  **A PL standardization to the four-sphere has not been established.**

## Claim matrix

| Claim | Decision | Evidence and exact boundary |
| --- | --- | --- |
| The source is the specified signed and oriented \(P(-2,3,7)\), with its marked meridian, longitude, and puncture | **PASS, under declared native PL semantics** | Direct rational-tangle wiring, exact peripheral arrays, an even combinatorial endpoint isomorphism with peripheral matrix \(+I\), and the explicit native handedness calibration. Intermediate native diagram construction and PL transport retain their stated software trust boundary. |
| The actual 103,744-pentachoron cover and 51,872-pentachoron orbit quotient are closed oriented PL four-manifolds | **PASS, under native manifold recognition** | The independent assembly checker regenerates the actual cells and all facet/covering maps. Both exact native exports have current validity, closedness, connectedness, orientability, and nonideal admissions. The cover's \(f\)-vector is \((4310,64796,216104,259360,103744)\), so \(\chi=2\). |
| The stored marked cover realizes the intended \(m=2,n=1\) construction | **PASS, with the explicit geometric theorem dependency** | `GEOMETRIC_CONSTRUCTION_ADMISSION.json` binds the corrected mixed profile, enlarged \(N^+\), source-positive puncture, complete winding, actual bridge/cap, and reviewed localization argument. The global framed normal chart is supplied by relative PL normal-neighborhood and boundary-extension results. An enumerated matrix atlas for that global chart is not claimed. |
| The full branch set has the required standard local models | **PASS within the reviewed geometric construction** | The review covers source vertex and edge models, products, relative equivariant ball replacements, refinement collars, cap cone vertices, and actual port gluing. Fixed-triangle normal circles alone would not establish this conclusion; the package explicitly uses the stronger geometric provenance. |
| \(\pi_1\) of the immutable cover is trivial | **PASS, exact independent replay** | All 64,796 edges and 216,104 triangles are included. The checked tree presentation has 60,487 generators. The trace has 60,487 one-occurrence Tietze eliminations and 13 retained-relator shortenings, with 470,791 rewritten relators checked, ending at the empty presentation. This does not depend on a native simplification. |
| Integral homology of the cover is \((\mathbb Z,0,0,0,\mathbb Z)\) | **PASS, theorem deduction from checked premises** | Trivial \(\pi_1\), admitted connected closed orientation, Poincare duality, the universal coefficient theorem, and independently counted \(\chi=2\). The group receipt applies directly to the cover; it is not an independent group calculation for the quotient. |
| The cover is a homotopy four-sphere | **PASS, theorem deduction** | Hurewicz and homological Whitehead, with the explicit argument below. A finite simplicial complex is a CW complex, so the required hypotheses hold. |
| The 44,728-pentachoron descendant has the same oriented PL type as the marked master | **PASS, unmarked** | All 3,307 full-link contractions and oriented relative four-ball replacements were independently replayed, including 1,190,001 full-face overlap checks. Marks and deck structure remain on the immutable master. Queue exhaustion is not global contraction-minimality. |
| The saved seed-1 endpoint has 2,072 pentachora | **PASS as a saved checkpoint; final replay pending at this snapshot** | The checkpoint contains two accepted macro blocks and a tight encoding. Its last numbered hash is `4f908feff00b9893668de1afeffa712d2311b995588bfdb23df9983d10f93316`. The older resume-prefix receipt verifies only the earlier 9,164-cell state. It cannot certify the later 2,072-cell block. |
| The cover is PL homeomorphic to standard \(S^4\) | **NOT ESTABLISHED** | No complete saved replay ends in the accepted standard-ball model. Successful replay to a 2,072-cell punctured endpoint would certify that partial simplification, while standardization would remain unresolved. |
| The bounded searches give evidence of exoticness, or permit the guarded closed \(X_0\) stage | **NOT ESTABLISHED / GATE REMAINS CLOSED** | A bounded search stopping without standardization has no such mathematical implication. The package correctly keeps `exoticness_inference` null and records no closed \(X_0\) construction. |

## Why the homotopy-sphere statement is justified

Let \(M\) be the admitted closed cover.  The exact certificate gives
\(\pi_1(M)=0\), and the integral homology calculation gives
\(H_2(M)=H_3(M)=0\).  Applying Hurewicz successively gives
\(\pi_2(M)=0\) and then \(\pi_3(M)=0\).  Thus the Hurewicz map
\(\pi_4(M)\to H_4(M)\cong\mathbb Z\) is an isomorphism.  Represent a generator
by \(f:S^4\to M\).  This map induces an isomorphism on every integral homology
group.  Both spaces are simply connected CW complexes, so the homological
Whitehead theorem makes \(f\) a homotopy equivalence.

This uses Hatcher, *Algebraic Topology*, Theorem 4.32 and Corollary 4.33,
pages 366–367, checked in the author's
[official text](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf).
It produces a homotopy equivalence, not a PL homeomorphism.

## What would close the PL standardization gate

The acceptance rule in the milestone is mathematically sufficient:

1. Retain the independently replayed closed contraction ancestry from the
   immutable master to the abstract 44,728-pentachoron descendant.
2. Remove the interior of the recorded embedded pentachoron.  Since it is
   an actual simplex of an abstract simplicial complex, its closure is a
   PL four-ball with standard PL three-sphere boundary.
3. Replay every guarded move and version-pinned native simplification macro
   from that exact punctured source, matching the complete numbered gluing
   hashes and final tight encoding.
4. Identify the endpoint with `Example4.ball()` by an explicit native
   combinatorial isomorphism, with the accepted ball signature `baa`.

Then the punctured descendant is a PL four-ball.  Gluing back the removed
simplex gives standard PL \(S^4\), since a PL homeomorphism of the common
three-sphere boundary extends over a PL ball by coning after subdivision.
This route needs no four-dimensional Schoenflies assumption.

Native macro replay still uses the pinned native PL-move semantics.  It
must continue to be described as a source-bound, version-pinned macro
certificate rather than a portable list of individual Pachner moves.
The code and README already make this distinction.

## Presentation corrections for the final frozen batch

These are final-status and scope clarifications; no mathematical failure was
found in the checked construction, group, homology, or contraction claims.

1. **Freeze the batch status after the final replay.** The README currently
   says standardization is in progress, and the milestone has
   `search_finished_for_this_batch=false` with empty replay arrays.  For a
   completed bounded-batch deliverable, update these fields to the actual
   completed status and bind the final seed-1 replay receipt if it passes.
   Keep `standardization.passed=false` unless the standard-ball gate is met.

2. **Distinguish saved checkpoints from replayed endpoints.** The current
   saved counts are seed 0: 2,252; seed 1: 2,072; seed 2: 9,708; seed 3:
   2,113.  The interrupted seed-2 full macro contributes no accepted new
   endpoint.  A resume-prefix replay proves only the blocks it names.  At
   this snapshot seed 1's old prefix receipt ends at 9,164 cells, so the
   later 2,072-cell result must await the complete replay now being run.

3. **Make the historical source-audit scope explicit.** The earlier
   `diagram/AUDIT_SOURCE.md` still has a heading about the remaining
   faithful-\(X_2\) gate and language saying global identification remains
   separate.  Its source-level caution was appropriate when written.  The
   later `GEOMETRIC_CONSTRUCTION_ADMISSION.json` now discharges geometric
   identification under the stated theorem trust basis.  Point readers to
   that later admission; PL standardization alone remains the outstanding
   requested gate.  This is a status clarification, not a reason to alter
   the immutable source or old receipt data.

4. **State the strongest justified result precisely.** Suggested wording:
   “The stored marked construction is a homotopy four-sphere, under the
   stated native and geometric theorem trust basis.  The direct fundamental
   group certificate and integral homology deduction pass.  The bounded
   search has not supplied a PL standardization.”  The homotopy claim may
   be added with the short proof above; it should not be presented as an
   endpoint recognition result.

5. **Refresh the final manifest and budget aggregate.** The current
   aggregate charges 5,437.0313 completed CPU seconds and reserves 1,200
   seconds for final replay, within the 14,400-second ceiling.  It counts
   four distinct seeds and explicitly includes interrupted/superseded work.
   Replace the reservation with the final measured replay cost.  The
   16-seed and four-CPU-hour values are ceilings, not exhausted-search
   claims.  Preserve the explicit invalidation of the first two RSS
   samples instead of advertising an all-run measured peak.

## Snapshot bindings

| Item | SHA-256 |
| --- | --- |
| Immutable marked master file | `af21a7dd107b765d5fef0a257dcb4df6d6df8689c594f584e4318d15c49aa709` |
| Master facets | `65de21dd16c5c11fa880e6f0c80a3cbf25a4780ba3d04ba7a6ab838ef37e5587` |
| Geometric admission receipt | `4dd03a21e6693930943832a44636bc61a590ead00e6c09440910e7ac710e053a` |
| Direct group replay receipt | `4259b4e7e794a259fe21088dbdd226352c66005792f23565e5d607bdc95d6c8a` |
| Closed contraction replay receipt | `0c9ebb6bdc29fc10753de47637d3a88f099c7b4666a5b92e1960bb0b0571d4b3` |
| README version read for this review | `23cca14d69dfd3d30c361334901c1259610cd686e80a95802ccbecd1e5eb9026` |
| Milestone version read for this review | `a9aa1f3829741a5e675ea698613097b5795b6601c396ebb1ff5fffe3cbb04d7f` |

Later final-status updates to README or the milestone should retain the
immutable source and proof bindings above.  A later complete native replay
receipt supersedes only the pending replay status in this review; it does
not turn a nonstandard-size endpoint into a standardization certificate.

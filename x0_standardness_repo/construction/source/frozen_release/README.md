# P(-2,3,7): marked even-twist control

This package continues the accepted trefoil-control work with the new knot
`J = P(-2,3,7)`, positive roll `n = 1`, and even twist `m = 2`.
It contains a concrete signed source, a marked branched double cover, the
actual roll and cap meshes, a closed 103,744-pentachoron control, and exact
verification records. The original marked master is kept separate from its
unmarked simplification descendants.

## Released batch status

The geometric source and assembly checks have passed. The identification
with the intended twist-roll construction includes an explicitly stated PL
normal-neighborhood theorem step. It is not advertised as a finite matrix
replay of a global normal-disk-bundle trivialization.

The independent group calculation has passed: the immutable marked source
has trivial fundamental group and integral homology `(Z,0,0,0,Z)`; it is a
homotopy four-sphere. The declared four-seed native batch is finished.
A complete native replay from the certified source reproduces the best
punctured endpoint at **2,072 pentachora**. The reference-ball endpoint was
not reached, so the status is
`CONTROL_NOT_STANDARDIZED_WITHIN_DECLARED_BATCH`.
**The PL-standardization gate remains open. No closed X0 candidate has
been constructed.** The source-level
`status` field inside the immutable master records its creation state; the
separate admission and milestone receipts record subsequent verification.
The authoritative combined status is `control/MILESTONE_STATUS.json`.
Earlier source-stage reports, including `diagram/AUDIT_SOURCE.md`, retain
their original scope. Its formerly pending geometric-identification gate
is discharged by `control/GEOMETRIC_CONSTRUCTION_ADMISSION.json` under the
stated theorem basis. Final review files record dated snapshots; their
transient search-status statements are superseded by the combined status.

## What is concrete

| Piece | Saved size or result |
|---|---:|
| Signed rational-tangle diagram | 12 crossings, one component |
| Initial equivariant Seifert/Montesinos cover | 7,440 tetrahedra |
| Certified compressed source | 628 tetrahedra |
| Source exterior after marked native transport | 3 ideal tetrahedra |
| Punctured tube with rectangular adaptor | 1,576 tetrahedra |
| Exterior with rectangular adaptor | 2,334 tetrahedra |
| Central roll collar | 696 tetrahedra |
| Even roll bulk | 73,696 pentachora |
| Actual puncture-to-cap bridge | 24,224 pentachora |
| m=2 cap | 5,824 pentachora |
| Closed marked control | 103,744 pentachora |
| Orbit-cell quotient | 51,872 pentachora |
| Replayed unmarked contraction descendant | 44,728 pentachora |
| Complete native replay after a recorded puncture | 2,072 pentachora |

The marked control has f-vector
`(4310, 64796, 216104, 259360, 103744)` and Euler characteristic 2.
Its quotient has f-vector `(2208, 32554, 108156, 129680, 51872)` and the same
Euler characteristic. Regina admits both as valid, connected, orientable,
closed, non-ideal four-manifolds, including its recognition of vertex links.
These facts do not themselves recognize a four-sphere.

The full fixed subcomplex is a sphere with f-vector `(106,312,208)`.
The separate twin sphere has f-vector `(10,24,16)` and meets the fixed
sphere in two vertices. Every recorded normal link of a fixed triangle is
a circle with free deck half-turn. The all-strata local branch proof also
uses the actual product, collar and cap geometry, including the cap cone
vertices; it does not infer local flatness from those triangle links alone.

## Signed source and peripheral marking

The diagram is constructed directly as the numerator closure of
`T(-1/2) + T(1/3) + T(1/7)`. The direct signed wiring agrees with the
Spherogram rational-tangle construction. The package records the PD code,
orientation, an exact polygonal embedding, a marked arc and puncture cube,
meridian, zero longitude, and a checked Wirtinger commutator certificate.
The diagram determinant is 1 and its Alexander polynomial is
`t^10 - t^9 + t^7 - t^6 + t^5 - t^4 + t^3 - t + 1`.

The branched source is built explicitly from equivariant Seifert pieces,
not selected by a knot-name lookup. Its quotient is recognized as S3, and
its branch exterior is compared with the specified diagram's marked
exterior. The selected exact combinatorial isomorphism preserves both
peripheral directions; its two transported normal-flow arrays agree
entry for entry with the diagram arrays. The decorated endpoint signature
is `dLQacccjnbs_BaRtb`.

One sign change is recorded openly. SnapPy's importer reverses the initially
extracted meridian and leaves the longitude unchanged, so the diagram chart
is `(u_old, -v_old)`. Its right-handed tetrahedron convention is opposite
the usual ordered-simplex convention. After calibrating both conventions,
the physical diagram ambient orientation agrees with the original master.
The positive puncture equator is `[56,150,58,148]`, with poles `[204,228]`.
The four-dimensional orientation is anchored to an unchanged source
tetrahedron times a positive time edge and propagated through every gluing.

## Roll, twist and cap

The positive roll is the actual longitudinal shear in the marked torus
collar. Exact rational overlays, affine inverses, face checks and a finite
source-profile transfer supply the central collar. All ten saved source
seam witnesses are located in actual four-dimensional facets, including
the wrap from time index 3 to time index 0. The full winding is +1.

The two even caps use
`beta_m(theta,t) = (theta - m*t/2, t)`. The m=0 and m=2 templates share
the same boundary mesh, but keep different full attachment maps. In the
eight-angle/four-time convention their lifted displacements are `[0,4]`
and `[-8,4]`. The m=2 cap has one negative upstairs turn and two negative
downstairs turns. Endpoint parity is not substituted for this path.

The bridge matches the actual 96-tetrahedron source boundary to the
1,792-tetrahedron shared cap boundary. Its 8,752-tetrahedron intermediate
refinement and both four-dimensional collars are checked by exact
barycentric maps, face incidences, volume partitions and deck maps.

The finite localized transfer uses the actual spherical family R:
`F2(x,q)=R_(1-q)(x)` and `Phi(x,q,t)=R_min(1-q,t)(x)`.
The data contain 1,792 spherical maps and 3,584 rational four-polytopes,
their affine inverses, face incidences and seam equations. The integer
lift correction from `theta+q-1` to `theta+q` is retained explicitly.

The global localization proof first applies a simultaneous longitudinal
and meridional profile transfer. It then uses the enlarged tube, including
the saved N adaptor, whose outer grid is disjoint from the puncture. This
separation is necessary: the original cut torus meets the puncture equator.
Relative PL normal-neighborhood and boundary-extension results supply a
framed pair chart on the enlarged tube. A joint PL normal-disk loop and
a space-and-time Alexander comparison identify its boundary loop with
the saved cap map. The full proof and its independent review state this
theorem dependency explicitly.

## Verification map

| Subject | Main files |
|---|---|
| Diagram | `diagram/pretzel_237_diagram.json`, `verify_diagram.py` |
| Signed source and exact marked isomorphism | `diagram/AUDIT_SOURCE.md`, `verify_marked_isomorphism.py`, `source_audit_certificate.json` |
| Initial cover and quotient | `seifert/` |
| Source compression and exact cuts | `coarsened/`, `peripheral/`, `verify_source_geometry.py` |
| Actual torus adaptor | `adaptor/`, `verify_torus_adaptor.py` |
| Roll and mixed profile | `central_roll/`, `verify_central_roll.py`, `verify_mixed_profile_transfer.py` |
| Global geometric identification | `central_roll/LOCALIZATION_LEMMA.txt`, `cap/REVIEW_BRANCH_AND_LOCALIZATION.txt` |
| Combined geometric admission | `control/GEOMETRIC_CONSTRUCTION_ADMISSION.json` |
| Cap, bridge and finite transfer | `cap/README_CAP.txt`, `cap/verify_caps.py`, `cap/verify_bridge.py`, `cap/verify_transfer.py` |
| Frozen marked control and native inputs | `control/marked_control.json.gz`, `control/cover_regina.json.gz`, `control/quotient_regina.json.gz` |
| Independent whole-assembly audit | `diagram/verify_control_assembly.py`, `diagram/control_assembly_verification.json` |
| Unmarked contraction prefix | `unmarked_contractions_v1/README_CONTRACTIONS.md`, `certificate.json.gz`, `independent_replay.json` |
| Direct master group and integral homology | `group_audit/README_GROUP_AUDIT.md`, `group_audit/independent_replay.json` |
| Native work and budget records | `native_control.py`, `supervise_native.py`, `control/native_budget_ledger.json` |
| Finished native batch and complete best-path replay | `NATIVE_BATCH_RESULT.md`, `control/NATIVE_BATCH_RESULT.json`, `control/standardization_seed1_replay.json` |
| Authoritative combined milestone | `control/MILESTONE_STATUS.json` |

The complete source and output hashes are in `MANIFEST_SHA256.json`.
The independent whole-assembly checker regenerates all 103,744 pentachora
without importing the producer, checks every cover and quotient facet
pairing and covering square, and verifies the entire fixed subcomplex,
the twin intersection, the source orientation anchor and source carriers.

## Certified unmarked prefix

The 3,307 edge contractions satisfy the full simplicial link condition in
every dimension. The independent replay checks each relative PL four-ball
replacement, identical oriented boundaries, and 1,190,001 full-face
overlap conditions. It yields f-vector
`(1003,25367,91458,111820,44728)` and preserves the oriented PL type.

The local proof uses a PL three-sphere vertex link and its three-ball
antistar; it does not assume a four-dimensional Schoenflies theorem.
The source and endpoint are linked by an explicit ordered trace. This
descendant is unmarked; the deck and surfaces remain on the frozen master.
Local queue exhaustion is not a claim of global contraction minimality.

## Independent fundamental group and integral homology

The full two-skeleton of the immutable marked master gives 60,487 generators
after a spanning tree and 216,104 triangle relations. The exact trace has
60,487 one-occurrence Tietze eliminations and 13 strict relator shortenings,
ending at the empty presentation. The independent standard-library replayer
rebuilds the skeleton and tree and checks 470,791 rewritten relators. It
imports neither producer, Regina nor a group algebra package.

There is a short explanation of the two-generator tail. Let `a=14941`,
`b=15470` and `T(x)=b^-1 x b`. Four retained original residual relators give
`T(a^2)=a`, `T(a^3)=a^3`, `T^3(a)=a`, and a relation which reduces to
`b^-1` when `a=1`. The first two imply `T(a)=a^2`, so `a^3=1`. The third
then gives `a^7=1`. Hence `a=1`, and then `b=1`. The exact words and their
indices are supplied in the group README and independently reviewed in
`cap/GROUP_TAIL_REVIEW.txt`.

Thus pi1 is trivial. For the admitted closed connected oriented four-manifold,
Poincare duality and the universal coefficient theorem give H3=0 and
torsion-free H2. The independently recomputed Euler characteristic 2 forces
H2=0. This proves integral homology `(Z,0,0,0,Z)` independently of the native
simplification attempt. It is not a PL standardization argument.

These results also show that the stored cover is a **homotopy four-sphere**.
Successive applications of Hurewicz give pi2=pi3=0 and then
pi4=H4=Z. A generator gives a map S4 to the cover inducing every integral
homology isomorphism. Both are simply connected CW complexes, so homological
Whitehead makes that map a homotopy equivalence. This conclusion has the
same admitted PL-manifold premise as the homology deduction; it does not
close the separate PL-standardization gate. See `REVIEW_FINAL_CLAIMS.md`
and Hatcher, Theorem 4.32 and Corollary 4.33, for this implication.

## Finished native batch

The best accepted path is now replayed in full. Starting from the closed
44,728-pentachoron descendant, seed 1 removes the embedded pentachoron at
index 8805 with vertices `[66,142,241,242,699]`. The 44,727-cell punctured
input reduces to 9,164 cells by a local macro; nine guarded 2-4 moves give
9,182 cells; a full native macro reduces this to 2,072. The separate replay
matches both completed macro hashes and the final tight encoding exactly.
It took 364.967 wall seconds and was charged 365.187 CPU seconds.

| Seed | Last saved completed endpoint | Complete-path verification |
| ---: | ---: | --- |
| 0 | 2,252 | Initial local prefix replayed; no separate replay of the later full block |
| 1 | 2,072 | Both blocks and all nine guarded kicks replayed successfully |
| 2 | 9,708 | Local-only saved path replayed during restoration; later full macro interrupted |
| 3 | 2,113 | Initial local prefix replayed; no separate replay of the later full block |

Seed 2's in-progress objective of 2,395 was never a saved completed macro
endpoint and is excluded from the accepted results. The native macro
certificate remains distinct from an exported elementary Pachner trace.
Its source, version, puncture, move legality, per-block policy, numbered
hashes and final encoding are all bound by the final replay receipt.

The charged native-control total is **5,802.218 CPU seconds (1.6117 CPU
hours)**, including unsuccessful and interrupted work, prefix production
and replay, the diagnostic allowance, and final replay. Four distinct seeds
were used. The original four-CPU-hour and sixteen-seed ceilings were not
exhausted. The detailed outcome, resource scope and reproduction caveats
are in `NATIVE_BATCH_RESULT.md`; all run charges are in
`control/native_budget_aggregate.json`.

The required success criterion remains a replay to the one-pentachoron
`Example4.ball()` with signature `baa`. The finished batch has no such
endpoint. Its bounded outcome supports no exoticness inference and leaves
the closed-X0 stage guarded.

## Reproduction and trust basis

The recorded environment is Python 3.12, Regina distribution 7.4.1
(API reports 7.4), SnapPy 3.3.2, Spherogram 2.4.1, SymPy 1.14.0, and
psutil 7.2.2. See `requirements.txt`. Most independent mesh, sign,
peripheral, contraction and assembly checks use only the standard library.
Producer regeneration and native checks require the listed packages.

Run checks on a working copy: some scripts refresh their verification
receipts. Representative commands, from this directory, are:

```sh
python verify_package.py
python diagram/verify_diagram.py
python diagram/verify_marked_isomorphism.py
python verify_source_geometry.py --mutations
python verify_torus_adaptor.py
python verify_central_roll.py --mutations
python verify_mixed_profile_transfer.py
python verify_twist_localization_inputs.py
python cap/verify_caps.py
python cap/verify_bridge.py
python cap/verify_transfer.py
python cap/verify_review_inputs.py
python diagram/audit_source_chain.py
python diagram/verify_control_assembly.py
python diagram/replay_unmarked4.py
python diagram/replay_control_group.py
python cap/verify_group_tail.py
```

The source identification trusts the pinned oriented native diagram
constructor, finite-vertex removal and simplification. Its input markings
and exact endpoint isomorphism are checked independently. A portable
elementary-move transcript for every intermediate native three-dimensional
triangulation is not supplied. The geometric identification also uses the
stated relative PL normal-neighborhood theorem. Native simplification
macros, if used, are version-pinned and replayed by numbered gluing hashes;
they are not relabelled as portable elementary Pachner transcripts.

The original limits remain 250,000 pentachora and 8 GiB, with at most
16 seeds and 4 cumulative CPU-hours in the first native batch. Records
include unsuccessful or deliberately superseded work. There is no
return to the earlier core42 campaigns and no closed X0 construction
before the X2 gate is resolved.

## Primary mathematical and software references

- Spherogram's explicit pretzel construction:
  https://snappy.computop.org/spherogram.html
- Kuhrman, *More Exotic RP2-knots and Homotopy Spheres*, sections 2.1–2.3:
  https://arxiv.org/html/2507.03798v1
- Rourke–Sanderson, *Block Bundles I*, section 4:
  https://webhomes.maths.ed.ac.uk/~v1ranick/papers/rs003.pdf
- Nevo, *Higher Minors and Van Kampen's Obstruction*, Theorem 1.4 and section 5.1:
  https://arxiv.org/pdf/math/0602531
- Hatcher, *Algebraic Topology*, Theorem 4.32 and Corollary 4.33:
  https://pi.math.cornell.edu/~hatcher/AT/AT.pdf

The sources explain the construction and PL preservation statements.
They do not substitute for this package's source, framing, assembly or
standardization checks.

# Fundamental group and integral homology of the immutable control

The complete two-skeleton of `control/marked_control.json.gz` has a verified
presentation that reduces to the empty presentation.  Thus the actual stored
closed control has

\[
\pi_1(M)=1.
\]

Combining this exact group certificate with the source-bound admission as a
connected closed oriented PL four-manifold and the independently recomputed
Euler characteristic \(\chi(M)=2\) gives

\[
H_i(M;\mathbb Z)=
\begin{cases}
\mathbb Z,& i=0,4,\\
0,& i=1,2,3.
\end{cases}
\]

This calculation uses the immutable 103,744-pentachoron master directly.  It
does not use a native simplification or the unmarked contraction descendant.
It certifies the group and integral homology of the stored control; a PL
standardization and the separate identification with the intended marked
geometric construction require their own certificates.

## Source and presentation

The master facet hash is

`65de21dd16c5c11fa880e6f0c80a3cbf25a4780ba3d04ba7a6ab838ef37e5587`.

Its exact source file hash is

`af21a7dd107b765d5fef0a257dcb4df6d6df8689c594f584e4318d15c49aa709`.

All vertices, edges, and triangles are obtained from the literal abstract
simplicial facets.  Their counts are

\[
(f_0,f_1,f_2)=(4310,64796,216104).
\]

A verified spanning tree has 4,309 edges.  Collapsing this tree gives one
generator for each remaining oriented edge, hence

\[
64796-4310+1=60487
\]

generators.  Edges are oriented from their smaller vertex label to their larger
one.  Each sorted triangle \([a,b,c]\) contributes the relation
\(e_{ab}e_{bc}e_{ac}^{-1}\), with tree edges omitted.  Every one of the 216,104
source triangles is included.  No quotient group relation is imposed.

The first stage performs 60,485 exact one-occurrence Tietze eliminations and
leaves two generators and 324 relators.  The second stage performs 13 strict
relator shortenings, each using a different retained relator, followed by the
last two generator eliminations.  The endpoint has no generators and no
relators.

For a relator \(A g^\epsilon B\) in which \(g\) occurs once, the elimination
solves

\[
g=(BA)^{-1}\quad(\epsilon=1),\qquad
g=BA\quad(\epsilon=-1),
\]

substitutes this exact word into every other relator, and removes the generator
and its defining relator.  Each shortening replaces a subword \(U\) with
\(V^{-1}\) using a retained relation \(UV=1\), allowing cyclic rotation and
inversion of that defining relation.  Both operations preserve the presented
group.  Free and cyclic reductions, empty-relator deletion, and duplicate
relator deletion preserve the same normal closure.

## A short independent explanation of the two-generator tail

Let \(a\) denote original generator 14941 and \(b\) original generator 15470.
Four of the surviving relators already force the group to be trivial.  The
indices below are zero-based indices in the prefix certificate's ordered
`relators` array.

| Index | Relator equal to the identity |
| ---: | --- |
| 47 | \(b^{-1}a^2ba^{-1}\) |
| 46 | \(b^{-1}a^3ba^{-3}\) |
| 262 | \(b^2a^{-1}b^{-3}ab\) |
| 26 | \(b^{-1}ab^{-1}ab^{-1}a^{-1}b^2a^{-1}\) |

Write \(T(x)=b^{-1}xb\).  The first two relations give

\[
T(a^2)=a,\qquad T(a^3)=a^3.
\]

Because conjugation is an automorphism,

\[
T(a)=T(a^2)^{-1}T(a^3)=a^2.
\]

Consequently \(a=T(a^2)=a^4\), so \(a^3=1\).  A cyclic rotation of relator 262
is \(b^3a^{-1}b^{-3}a\); it gives \(T^3(a)=a\).  On the other hand,
\(T(a)=a^2\) gives \(T^3(a)=a^8\).  Thus \(a^7=1\), and

\[
a=a^7(a^3)^{-2}=1.
\]

With \(a=1\), relator 26 reduces to \(b^{-1}\), so \(b=1\).  The other
surviving relators are then automatically satisfied.  This short algebraic
argument was obtained independently of the logged shortening search; the
machine certificate replays the actual 15-move tail instead.

## Independent replay

`diagram/replay_control_group.py` uses only the Python standard library.  It
does not import the producers, `group_cells`, Regina, or a group algebra
package.  It performs the following source-bound checks:

1. Reconstruct the entire source two-skeleton from the immutable facets and
   compare it with the stored skeleton.
2. Check that the selected tree is acyclic and spans all source vertices, that
   every remaining edge has its declared generator, and that every triangle
   contributes its exact oriented boundary word.
3. Reconstruct the evolving presentation and replay all 60,487 generator
   eliminations and all 13 retained-relator shortenings, checking every
   defining relation and every substitution.
4. Compare both intermediate and final presentations exactly with their
   certificate endpoints.  The final presentation is empty.
5. Recompute the source's complete \(f\)-vector, closed codimension-one
   incidence, and oriented top-dimensional boundary chain.  Bind the separate
   closed PL-manifold admission to the exact source geometry.

The replay checked 470,791 rewritten relators.  The source \(f\)-vector is

\[
(4310,64796,216104,259360,103744),
\]

whose alternating sum is 2.

The group conclusion itself needs no manifold recognition: a connected
finite simplicial complex has the fundamental-group presentation above.  The
integral homology deduction additionally uses the already admitted closed
oriented PL-manifold premise.  That admission remains a native mathematical
trust boundary; it is not replaced by a homology test in this audit.

## Integral homology deduction

Trivial fundamental group gives \(H_1(M;\mathbb Z)=0\).  The universal
coefficient theorem then gives \(H^1(M;\mathbb Z)=0\), and Poincare duality gives
\(H_3(M;\mathbb Z)=0\).

For degree two, the universal coefficient sequence has zero
\(\operatorname{Ext}(H_1(M;\mathbb Z),\mathbb Z)\) term.  Therefore
\(H^2(M;\mathbb Z)\) is free abelian, and Poincare duality identifies it with
\(H_2(M;\mathbb Z)\).  Thus \(H_2\) has no torsion.  Finally,

\[
\chi(M)=2+\operatorname{rank}H_2(M;\mathbb Z)=2
\]

forces \(H_2(M;\mathbb Z)=0\).  Connectedness and orientation supply
\(H_0\cong H_4\cong\mathbb Z\).

## Files and resources

| File | Contents |
| --- | --- |
| `source_binding.json` | Exact master, skeleton, presentation, and manifold-admission bindings |
| `skeleton.json.gz` | All source vertices, edges, and triangles |
| `initial_presentation.json.gz` | Spanning tree, edge generators, and every initial relator |
| `tietze_certificate.json.gz` | 60,485-move prefix and exact two-generator endpoint |
| `residual_tietze_certificate.json.gz` | 13 shortenings and two final eliminations |
| `producer_summary.json` | Prefix counts, limits, and resource use |
| `residual_producer_summary.json` | Short-tail counts and resource use |
| `independent_replay.json` | Successful source-to-empty-presentation replay receipt |

| Stage | Wall seconds | CPU seconds | Peak resident bytes |
| --- | ---: | ---: | ---: |
| Skeleton, presentation, and prefix producer | 10.7190 | 10.7141 | 326,676,480 |
| Residual producer | 0.0939 | 0.0939 | 41,627,648 |
| Independent replay | 7.8136 | 7.8077 | 474,210,304 |

The producers together used less than 11 seconds, within the requested
90-second producer budget.  Each process enforced an 8 GiB address-space limit.

To repeat the independent check from the project parent directory:

```sh
python j_source_v1/diagram/replay_control_group.py
```

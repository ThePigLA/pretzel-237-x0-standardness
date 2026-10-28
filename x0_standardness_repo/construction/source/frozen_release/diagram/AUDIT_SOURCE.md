# Audit of the signed source and attachment interfaces

The saved source is identified with the **specified signed, oriented diagram of
\(P(-2,3,7)\)**, subject to the native construction and move semantics stated
below. The exact endpoint map preserves both peripheral directions. The
diagram-positive meridian is the negative of the initially extracted meridian.
The physical ambient orientation agrees with the canonical Seifert master.

The audit corrected one metadata error: the imported tetrahedra have ordered
simplex orientation **minus** the master, but SnapPy's `right_handed` convention
is itself **minus** the usual ordered simplex convention. Their product is
plus. No coordinates, periods, simplices, puncture cycles, or cap winding were
changed by that correction. The dependent adaptor and cap-bridge hashes have
been rebound and their stored-data gates have been rerun.

The machine record is [source_audit_certificate.json](source_audit_certificate.json).
It contains exact local witnesses, the selected isomorphism, both chart checks,
canonical source links, primary-source hashes, and file digests. Reproduce it
with Python's standard library:

```sh
python diagram/audit_source_chain.py
python diagram/verify_marked_isomorphism.py
```

## Diagram and native peripheral directions

The source is the numerator closure of the ordered rational-tangle sum
\(T(-1/2)+T(1/3)+T(1/7)\). The direct signed wiring, PD code, knot orientation,
spatial embedding, meridian and zero longitude are independently checked in
`verify_diagram.py`. This construction is also the explicit pretzel example in
the [official Spherogram documentation](https://snappy.computop.org/spherogram.html).
No knot-name lookup or invariant-only identification is used.

The planar orientation is fixed by viewing the diagram from positive z, with
overpasses at positive z. At its marked arc the tangent is \((-1,0,0)\), and
the positive normal basis is \((e_y,-e_z)\). The determinant of those three
vectors is +1. Thus the recorded meridian is the right-hand positive meridian,
and the longitude follows the chosen knot direction.

SnapPy's pinned
[link-complement constructor](https://raw.githubusercontent.com/3-manifolds/SnapPy/42217b3867ce3a00d1125a1b6860db335684c333/kernel/kernel_code/link_complement.c)
defines its longitude in the forward link direction, chooses its meridian by
the right-hand rule, and corrects the longitude by the signed crossing sum.
Spherogram's KLP conversion preserves the oriented incoming/outgoing strands
and sends the negative crossings here to the kernel's counterclockwise case.
The crossing sum is -12, so the native zero-framing correction agrees with the
diagram's \(\lambda=W\mu^{12}\). This identifies peripheral directions and
framing up to the usual simultaneous basepoint conjugacy.

## Exact exterior comparison

The audit reads the saved `.tri` files directly. It does not reimport them
through SnapPy while checking them. It enumerates all possible images of
tetrahedron 0 and its four vertices, propagates each candidate over every face,
and verifies all gluing squares. Regina independently returned the same two
isomorphisms.

The selected map from `matched_marked_exterior.tri` to `diagram_exterior.tri` is:

| Source tetrahedron | Target tetrahedron | Images of vertices 0,1,2,3 |
|---|---:|---|
| 0 | 0 | 1,2,0,3 |
| 1 | 2 | 0,3,1,2 |
| 2 | 1 | 2,1,3,0 |

Every permutation is even. Its peripheral matrix, with columns the images of
source \(M,L\) in target \((M,L)\), is the identity. In fact both transported
normal-flow arrays equal the diagram arrays entry for entry. The more general
integer-potential test also succeeds: both potential witnesses on the six
global cusp vertices are zero. The second isomorphism has peripheral matrix
minus the identity; it is recorded but is not the selected marking.

## The sign correction and ambient orientation

The 756-tetrahedron coned source retains its original integer peripheral
arrays. The native importer changes exactly ten nonzero meridian entries by
negation and leaves every longitude entry unchanged. Exact rational integration
of the normal dual chains through the source torus chart gives:

| Curve | Raw source periods \((u_{old},v_{old})\) | Native source periods |
|---|---|---|
| Meridian | (0,1) | (0,-1) |
| Zero longitude | (1,0) | (1,0) |

The kernel uses positive entries for flow **entering** a cusp triangle. It
normalizes the pairing to \(L\cdot M=+1\), negating M when necessary. The audit
computes that pairing exactly as -1 before import and +1 after import.
Consequently the diagram chart is

\[
u=u_{old},\qquad v=-v_{old}.
\]

There is an independent geometric orientation check. The source N-cocycle
evaluates to +1 on the oriented fixed-core cycle beginning with the puncture
edge \([204,228]\). Thus this edge agrees with positive longitude. For each
successive edge of the raw equator \([56,148,58,150]\), the ordered join
tetrahedron \([204,228,e_i,e_{i+1}]\) has sign -1 against the master orientation.
Reversing the equator to \([56,150,58,148]\) makes all four signs +1. This is
the positive normal-disk orientation in the master.

The apparent conflict with the imported tetrahedron signs is resolved by the
kernel convention. Its documented
[cusp corner orders](https://raw.githubusercontent.com/3-manifolds/SnapPy/42217b3867ce3a00d1125a1b6860db335684c333/python/snap/peripheral/link.py)
are outward-positive for a kernel right-handed tetrahedron. Place a usual
positive simplex at \(0,e_1,e_2,e_3\) and truncate at one quarter of each edge.
For each of the four kernel corner orders, the determinant of an outward
vector followed by its two tangent vectors is exactly \(-1/64\). This proves
that kernel right-handedness is opposite standard ordered-simplex orientation.
All 584 compact imported tetrahedra have standard ordered-simplex sign -1
relative to the master. Combining the two signs gives the physical master
orientation, followed by the even endpoint isomorphism.

For positive longitude and positive normal angle, the local ambient orientation
is \(du\wedge dr\wedge dv\). The outward tube boundary therefore has orientation
\(-du\wedge dv\), while the exterior boundary has \(+du\wedge dv\). The corrected
chart has exactly these signs on every actual source triangle. Requiring
\(+du\wedge dv\) on the tube boundary would introduce an erroneous extra
ambient reversal.

## Puncture, torus adaptor and cap

The puncture is the actual star of the fixed edge \([204,228]\), with its eight
boundary triangles and four-edge normal equator. Its positive cycle is
\([56,150,58,148]\); the deck map shifts this cycle by two positions. Pole order
is \([204,228]\), consistent with the oriented core.

Both upstairs and downstairs chart coordinates and lattice periods receive
the same meridional sign change. The covering remains
\((u,v)\mapsto(u,2v)\), and the deck remains a half-turn modulo the period.
The source-bound adaptor uses this corrected chart. Its independent gate
checks exact affine partitions, face compatibility, deck action, source
provenance and all relative coarsening steps. Its final pieces have 1,576 N
tetrahedra and 2,334 E tetrahedra. Their branching orders should not be confused
with an ambient orientation assignment in a later four-dimensional assembly.

The cap bridge maps successive vertices of the positive puncture cycle to
increasing \(\theta=i/4\), preserves time, and maps the two ordered poles to
cap hemispheres 0 and 1. Its exact barycentric refinement/collar gate and
native admission pass. The source port has 96 tetrahedra and the cap port has
1,792. This is a concrete collar interface.

With that positive angular coordinate the caps retain
\(\beta_m(\theta,t)=(\theta-mt/2,t)\). The complete integer displacements remain
\([0,4]\) for m=0 and \([-8,4]\) for m=2. The latter is one negative upstairs
turn and two negative downstairs turns. It has not been reduced to parity.

## Trust boundaries and the remaining faithful-X2 gate

The source identification uses the pinned native semantics of the oriented
diagram-to-exterior constructor, finite-vertex removal, and simplification.
The input arrays and final marked isomorphism are exact and independently
checked. A portable elementary-move ancestry for every native intermediate
state is not supplied.

The endpoint map identifies the oriented peripheral slopes, and therefore
extends over meridional filling with the marked knot direction and framing.
An explicit based Wirtinger-generator isomorphism or pointwise identification
of the source puncture with the diagram's chosen cube is not exported. Those
local markings agree up to oriented tubular isotopy; the saved source has its
own explicit puncture and positive normal cycle.

The actual neighborhood is a barycentric neighborhood of the fixed circle.
An equivariant regular-neighborhood theorem, together with the positive
meridian framing, can justify a disk-bundle trivialization up to isotopy. A
finite global equivariant disk-bundle atlas is a stronger deliverable. Without
such an atlas or an explicit accepted theorem boundary, a boundary rotation or
cap formula alone does not certify a pointwise global source twist.

Faithful X2 still requires the source-relative roll/twist identification on the
actual cells, the corresponding mapping-torus and cap assembly, a fixed
orientation anchor, closed native admission, and the requested standardization.
A source control or successful homology calculation does not discharge those
gates. This audit makes no claim that X2 or X0 has been standardized.

## Subsequent finite control assembly check

The separate [control assembly certificate](control_assembly_verification.json)
now checks the actual 103,744-pentachoron control, with f-vector
(4,310, 64,796, 216,104, 259,360, 103,744), directly against the saved input
pieces. Its checker imports no assembler. It regenerates every ordered
product simplex, verifies all cover and orbit-quotient face maps, checks the
full fixed subcomplex including possible isolated lower-dimensional pieces,
locates all ten source seam witnesses in actual pentachora, and validates the
source-times-positive-time orientation anchor. The branch sphere has
f=(106,312,208); the twin has f=(10,24,16) and meets it in exactly two cap
vertices. Current closed native admission receipts are bound to both the
103,744-cell cover and the 51,872-cell quotient. These results establish the
assembled control and its markings; global source-twist identification and
standardization retain their separate status.

```sh
python diagram/verify_control_assembly.py
```

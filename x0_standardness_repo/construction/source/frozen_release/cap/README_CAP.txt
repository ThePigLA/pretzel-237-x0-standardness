EXPLICIT EVEN-TWIST CAP TEMPLATES
================================

Scope
-----
This directory constructs two finite, marked cap models and a common boundary
triangulation, followed by an exact source-bound boundary bridge and finite
local F2/Phi transfer maps.  It advances
the accepted ANALYTIC_TRANSFER_ONLY result to actual cap and spherical boundary
meshes.  The bridge is bound to the actual P(-2,3,7) puncture labels and supplied
diagram-positive meridian convention.  The closed bulk assembly and its X2
admission/standardization are separate parent-package work; this directory
does not construct a closed X0.

The admitted parameters are m=0 and m=2.  No frozen trefoil builder or verifier
is imported, modified, or bypassed.  The accepted integer angular formula is
the input convention:

    beta_m(theta,t) = (theta - m*t/2, t).

Theta is the upstairs angular coordinate of period one.  The deck involution
adds 1/2.  The downstairs angle is 2*theta.  Thus beta_2 has one negative full
upstairs turn and downstairs winding -2; beta_0 has zero winding.  These are
marked attachment distinctions, not a claim of different unmarked manifolds.

Files and actual results
------------------------
data/shared_boundary.json.gz
    One PL S2 x S1 triangulation with f=(432,2224,3584,1792).
    Its equatorial torus has f=(416,1248,832).
    It includes the 256-tetrahedron coarse source boundary and its coordinates.

data/cap_m0.json.gz and data/cap_m2.json.gz
    Each cap has 5,824 pentachora and f=(634,5272,14272,15456,5824).
    Both have exactly the same outer boundary facets and vertex numbering.
    Each has a complete deck permutation, native gluing table, orientation,
    exact rational refinement/collar charts, two fixed disks, and a twin sphere.

    Fixed locus: two disks, combined f=(18,32,16).
    Twin sphere: f=(10,24,16), meeting the fixed locus at its two cone vertices.

data/puncture4_to8.json
    An exact deck-equivariant uniform subdivision from a suspension over a
    four-edge equator to one over an eight-edge equator.  Every old triangle
    is split into two.  This is a refinement, not a quotient or covering map.

verification.json
    Independent stored-data verification, including raw native admission.
    It does not import or call build_caps.py.

native_copy_m0.json, native_copy_m2.json, native_copy_summary.json
    Bounded native cross-checks on unmarked simplification copies.  Copies
    reduced to 14 and 13 pentachora respectively and have H1=0, H2=Z, chi=2,
    valid links as recognized by Regina, and one boundary component.
    simplify() is trusted here.  These are not portable elementary-move or
    standardization certificates; the original cap data are not modified.

data/source_cap_bridge.json.gz, bridge_verification.json
    The completed source-bound bridge has 24,224 pentachora.  Its source and
    cap ports, common refinement, deck action and native admission are checked.

data/finite_F2_Phi_transfer.json.gz, transfer_verification.json
    A finite rational polyhedral certificate for the local puncture-collar
    transfer, using the actual beta2 spherical family.  There are 1,792 Gamma
    tetrahedra and 3,584 Phi four-dimensional cells.  Affine maps, inverses,
    vertex sets, halfspaces, facets, incidence tables and exact seam identities
    are supplied.  The canonical simplicial realization is specified below.

REVIEW_BRANCH_AND_LOCALIZATION.txt, geometric_review.json
    Independent review of the source-global localization proof and of the
    local branch models at all dimensions.  It includes the corrected enlarged
    neighborhood N-plus and an independent saved-matrix check of the six-cell
    simultaneous u/v source-profile transfer.  The remaining normal-chart
    existence theorem dependency is explicit.

GROUP_TAIL_REVIEW.txt, group_tail_review.json
    Independent exact replay of the small residual fundamental-group tail:
    13 retained-relator shortenings and two generator eliminations reach the
    empty presentation.  A four-relator algebraic explanation is also given.
    The full source-to-prefix replay is owned by the parent group_audit/.

The raw 5,824-simplex homology call was stopped because the simplified-copy
calculation answers the same invariant question much more cheaply.  The final
native admission routine deliberately does not calculate raw homology.

Construction and topological meaning
-----------------------------------
Start with an eight-edge angular circle E and a four-edge time circle T.
Let A triangulate E x T by the standard diagonal in each rectangle.  Let B be
its exact shear image (theta,t) -> (theta-t,t).  All images are kept in rational
universal angular/time coordinates before taking the quotient.  Intersect A
and B, retain all collinear edge points, and cone each convex intersection
polygon from its rational barycenter.  This gives a deck-invariant common
refinement R of both torus triangulations.

The sphere is the suspension of E, with fixed north/south poles.  For each
hemisphere and each time interval, triangulate D2 x I by coning its boundary
2-sphere from a fixed interior axis point.  The side annulus has the selected
torus subdivision; the end disks are cones on E.  Successive interval blocks
share those end disks.  The two hemispheres share their equatorial torus.
This supplies an explicit PL product S2 x S1 and compatible refinements.

The spherical extension of beta_m fixes the pole and interval-axis vertices,
and is linear on the listed refined tetrahedra.  For every boundary tetrahedron
the data store its parent coarse-source tetrahedron and a 4-by-4 matrix of
exact barycentric vertex coordinates for beta_m^{-1}.  The matrices give
actual maps on entire simplices, with compatible face restrictions.  Their
local inverses are exact matrix inverses.  On the equator the maps agree with
theta -> theta-m*t/2.  They preserve time and extend over the disk fibers by
the explicit coned interval blocks.  This is a finite representative of the
specified rotating-sphere attachment class.

To fill, begin with H=E x D2, triangulated by staircase products.  Attach an
explicit identity refinement collar from the coarse torus A to beta_m^{-1}(R).
Each original triangle prism is a convex prism.  Its side rectangles are
triangulated by invariant face-center fans; the full prism is coned from an
interior center.  Exact determinants verify its volume partition.  This
preserves H as a PL solid torus, with its original marked core E.

For each sphere hemisphere, the union

    (D2 x S1) union_(E x S1) (E x D2)

is the boundary of D2 x D2, hence a PL 3-sphere.  Cone this explicitly listed
3-sphere to a new fixed vertex.  The result is a PL 4-ball.  The two resulting
4-balls are glued along their common H, producing S2 x D2.  The chosen boundary
marking is beta_m; the outer triangulation is the same R-based spherical
boundary for both cases.  The cap's product/topological identification rests
on this recorded product, refinement-collar and cone construction, rather than
on its homology alone.

The deck half-turn extends over every piece.  It has no simplex inversions.
The fixed sets are the two coned axis circles, hence two disks; the normal
circle of every fixed triangle is checked and has a free deck half-turn.  The
union of the two cones on the original core E is the recorded twin sphere.

The actual winding certificate
------------------------------
The inverse charts recover the source angle-zero time loop as a path of
actual edges in the common equator.  Its universal image has total normalized
displacement (0,1) for m0 and (-1,1) for m2.  In the eight-angle/four-time units
these are [0,4] and [-8,4].  The m2 path uses 16 refined edges; the m0 path uses
8.  Identical endpoint labels do not replace these lifted paths by constants.

The verifier binds the entire spherical barycentric map to the angular map.
Replacing all m2 spherical refinement matrices by the m0 matrices is rejected,
even though each matrix collection separately defines a valid refinement.

Attachment interface for the parent construction
------------------------------------------------
1. The first len(shared_boundary.vertices) vertex ids of both cap files are
   literally the shared boundary ids.  Their boundary facet lists are equal.
   A certified map from the bulk boundary to this shared boundary therefore
   suffices for gluing either cap; no endpoint relabelling is guessed.

2. A shared equatorial vertex is tagged ["E",theta,t] with rational strings
   modulo one.  ["P",p,j] is pole p at time j/4; ["C",p,j] is the fixed axis
   point in hemisphere p and time slab j.  p=0 and p=1 are the two poles.
   The labels ["C",p,j] are auxiliary boundary-mesh vertices, not source
   branch vertices.  The time coordinate at C is (j+1/2)/4.

3. shared_boundary.coarse_source_boundary gives the small 256-tetrahedron
   boundary before the torus overlay.  cap.boundary_to_source_refinement has
   one row per sorted cap.boundary tetrahedron.  row.parent indexes a coarse
   tetrahedron.  row.barycentric lists the images of the four boundary vertices
   in their sorted order, using the parent's sorted vertex order.  This is
   beta_m^{-1}: shared refined target -> coarse cap-source boundary.

4. For a source puncture that is the suspension of a four-cycle, identify its
   ordered equator with e0,e1,e2,e3 and its fixed poles with N,S in
   puncture4_to8.json.  Old e_j is fine f_(2j), and fine f_(2j+1) is the
   midpoint of edge e_j e_(j+1).  Fine f_i has angle i/8.  The bisection map
   commutes with deck e_j -> e_(j+2), f_i -> f_(i+4).  The source certificate
   must supply the ordered cycle and pole correspondence and verify its
   meridian orientation; reflecting the cycle changes the winding convention.

5. A standard staircase product on the puncture sphere times the time circle
   is a different triangulation from the coarse boundary used here.  The
   completed source-bound bridge described below supplies the exact missing
   common-refinement collar.  Matching only boundary vertex sets is insufficient.

6. Other angular cycles are supported by --angles N --times T, provided N is
   even, N>=4, T>=3 and T divides N.  Uniform edge subdivision can arrange this.
   The saved default is N=8,T=4.  General integer twists are not admitted by
   this implementation; its two explicit cases are 0 and 2.

Completed source-bound boundary bridge
-------------------------------------
build_bridge.py and verify_bridge.py construct/check
data/source_cap_bridge.json.gz.  The input is the parent source's
peripheral/puncture_product_interface.json and its cut_source.json.gz.
The interface, its byte hash, the source-cut byte hash and puncture canonical
hash are bound in the output.  A later metadata-only sign-audit correction was
rebound without changing any geometry; bridge_metadata_rebind.json records it.

The actual positive angular cycle is [56,150,58,148], in that order, and the
fixed poles are [204,228].  This reverses the initially extracted normal cycle
[56,148,58,150], as prescribed by the parent diagram-meridian certificate.
The ordered source triangles use min(v,deck[v]) ranks, with poles last.  Their
staircase product with the cyclic four-edge time circle has 96 tetrahedra.

The geometry uses the rational diamond with equatorial vertices
(1,0),(0,1),(-1,0),(0,-1), each hemisphere filled by its disk.  Time is linear
on each interval.  Exact tetrahedron intersections between the source product
and the shared cap boundary give a common polytopal refinement.  Each convex
intersection cell is triangulated by a pulling order on deck orbits.  No cell
contains two distinct members of a deck orbit, so this local order is strict
and equivariant.  The common triangulation has

    f=(1728,10480,17504,8752).

Complete rational barycentric maps to both coarse boundaries are recorded.
Every parent tetrahedron has exactly unit total child volume, and shared-face
maps agree.  All 1,728 common vertex links are checked 2-spheres.

Two exact refinement collars join source -> common <- cap.  For each coarse
simplex prism, pulling from its least bottom vertex cones the top subdivision
and the opposite face prism; induction on faces gives a compatible product
triangulation.  This avoids expanding a generic center-coned prism collar.
The final bridge has 24,224 pentachora:

    source collar:  9,568 pentachora
    cap collar:    14,656 pentachora
    total f=(2184,19152,54248,61504,24224), chi=0.

The bridge is PL S2 x S1 x I.  Its deck-fixed locus consists of two annuli,
with f=(40,96,56).  The independent verifier checks all common-refinement maps,
all 1,888 coarse-prism volumes, boundary ports, deck compatibility, and all 56
fixed-triangle normal circles.  Three integration-specific corruption controls
are rejected.  Raw Regina admission is valid, orientable, connected, has two
boundary components and chi=0.  No large raw homology computation is used.

EXACT PORT API:

  bridge.ports.source.vertices[i]
      is the bridge vertex id corresponding to
      bridge.source_boundary.vertices[i] = [cut_vertex, time_index].
      Map this to the actual bulk product vertex with those source labels.

  bridge.ports.cap.vertices[i]
      is the bridge vertex id corresponding to shared_boundary vertex i,
      also cap_m2 boundary vertex i.  This is a literal index match.

  All bridge vertices outside the ports receive fresh ids.  All cap interior
  vertices receive fresh ids.  Transport the recorded deck permutations with
  these same identifications.  The whole assembly still requires its own
  boundary-incidence, orientation, deck, branch, and native admission checks.

The bridge preserves diagram-positive theta and time.  It is independent of
the cap choice and does not change beta2's unwrapped [-8,4] displacement.
Its existence is a source-bound geometric refinement certificate; it does
not substitute for the parent proof that the supplied cycle is the positive
diagram meridian.

Bridge reproduction:

    python build_bridge.py
    python verify_bridge.py --native

Finite local F2/Phi transfer
---------------------------
The transfer uses exactly the checked beta2 spherical map, not just its
equatorial winding formula.  Cut the time circle of its graph at s=0=1 to
obtain Gamma: a tetrahedral PL map (x,s) -> (R_s(x),s) on S2 x [0,1].  It has
442 vertices and 1,792 tetrahedra.  Each tetrahedron records its affine matrix
and inverse in the appropriate hemisphere's rational diamond coordinates.
The source points are read from cap_m2's full barycentric inverse charts;
the target points are read from the shared spherical boundary.  This binds
the entire family R_s to the actual attachment, including the hemisphere
interiors.  The dependency byte hashes and corrected source interface are
recorded in the transfer file.

R_0 and R_1 are pointwise identity.  Every R_s fixes the two poles and commutes
with the deck half-turn.  Each hemisphere's radius r=|x|+|y| is preserved:
every source and target tetrahedron lies in one radial quadrant, so r is
affine there, and all 1,792 coefficient identities are checked exactly.
The family is a radial-preserving disk isotopy.  No assertion is made that
it is a full circle group action or uniform rotation at each interior radius.

On an added puncture collar S2 x [0,1]_q, with q=0 at the puncture and q=1
at the bulk end, set

    F2(x,q) = (R_(1-q)(x),q),
    Phi(x,q,t) = (R_min(1-q,t)(x),q,t).

Both collar ends are fixed by F2, so it extends by identity outside the collar.
Phi is identity at t=0 and at q=1; its q=0 restriction is the actual R_t cap
loop; its t=1 restriction is F2.  These assertions hold on entire cells,
not merely at sampled parameter values.

Here is the finite construction.  Split the parameter square at q+t=1.
For each Gamma tetrahedron K in (x,y,s) coordinates, form two convex cells:

    A_K = {(x,y,q,t): (x,y,t) in K, 0 <= q <= 1-t},
    B_K = {(x,y,q,t): (x,y,1-q) in K, 1-q <= t <= 1}.

They are affine extrusions with at most eight vertices.  A_K uses s=t;
B_K uses s=1-q.  Substitute these affine expressions into the exact Gamma
matrix.  The resulting 5-by-5 homogeneous affine matrix preserves q,t and
has its recorded exact inverse.  The target cell is the same extrusion of
the image Gamma tetrahedron.  On the common branch face q+t=1 the two
matrices differ by (q+t-1) times the spatial s-coefficient of Gamma, hence
agree identically.  The verifier checks this coefficient equation for every
Gamma tetrahedron, as well as every final-time F2 restriction.

The output has 1,306 vertices, 3,584 four-dimensional cells, and 12,544
three-dimensional facets.  Every interior facet has two incident cells.
Boundary facets are 1,792 at q=0, 1,792 at t=1, 16 at t=0, and 16 at q=1.
Facet completeness is checked independently from the four extruded Gamma
faces and the two extrusion ends.  The stored halfspaces and source/target
face posets are checked separately.  Coverage follows from the checked
Gamma triangulation and the two exhaustive parameter-square branches;
face agreement follows from the affine maps' shared vertex values.
Consequently these data define a PL homeomorphism of S2 x I x I, preserving
the last two coordinates, with the stated boundary restrictions.

To materialize a simplicial map, take the barycenter of every nonempty face
of each convex cell, identifying shared faces by their recorded vertex sets,
and take one simplex for every strictly increasing maximal face chain.
Faces of lower dimension are obtained by intersecting the listed facets.
An affine map takes the vertex-average barycenter of a face to that of its
image; the source and target chains therefore define compatible simplices.
This construction is canonical and deck-equivariant.  The larger expanded
4-simplex lists are deliberately not materialized in this local certificate.
The finite cell lists and exact inverse affine maps already determine them.

The integer lift is retained explicitly.  On the equator,

    Gamma expression for F2: theta + q - 1,
    accepted F2 lift:         theta + q,
    Phi branch A:             theta - t,
    Phi branch B:             theta + q - 1.

The two F2 lifts differ by the constant lattice translation +1, and the
recorded bulk D2 lift contributes +1.  Thus at t=1 the full lifted seam is

    1 + (theta + q - 1) = theta + q.

The cap boundary still has upstairs winding -1, downstairs winding -2 and
displacement [-8,4].  Dropping the full bulk turn is rejected even though
the physical endpoint map is unchanged.  The other transfer corruption
controls change an affine coefficient and omit a four-dimensional cell;
both are rejected.

SCOPE OF THE GLOBAL IDENTIFICATION: this file gives exact finite PL maps on
the added source-bound puncture collar, and identity extension to its
complement.  Identifying that relative F2 class with the original
torus-supported meridional twist uses the oriented normal-disk-bundle
localization argument in the parent package.  The reviewed argument first
applies the exact mixed source-profile transfer and enlarges N by its saved
torus adaptor: the actual puncture touches the original cut torus, but lies
strictly inside this N-plus.  The outer grid and puncture can then both be
prescribed in a global framed chart.  The actual regular neighborhood and
positive normal meridian have been certified there, but this directory does
not export a global normal-disk-bundle fiber atlas on N-plus.  In particular,
the local finite certificate is not described as a fully meshed source-global
torus-to-puncture isotopy.  The separately certified roll support lies away
from this added puncture collar.

Transfer reproduction:

    python build_transfer.py
    python verify_transfer.py

Reproduce
---------
Python 3's standard library is sufficient for generation and offline checks:

    python build_caps.py --output data
    python verify_caps.py --data data --output verification.json

The four-to-eight puncture certificate is a separate static input.  Generation
does not overwrite it.  With Regina 7.4.1 installed:

    python verify_caps.py --native --output verification.json
    python native_copy_homology.py --seconds 45

The six corruption controls cover parity-only winding, parity-only spherical
map, wrong twist binding, corrupted deck permutation, missing cap facet, and
corrupted rational boundary coordinates.  Native admission checks the original
facets' validity, orientability, connectedness, one boundary component and
Euler characteristic.  The final receipts carry the exact executed results.

# Explicit equivariant Seifert/Montesinos source

## Construction

The producer builds geometry from the signed pretzel parameters `(-2,3,7)`.
It does not ask Regina for a named Seifert manifold.

Let `P` be the double of a hexagon along the three alternating edges `01`,
`23`, and `45`. The other three edges are separately subdivided in the two
halves. This gives a pair of pants with 14 vertices and 18 triangles. The
reflection exchanges the two hexagons and fixes their three shared edges.
Each oriented boundary circle has four edges. The `+` hexagon has its displayed
orientation and the `-` hexagon has the opposite orientation.

Cross `P` with a four-edge oriented circle. The involution is reflection on
both factors. On each boundary torus it is `(x,y) -> (-x,-y)`, where `x`
follows the oriented pants boundary and `y` is the fibre coordinate.

For each signed parameter `p`, use a four-triangle oriented disk times a
four-edge circle. Its involution is conjugation on the disk and reflection
of the circle. Attach its boundary using

    (x,y) = A_p (u,v),   A_p = [[p,1],[1,0]],

where `u` is the positive meridian angle and `v` is the solid-torus longitude.
The determinant is -1, so these maps reverse the induced boundary orientation.
They commute with the boundary involution because they are integral linear
maps commuting with `-I`.

The columns are the explicitly retained winding vectors `(p,1)` and `(1,0)`.
The meridional slopes are thus `-2 q_0+h`, `3 q_1+h`, `7 q_2+h`. In positive
denominator convention the Seifert pairs are `(2,-1),(3,1),(7,1)`.

The matrix decomposition `A_p = [[1,p],[0,1]] [[0,1],[1,0]]` records the
individual signed shear and coordinate exchange. In the corresponding
Montesinos convention these are rational tangles of slopes `1/p`.

For background on the rational-tangle/Seifert-cover construction and its
conventions, see Jeffrey Meier, *Small Seifert fibered surgery on hyperbolic
pretzel knots*, Sections 2.3--2.4, especially Proposition 2.3:
https://msp.org/agt/2014/14-1/agt-v14-n1-p16-p.pdf

The external oriented-diagram comparison is a separate gate. A matching
unmarked Seifert manifold or homology result alone does not settle branch
knot handedness or its peripheral signs.

## Finite triangulation and certificate data

The common refinement of each boundary pair is computed by exact rational
polygon intersection in the universal covering plane, retaining lattice
translations. Every polygon is subdivided using its edge midpoints and
interior centroid. This rule is invariant under the involution. The parent
product faces receive the corresponding edge subdivisions; interior faces
are coned from their centers. Each triangular prism is coned from its own
center to this triangulated boundary.

The final source has 1,256 vertices and 7,440 tetrahedra. There are 240
common-refinement polygons and one connected ramification circle with 18
vertices. The quotient has 3,720 tetrahedra. Each source tetrahedron has four
distinct deck-vertex orbits. Vertex ranks `min(v,deck[v])` are therefore
deck-invariant and distinct within each tetrahedron.

`data_oriented/source.json.gz` contains the most complete master. In addition
to facets, deck involution and source labels, it records the surface blocks,
reflections, boundary cycles, exact common polygons, explicit paths for the
two filling-basis columns, and all tetrahedral orientations.

`data/master_orientation.json` is an orientation supplement to the earlier
master with the same facet hash. These signs are computed from the oriented
product identity

    d(triangle x I) = d(triangle) x I
                     + triangle x {1} - triangle x {0},

the exact boundary-polygon orientations, and positive coning. They are not
chosen by an arbitrary BFS seed. The producer verifies every gluing parity
equation and every orientation-preserving deck equation. Each sign is relative
to the increasing vertex order of the corresponding listed tetrahedron.

`data/gluing_paths.json` retains closed refined-edge paths with their universal
coordinate lifts. These represent the rational-tangle filling meridian and a
chosen filling longitude, **not** the meridian and zero-longitude of the full
pretzel knot. The latter must be reconstructed on a regular-neighborhood
boundary of the ramification circle.

## Native cross-check

Regina 7.4.1 admits both exported raw triangulations as valid, closed,
connected and orientable. On unmarked copies, native simplification gives:

- Quotient: one tetrahedron, `bkaagj`, `isSphere() = true`, and `H1 = 0`.
- Cover: nine tetrahedra, `jLAMLLQbcbdeghhiixxjxxjqisj`, and `H1 = 0`.

The source-derived abelian relations have matrix

    [[ 1, 1, 1, 0],
     [-2, 0, 0, 1],
     [ 0, 3, 0, 1],
     [ 0, 0, 7, 1]],

of determinant 1, agreeing with the native homology. This is a cross-check,
not a substitute for the marked construction.

The native computations use Regina's trusted simplification and recognition
algorithms on copies. They do not provide a new portable elementary-move
certificate and do not transport the marks to the simplified copies.

## Scope

This stage constructs and admits a three-dimensional marked branched-cover
pair, subject to the stated source-interpretation and external sign-comparison
gate. It constructs no roll mapping torus, cap attachment in dimension four,
`X2`, `X0`, or four-dimensional standardization certificate.

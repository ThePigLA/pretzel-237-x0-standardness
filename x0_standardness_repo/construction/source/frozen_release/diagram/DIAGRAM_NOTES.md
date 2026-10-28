# Explicit marked diagram for J = P(-2,3,7)

This directory supplies a concrete oriented diagram, a rational PL knot embedding,
a marked trivial-arc puncture, and diagrammatic peripheral words. It does not by
itself construct a triangulated branched cover or transport these marks to one.

## Source and handedness

The source is the numerator closure of the left-to-right sum

    T(-1/2) + T(1/3) + T(1/7).

The official SnapPy/Spherogram tutorial explicitly calls this construction the
(-2,3,7) pretzel:

https://snappy.computop.org/spherogram.html

The generator uses that rational-tangle expression, without a knot-name or census
lookup. A separate direct construction stacks 2, 3, and 7 crossings in three
rectangles, connects corresponding top and bottom endpoints between rectangles,
and adds the two outer numerator-closure arcs. The recorded rotations by 0 or
180 degrees give a literal crossing-by-crossing identification of these two
constructions. The 12-crossing PD code is retained in full.

This convention has twelve crossings of Spherogram sign -1, hence writhe -12.
It is frozen as supplied; the generator does not silently mirror it. A paper using
a different sign convention for pretzels must be matched by comparing the actual
signed diagram, not by matching its determinant or unmarked covering manifold.

The raw crossing slots are numbered counterclockwise, with 0--2 under and 1--3
over. For a positive box the physical slots (SW,SE,NE,NW) are (0,1,2,3). For the
negative box they are (3,0,1,2). This fixes the sign convention independently of
textual names.

## Orientation and geometric mark

The outside top arc runs from right to left. Its middle point is

    p = (6,2,0), tangent = (-1,0,0).

The marked arc is PD strand 22, entering the upper-left port of the -2 box. The
positive normal basis is (e_y,-e_z): tangent followed by this ordered pair is
positively oriented in the ambient (x,y,z) coordinates.

The puncture is the closed PL cube with centre p and half-width 1/4. The full knot
polygon has 89 distinct rational vertices. An exact rational segment check shows
that this cube meets the knot only in the straight arc with oriented endpoints

    (25/4,2,0) --> (23/4,2,0).

Consequently it is a specified ball/trivial-arc pair. This is geometric diagram
data; locating its image in a subsequently constructed triangulation remains a
separate transport obligation.

The positive meridian is the small diamond in the normal plane recorded in JSON.
It winds positively in the ordered (e_y,-e_z) basis and has linking number +1.
Its basepoint is (6,17/8,0). The diagram group uses the usual stems into the upper
half-space from the source diagram and labels the marked arc's meridian x1.

## Preferred longitude and finite group certificate

Number the under-to-under Wirtinger arcs successively from the marked point.
At an oriented crossing of sign epsilon, use the relation

    over^(-epsilon) * incoming * over^(epsilon) * outgoing^(-1) = 1.

This is the relation convention implemented in Spherogram's knot_group method;
the verifier independently derives each incoming, outgoing, and over arc from
the explicit oriented wiring.

The blackboard push-off word, in this based convention, is

    W = x8^-1 x10^-1 x11^-1 x12^-1 x8^-1 x9^-1
        x2^-1 x6^-1 x2^-1 x3^-1 x4^-1 x5^-1.

It follows the chosen orientation once and has linking number equal to the
writhe, -12. Therefore the preferred zero longitude is

    mu = x1,
    lambda = W x1^12.

The positive meridian correction is part of the source record. It must not be
dropped just because the abelianization of another presentation is known.

There is an explicit certificate that mu commutes with lambda. Let b_i be the
signed over-meridian at the i-th underpass, let g_i be its incoming meridian, and
let R_i = b_i^-1 g_i b_i g_(i+1)^-1. With W_i = b_0...b_(i-1), the product

    product_i (W_(i+1) R_i W_(i+1)^-1)

freely reduces to mu W mu^-1 W^-1, which also freely equals [mu,lambda]. The JSON
contains every conjugator and relator; the checker performs the literal free
reductions. Commutation and zero abelianization are algebraic checks. The actual
peripheral interpretation additionally uses the explicit diagrammatic push-off
and framing prescription given above.

## Determinant and Seifert-cover interface

The independently rebuilt Fox-colouring matrix has a reduced determinant -1,
so the knot determinant is 1. A Spherogram-produced 10-by-10 Seifert matrix gives
the additional polynomial cross-check

    Delta(t) = t^10 - t^9 + t^7 - t^6 + t^5 - t^4 + t^3 - t + 1.

The standard-library verifier checks the matrix polynomial by eleven exact
integer evaluations: both sides have degree at most ten. In particular,
Delta(-1)=-1. These are cross-checks of the source diagram, not identification by
invariants.

For the Montesinos double-cover construction, the intended convention is an
oriented pair of pants times an oriented fibre h, with base boundary classes
c1+c2+c3=0, and the three filling meridians

    2 c1 - h,    3 c2 + h,    7 c3 + h.

Thus the filling pairs (alpha,beta) are (2,-1), (3,1), (7,1); the corresponding
slopes beta/alpha are exactly the source slopes. The group convention is h
central, ci^alpha_i h^beta_i=1, and c1 c2 c3=1. Its abelianization matrix, with
columns (c1,c2,c3,h), is

    [[2,0,0,-1], [0,3,0,1], [0,0,7,1], [1,1,1,0]],

whose determinant is +1. The rational sum is -1/42. In the normalized convention
q1 q2 q3=h^b this becomes b=-1 and pairs (2,1),(3,1),(7,1).

The Montesinos/Seifert covering relation and convention comparison are discussed
in Bonahon and Siebenmann, Geometric Splittings of Classical Knots and the
Algebraic Knots of Conway, section 11.2 (printed page 133):

https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

This directory specifies the filling interface exactly. It does not supply the
cellwise quotient map identifying a constructed Seifert triangulation with the
signed source diagram. That map and its peripheral/puncture transport remain
required; a named Seifert manifold with H1=0 would not discharge them.

## Run

The saved data can be checked with Python's standard library alone:

    python verify_diagram.py --receipt diagram_verification.json

Regeneration requires the pinned environment with Spherogram 2.4.1 and SymPy
1.14.0:

    python generate_diagram.py
    python verify_diagram.py --with-library --receipt diagram_verification.json

The checker does not import the generator. It checks the rotation system has
14 faces and Euler characteristic 2; the source is one component; the PD and
direct wiring agree; the 89-vertex spatial polygon is embedded; there are exactly
the intended twelve projected crossings; the cube puncture is clean; the
meridian orientation and longitude correction are present; and the explicit
commutator certificate and determinant calculations pass. Four modifications of
this new source (wrong slope, missing longitude correction, wrong puncture size,
and altered over/under chart) must be rejected. These do not rerun any accepted
trefoil control.

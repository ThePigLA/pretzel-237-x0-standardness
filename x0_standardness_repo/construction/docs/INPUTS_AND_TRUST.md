# Construction inputs and verification scope

## What is specified mathematically

The example is fixed by the signed pretzel parameters `[-2, 3, 7]`, the
associated rational tangles, the oriented peripheral conventions, one
longitudinal roll, and the zero-twist cap. The builders use exact combinatorial
and rational-coordinate descriptions of the source, cut, collar, and cap.
The four-dimensional product grid is `[3, 4, 4]` in the accepted master.

The original oriented source has 7,440 tetrahedra. The fixed source
compression produces 628 tetrahedra. The master has the following parts:

| Part | Pentachora |
| --- | ---: |
| Exterior product `E` | 37,344 |
| Neighborhood product `N` | 25,216 |
| Longitudinal roll collar `central` | 11,136 |
| Source-to-cap bridge | 24,224 |
| Zero-twist cap | 5,824 |
| **Total** | **103,744** |

## Supplied construction choices

The recipe includes explicit finite choices made during development:

- `coarsened/coarsening_certificate.json.gz` records the source compression.
  The source endpoint is reconstructed by replaying all local replacements.
- `adaptor/coarsening_certificate.json.gz` records the compact adaptor mesh.
  The adaptor is first generated and then its compression is replayed.
- The raw peripheral chart and marked native transport fix the intended
  meridian, zero longitude, signs, and orientation conventions. The default
  build uses this fixed transport; the native audit repeats the source-bound
  native calculation and checks its marked endpoint.
- `peripheral/puncture_product_interface.json` fixes the small spherical
  product interface and its time subdivision. The geometric bridge checker
  independently binds its cells, orientations, and ports to the source and cap.

These inputs specify a reproducible triangulation. They are not a supplied
answer to the final four-dimensional assembly: the assembler generates its
facets, deck action, markings, and cover/quotient gluing tables from the pieces.

## What the fresh checks establish

The construction audit has separate checks for the signed diagram, source
compression, source orientation and clipping, exact marked exterior endpoint,
torus adaptor, roll collar, cap, bridge, source/interface consistency, zero-cap
identity, and independent marked assembly. Additional original-source/control
geometry checks are retained and identified in the stage plan.

The source-compression replay is a separate stage. Orientation transport alone
does not prove that a compression preserves PL type. Its replay checks the
relative three-ball replacements and protected fixed locus. Likewise, the
adaptor's compression certificate is replayed, not inferred from its endpoint.

The final assembly audit checks the complete fixed subcomplex, the covering
face squares, the deck and orientation data, both marked surfaces, the
attachment maps, and all finite cells. The optional link replay establishes
that every vertex link of the cover master is a PL three-sphere.

The reference cover and quotient native-admission receipts remain explicitly
inherited evidence in a standard-library audit. The native audit regenerates
those admissions. A receipt is used only when its input hash matches the
actual gluing table being checked.

## Native semantics and geometric arguments

The source identification retains the documented semantics of the explicit
Spherogram diagram construction and SnapPy/Regina import, finite-vertex
removal, simplification, and peripheral transport. A fresh native rerun and
an exact marked endpoint comparison check this computation, but they do not
expand its internal transformations into a supplied elementary move trace.

The regular-neighborhood, collar, mapping-torus, and cap identifications also
use the mathematical arguments in the construction dossier and paper. The
finite checkers validate their stated meshes, maps, signs, winding numbers,
and attachment hypotheses. Native recognition or cell counts alone do not
replace those identifications.

No part of this module is presented as proof-assistant verification. The
construction audit reports the operations it freshly performed and the
source semantics it retains.

## Historical bindings

The unchanged assembler requires the old source-bound X2 acceptance receipt
and the hash of the original marked control. Those are retained as operational
provenance. The X0 cells use the direct zero-twist cap; the control's
standardization is not a mathematical premise of that construction.

Some original receipts describe the project's then-current status. They are
preserved byte-for-byte because subsequent receipts bind their hashes. The
top-level documentation describes this module's current scope. The later
standardization certificates are outside this construction directory.

## Fresh outputs versus original evidence

The generation and audit drivers use fresh working trees. They do not rewrite
the imported source, reference master, or historical evidence. The original
and fresh receipts have different roles and are kept separate.

The bridge producer stores its run time inside a geometry object. To preserve
the fixed recipe's downstream hashes, the build explicitly retains the
historical `elapsed_seconds` field after comparing all freshly generated
geometric fields. Its own measured time is recorded in the build report.
This is a serialization/provenance convention, not a geometric assumption.

The source inventory records every imported file's original path, byte count,
and SHA-256. New driver code and validation records are covered by the new
release manifest. The outer original package manifest is retained for provenance;
the new module omits the later X0 reduction and terminal-standardization
packages. The frozen source tree retains historical control material needed
to preserve its original bindings.


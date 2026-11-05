# Adding this module to the existing repository

The ZIP contains one directory, `motoo_z/`. Place it inside the existing
`x0_standardness_repo/` directory, alongside `construction/`, `reduction/`, and
`paper/`.

This module contains the new Motoo construction and proof data. It is designed
as an addition to the repository; it does not replace an older or newer paper
draft. The reference X0 construction/reduction modules remain their own proof
route.

Add a row to the root README's directory table:

```markdown
| [motoo_z/](motoo_z/README.md) | Katie construction, simplification searches, and complete certificates for the original-drawing Motoo route | `motoo_z/build.py`, `motoo_z/search.py`, `motoo_z/verify.py` |
```

Suggested accompanying text:

> The Motoo module reconstructs the 24,170-pentachoron closed triangulation
> from the marked diagram via Katie and explicit closing data. It includes
> simplification search programs and the complete recorded 781,277-move
> connection from the closed triangulation to the explicit standard sphere.
> Construction, terminal-only verification, and full verification have separate
> result labels. See the module README for its pinned environment and software
> trust boundaries.

The module has its own `MANIFEST_SHA256.json`. Its largest file is the original
32,810,768-byte long-move archive under `certificates/`; no move chunks are
omitted. Fresh verification outputs and the complete 19-link recognition
proofs from the preparation run are included under `verification/`.

If the main repository's manifest is expanded to include the module, generate
that updated root manifest as part of the actual repository integration. The
module manifest is already complete for the files in this ZIP.

The original confidentiality and source notices remain in the imported
material. This package preparation did not upload or publish any changes to
GitHub and did not assign new licensing or citation metadata.

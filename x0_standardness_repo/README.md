# X0 standardness: construction, reduction, and proof

This repository accompanies the proof that the double branched cover

$$
X_0=\Sigma_2\bigl(S^4,\tau^0\rho^1\,P(-2,3,7)\bigr)
$$

is PL homeomorphic to the standard four-sphere, and hence has the standard
compatible smoothing. It contains the marked triangulation construction,
the complete recorded reduction, the mathematical checkers, and the paper.
The standardization is nonequivariant: the reduction does not transport the
deck involution or marked surfaces through the move sequence.

**Status:** draft prepared for review. Authorship, affiliations,
acknowledgements, licensing, repository citation metadata, and DOI decisions
remain open. The manuscript's scholarly references are included; publication
metadata will be added after those decisions are made.

## Start here

| Directory | Role | Main entry point |
| --- | --- | --- |
| [paper/](paper/README.md) | Mathematical statement, proof, and computational dependencies | [Manuscript PDF](paper/x0_standardness_proof_note.pdf) and [LaTeX source](paper/x0_standardness_proof_note.tex) |
| [construction/](construction/README.md) | Generate the marked complex and audit its geometric identification | `construction/build.py`, `construction/audit.py` |
| [reduction/](reduction/README.md) | Replay the reduction of that exact complex to the standard sphere | `reduction/verify.py` |

The construction produces a closed marked simplicial complex with **103,744
pentachora and 4,310 vertices**. The reduction checks its vertex links,
contracts and simplifies it, punctures and caps through explicitly checked
interfaces, and connects the resulting small triangulation to the explicit
double of a four-simplex. The latter is the two-pentachoron standard sphere
with signature `cPkbbbbaaaaaaaa`.

Reproduction follows fixed construction choices and recorded certificates.
It checks those choices and transformations afresh; it does not repeat the
searches that originally found them. Each module contains its source code,
input identities, logs, and validation records.

## Reproduce the construction

Run commands from this directory. Use Python 3.12, and choose a new output
directory for every run.

```sh
python3 -E -B construction/build.py --output runs/x0-build
```

The default geometry build regenerates the source, replays its compression,
reconstructs the peripheral adaptor, roll collar, caps and bridge, and
assembles the marked master. It compares all fields of the output objects
with the accepted references. The build uses Python's standard library.

The generated master and its native gluing table are in
`runs/x0-build/candidate/`. Audit the finite construction checks with:

```sh
python3 -E -B construction/audit.py \
  --candidate runs/x0-build/candidate \
  --output runs/x0-construction-audit
```

For the complete construction audit, including fresh native source transport,
native manifold admissions, and all vertex-link certificates, install the
construction dependencies in a separate environment:

```sh
python3.12 -m venv .venv-construction
.venv-construction/bin/python -m pip install \
  -r construction/environment/requirements.txt
.venv-construction/bin/python -E -B construction/audit.py \
  --candidate runs/x0-build/candidate \
  --native --links \
  --output runs/x0-construction-audit-full
```

The [construction README](construction/README.md) distinguishes fresh checks,
fixed recipe inputs, and inherited evidence when native checks are omitted.

## Replay the complete reduction

Full reduction requires **CPython 3.12, Linux x86-64, and Regina 7.4.1 with
the recorded native engine**. The driver checks that environment before
starting the macro segment. On Linux, install the exact wheel:

```sh
python3.12 -m venv .venv-reduction
.venv-reduction/bin/python -m pip install \
  --require-hashes --only-binary=:all: --no-deps \
  -r reduction/environment/requirements.txt
.venv-reduction/bin/python -E -B reduction/verify.py \
  --candidate runs/x0-build/candidate \
  --output runs/x0-reduction
```

On a Mac, use the digest-pinned Linux container described in the
[runtime instructions](reduction/environment/README.md). Those instructions
cover Intel and Apple Silicon, retain outputs on the host, and run with the
host user's numeric IDs. The matching native Linux engine has passed the
complete replay; the container recipe itself has not been built in the
preparation environment, where Docker was unavailable.

The candidate's master and gluing-table hashes must equal the released inputs.
Omitting `--candidate` replays the exact reference master included in
`reduction/`. The driver records which input was selected and checks every
stage's fresh output against the next stage's input.

A complete successful run reports:

```text
PASS_MASTER_PL_HOMEOMORPHIC_TO_STANDARD_S4
```

Its main receipt is `runs/x0-reduction/REDUCTION_REPLAY.json`. The directory
also contains the replayed intermediate endpoints, capping maps, independent
terminal receipts, logs, and `terminal/standard_sphere.json`.

### Smaller checks

The elementary reduction uses only the standard library on a POSIX system:

```sh
python3 -E -B reduction/verify.py --scope elementary \
  --candidate runs/x0-build/candidate --output runs/x0-elementary
```

The small terminal paths can be checked independently on Mac or Linux with
Regina 7.4.1:

```sh
python -E -B reduction/verify.py --scope terminal \
  --output runs/x0-terminal
```

These partial scopes have distinct result labels. Neither reports that the
full master-to-sphere chain was replayed. To check the released input bytes
without mathematical replay:

```sh
python3 -E -B construction/common.py
python3 -E -B reduction/verify.py --integrity-only
```

**Keep assertions enabled for mathematical replay.** Do not use `-O` or `-OO`
for the construction and reduction entry points: some unchanged elementary
checkers use assertions for mathematical preconditions. Those drivers reject
optimized execution and refuse to overwrite existing run directories.

## Verification scope

The geometric identification belongs to `construction/`; the PL-type
certificate belongs to `reduction/`. The elementary reduction checks all
4,310 vertex links, 3,308 prefix contractions, and 629,874 subsequent
contractions or bistellar flips. The terminal certificate contains 8,436
individually replayed moves, with a separate implementation checking all
594 admissible 2–0 operations.

The native punctured stage contains 5,173 `simplify()` blocks. Every recorded
block output is checked exactly, but its internal simplify moves are not
expanded into individual records. Regina's move implementations remain a
declared dependency. The terminal independent 2–0 checker reconstructs the
facet pairings without invoking Regina's 2–0 functions, while still sharing
its signature codec, skeleton queries, and constructors. The whole proof is
not a pure Pachner transcript or a proof-assistant formalization.

See [construction assumptions](construction/docs/INPUTS_AND_TRUST.md),
[reduction assumptions](reduction/docs/CHAIN_AND_TRUST.md), and the paper
for the precise mathematical and computational boundaries.

## Recorded validation

Fresh relocated checks completed on 28 October 2026 with Python 3.12.14 on
Linux x86-64:

| Validation | Measured time | Record |
| --- | ---: | --- |
| Construction geometry build | 139.6 seconds | [Construction validation](construction/verification/release_validation.json) |
| Complete 21-stage construction audit | 125.6 seconds | [Construction validation](construction/verification/release_validation.json) |
| Fresh construction rebuild for reduction | 128.7 seconds | [Reduction validation](reduction/verification/release_validation.json) |
| Complete six-stage reduction replay | 306.7 seconds | [Reduction validation](reduction/verification/release_validation.json) |
| Two-pass paper build and reference-text check | 1.3 seconds | [Paper validation](paper/verification/release_validation.json) |

Both construction builds reproduced the accepted outputs byte-for-byte.
The records bind exact code, inputs, and fresh results, and include targeted
rejection checks for optimized execution and altered inputs. These timings
describe the recorded machines and runs; performance elsewhere will vary.
The rebuilt paper has 13 pages, no typesetting warnings, and identical text
and rendered page pixels to the preserved Revision 4 reading copy. This
paper-build check does not rerun the mathematical certificates.

## Manuscript and release preparation

The [paper directory](paper/README.md) contains the current manuscript,
editable source, bibliography, and build instructions. Existing mathematical
claims and proof dependencies are preserved while the repository is assembled.

With Python 3 and pdfLaTeX installed, compile the paper into a new directory:

```sh
python3 -E -B paper/build.py --output runs/paper-build
```

The paper source embeds its bibliography, so this build does not need BibTeX
or the topology software. With Poppler installed, add `--verify-reference`
to require matching page counts and normalized text against the supplied
reading copy. The [repository map](paper/REPOSITORY_MAP.md) connects
Appendix A's original C1-C4 archive designations to the current module paths.

The current snapshot is prepared for review. It does not assign a new license
on behalf of contributors, and it contains no provisional DOI or repository
citation entry. Existing source notices and exact provenance are retained.

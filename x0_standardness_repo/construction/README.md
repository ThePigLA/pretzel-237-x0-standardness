# Construction of the marked X0 triangulation

This is the first repository component for the X0 proof. It contains readable
generation code, the fixed construction choices, independent geometric
checkers, and the accepted marked reference triangulation.

The fixed example is the double branched cover of the zero-twist, one-roll
spin of the signed pretzel knot P(-2,3,7). The output is the closed marked
four-dimensional simplicial complex with **103,744 pentachora and 4,310
vertices**. Its deck involution, orientation, ramification sphere, twin
sphere, and attachment maps are part of the output.

The construction and its geometric identification are the scope of this
directory. The subsequent reduction to the standard sphere belongs to the
next repository component.

## Build the triangulation

Use Python 3.12 with assertions enabled. From the repository root:

```sh
python3 -E -B construction/build.py --output runs/construction-build
```

The default `--level geometry` performs the following fresh computation:

1. Generate the oriented Montesinos source from the explicit parameters.
2. Replay the recorded source compression, checking its local moves.
3. Reconstruct the source cut and transport its orientation.
4. Check the fixed peripheral chart and apply the recorded marked transport.
5. Reconstruct the torus adaptor and replay its recorded compression.
6. Generate the longitudinal roll collar, the caps, and the source-cap bridge.
7. Assemble the marked X0 master and its cover and quotient gluing tables.
8. Compare every field of the three output objects with the accepted references.

The build uses the Python standard library. It makes a fresh working copy of
the source tree and refuses to overwrite an existing output directory. It
leaves the supplied source and reference data unchanged.

The principal outputs are:

| File inside the output directory | Meaning |
| --- | --- |
| `candidate/marked_candidate.json.gz` | Facets, labelled vertices, deck map, orientation, marked surfaces, parts, and attachment maps |
| `candidate/cover_regina.json.gz` | Literal facet-pairing table of the 103,744-pentachoron cover |
| `candidate/quotient_regina.json.gz` | Literal orbit-cell quotient, with 51,872 pentachora and covering maps |
| `candidate/assembly_receipt.json` | Fresh assembler receipt |
| `components/` | Intermediate data; generated or retained origin is recorded in the build report |
| `build_report.json` | Stage results, exact comparisons, input identities, and timing |
| `logs/` | Complete stage output |

For a shorter assembly-only run from the supplied geometric components:

```sh
python3 -E -B construction/build.py --level assembly --output runs/assembly-only
```

That option explicitly uses the recorded intermediate meshes. Use the default
geometry level for the construction reproduction.

## Audit the construction

The independent audit recomputes geometric checks in separate working copies.
It also checks the source compression, exact marked exterior endpoint,
clipping and orientations, adaptor compression, roll, cap and bridge,
attachment data, and independent X0 assembly.

```sh
python3 -E -B construction/audit.py --output runs/construction-audit
```

To audit the freshly generated outputs rather than the supplied reference:

```sh
python3 -E -B construction/audit.py \
  --candidate runs/construction-build/candidate \
  --output runs/generated-construction-audit
```

The optional `--links` check replays the certificate for every one of the
4,310 vertex links. The optional `--native` checks additionally rerun the
library-based signed-diagram construction, source manifold admission,
marked exterior transport, and cover/quotient manifold admission.

Install the native dependencies in an isolated environment:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r construction/environment/requirements.txt
python -E -B construction/audit.py \
  --candidate runs/construction-build/candidate \
  --native --links --output runs/construction-audit-full
```

The five primary dependencies are pinned. `environment/requirements-linux-tested.txt`
also records the complete package versions used for this release's Linux
validation. `environment/tested_environment.json` records the interpreter,
platform, and hashes of the installed native topology modules. These are an
environment record, not an archived container image.

Each audit stage retains its log and fresh receipt. A nonzero exit code or a
failed required result makes the audit fail. An old successful receipt is not
accepted as a fresh stage result. `--plan` lists the stages without executing
them. Do not use `python -O`: several inherited checkers use assertions as
mathematical precondition checks.

## Fixed recipes and exact reproduction

The source and adaptor compression were originally found by time-limited
searches. Their finite recorded choices are now recipe inputs; reproduction
replays those choices and checks their legality. It does not repeat the
searches. The marked native transport and the small peripheral/interface
records are also explicit supplied inputs. Their geometric and source
bindings are audited separately. See [docs/INPUTS_AND_TRUST.md](docs/INPUTS_AND_TRUST.md).

One legacy bridge model contains an elapsed-time field inside data that the
master hashes. The builder regenerates the bridge geometry, checks every
other field, and retains that field as **historical recipe metadata**. The
fresh elapsed time is reported separately. JSON formatting and the gzip
platform byte are reproduced from the recorded formats; no reference master
is copied into the generated output. Both parsed-object equality and compressed
file equality are reported. Different compression libraries may produce
different bytes even when the mathematical objects are identical.

## Directory guide

| Path | Contents |
| --- | --- |
| `build.py` | Coordinated generation from the fixed construction recipe |
| `audit.py`, `audit_support.py` | Fresh independent construction audit |
| `common.py` | Explicit integrity checks for the imported inputs |
| `source/frozen_release/` | Original source generators, local-move checkers, exact geometric data, and native transport evidence |
| `source/scripts/` | Marked X0 assembler, independent assembly checker, native admission, and historical admission binding |
| `source/review/` | Independent zero-cap geometry checker and its recorded evidence |
| `source/portable_links/` | Optional complete portable vertex-link replay |
| `reference/x0/` | Accepted output objects and their original receipts |
| `provenance/` | Original package manifest and exact mapping of imported files |
| `environment/` | Dependency specifications and the tested native environment |
| `verification/` | Fresh results and release-validation record |

`source/frozen_release/` preserves the original relative path layout because
the independently checked scripts and their receipts bind those exact files.
The source is fully expanded and readable. This is not a wrapper around an
opaque binary or a missing private archive.

## Exact reference identity

Marked master file SHA-256:

```text
251d7ad38235869fa268e9ee6ec9da2f397cf1ca71db267d0bf0ca9b434c5e4c
```

Canonical facet-list SHA-256:

```text
f40b43c09373dff94cde542a3076a8c2f482b1632726b63eed2c5610ddd1e92b
```

To check the supplied inputs without doing a mathematical replay:

```sh
python3 -E -B construction/common.py
```

## Release status

The complete default build and the 21-stage audit with `--native --links`
passed from a fresh relocated copy of this directory on 28 October 2026.
The measured times were 139.6 seconds for generation and 125.6 seconds for
the audit, using Python 3.12.14 on Linux x86-64. All three generated output
files were byte-identical to the accepted references. The link replay checked
all 4,310 links and 529,460 individual three-dimensional flips.

The complete stage reports and logs are in `verification/build/` and
`verification/audit/`. `verification/release_validation.json` summarizes
the exact code/input/output bindings and three negative checks: both drivers
reject optimized Python, and the audit rejects a changed roll parameter before
running any stages. These times describe this validation run, not a guarantee
for other machines.

This is the construction module prepared for repository review. Final
project authorship, acknowledgements, and licensing are to be agreed before
public publication. Existing source notices and provenance are retained;
this package does not grant a new license on behalf of its contributors.


# Motoo diagram construction, simplification, and certificates

This module packages the **188-crossing original-drawing construction**, denoted
`Z-circle` in Revision 6.1 of the paper. Its closed triangulation has **24,170
pentachora**. The recorded paths connect that exact triangulation to the explicit
standard sphere made from two four-simplices with identity facet pairings.

The module contains the construction and search programs, the fixed construction
data, and the complete reduction certificates. It combines `published_fig2_S4.zip`
and `published_fig2_closed_path_moves.tar.xz` with the required code and cap data
recovered from `motoo_figure2_simplification_v1.zip`. Imported files retain their
original bytes. Their origin is recorded in `provenance/`.

**This is the original-drawing route in the main proof.** The separate redraw
route starts with 25,200 pentachora and has different endpoints. The paper's
geometric argument relates the two presentations; this module does not establish
that relation by comparing their triangulation signatures.

## What is here

| Task | Entry point | Result |
| --- | --- | --- |
| Run Katie and construct the closed triangulation | `build.py` | Raw Katie output, finite truncation, boundary layers, reconstructed cap, and closed triangulation |
| Search from a small triangulation | `search.py small` | Recorded move paths; a result is checked before being reported as standard |
| Search from a large closed triangulation | `search.py large` | Chunked guarded-move journal for separate replay |
| Replay the released proof | `verify.py` | A terminal-only or full closed-triangulation standardization receipt |
| Check archive and endpoint bindings | `audit_bindings.py` | Structural checks, explicitly separate from move replay |
| Compare the original marked PD with the Katie-position PD | `audit_diagram.py` | Positive marked-exterior comparison using SnapPy |

Run the commands below from this directory. Each run needs a new output location.
Keep Python assertions enabled when invoking the public entry points.

## Environment

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r environment/requirements.txt
```

The package pins Regina 7.4.1. The terminal replay and small search use its Python
interface. The optional recognition accelerator for large runs additionally
requires Linux x86-64, the pinned CPython 3.12 Regina wheel, and a C++ compiler.
It checks the engine's exact hash and refuses a different engine. See
`environment/README.md` for details.

## A. Run Katie and construct the triangulation

On a compatible Linux x86-64 machine, the supplied installer obtains the pinned
Katie packages into this module's private `legacy/scratch/` directory. It does
not install system packages or require sudo:

```sh
.venv/bin/python legacy/code/setup_katie_linux.py
KATIE_PREFIX="$PWD/legacy/scratch/katie/prefix"
LD_LIBRARY_PATH="$KATIE_PREFIX/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  .venv/bin/python build.py --katie "$KATIE_PREFIX/usr/bin/katie" --output runs/build
```

Alternatively, provide an existing Katie 7.4.1 executable with `--katie`. The
construction uses the supplied marked PD code and framings, regenerates the
orientation convention, runs Katie, and checks each labelled interface against
the recorded construction. The raw Katie triangulation is a handlebody with an
ideal boundary vertex; it is not yet the closed sphere. The builder then performs
the explicit truncation, native simplification, fixed boundary layering, cap
reconstruction, and gluing.

`construction/README.md` gives the detailed stages and describes the distinct
mode that starts from the supplied raw Katie triangulation. Its result label
records that Katie itself was not rerun in that mode.

## B. Run a new simplification search

For a small triangulation, the default is the **eight-pentachoron endpoint of
this module's original-drawing construction**:

```sh
.venv/bin/python search.py small --output runs/small-search \
  --hours 12 --workers 4 --seed 1
```

Use `--start ISOSIG` or `--start-file FILE.rga` to choose another input. Both
root and nested Regina triangulation packets are supported. Resume a saved
small search with:

```sh
.venv/bin/python search.py small --resume runs/small-search \
  --hours 12 --workers 4
```

For a large closed triangulation, use the sampled move producer:

```sh
.venv/bin/python search.py large \
  --start-file runs/build/closed/closed_candidate.rga \
  --output runs/large-search --hours 12 --seed 11 --recognition-helper
```

The large producer saves a chunked journal. That output requires independent
replay; a producer receipt alone is not a verified standardness result. Its
historical implementation does not support resuming a large search, so this
interface does not promise that feature.

Searches are stochastic and may follow different paths. Reproduction of the
released proof uses the fixed certificates below and does not depend on finding
the path again. A failed bounded search does not determine the PL type.

## C. Replay the certificates

### Small terminal certificate

```sh
.venv/bin/python verify.py --scope terminal --output runs/terminal-replay
```

This checks all **430 moves** connecting the exact eight-pentachoron endpoint
to the explicit standard sphere. It verifies the move schema, every move guard,
all block endpoints, and the standard anchor. It does not report a replay of
the 24,170-pentachoron prefix.

### Complete closed-triangulation reduction

```sh
.venv/bin/python verify.py --scope full --recognition-helper \
  --source runs/build/closed/closed_candidate.rga \
  --output runs/full-replay
```

Omit `--source` to use the exact released closed triangulation. The supplied
source must match the certificate's labelled triangulation and complete
uncompressed RGA payload. Gzip compression bytes can differ between runtimes;
that case uses an explicit identity check to connect the fresh source with the
archival input of the unchanged legacy verifier. Both file hashes and the
identity check are recorded. The driver safely extracts the long
journal, replays all **780,847 moves**, binds its actual eight-pentachoron endpoint
to the terminal certificate, and verifies the additional **430 moves**. There
are **781,277 moves** in the composed certificate.

Without `--recognition-helper`, the full verifier uses Regina's original
recognition routines. That route can be substantially slower. The helper changes
only positive three-dimensional sphere/ball recognition: each accepted shortcut
is supported by an independently replayed group certificate, with Regina's
original routine used otherwise. The four-dimensional move guards are unchanged.

### Structural inspection

```sh
.venv/bin/python audit_bindings.py --output runs/bindings.json
```

This checks archive hashes, all journal chunk hashes and index/hash/size
transitions, endpoint encodings, and composition with the terminal certificate.
Its success result explicitly does **not** claim move legality or standardness;
those require the replay command.

### Optional marked-diagram comparison

The original 60-crossing marked PD is included from the reviewed paper package.
With the optional SnapPy dependency installed, compare it to the 188-crossing
Katie input:

```sh
.venv/bin/python -m pip install -r environment/requirements-diagram.txt
.venv/bin/python audit_diagram.py --output runs/marked-diagram.json
```

The check verifies the actual PD-to-component correspondence, equal decorated
exterior signatures, and the explicit component permutation with peripheral
matrices equal to plus or minus the identity. Framings and dotted roles must
agree under that permutation. This is a native SnapPy marked-exterior comparison;
it does not reread the original picture or reprove the paper's separate Kirby
slide argument.

## Certificate locations

| File | Meaning |
| --- | --- |
| `reference/published_fig2_S4/diagram/fig2_published_katie.json` | Marked 188-crossing Katie input |
| `reference/published_fig2_S4/katie/raw_katie.rga` | Recorded 1,528-pentachoron ideal output |
| `reference/published_fig2_S4/closed/closed_candidate.rga` | Recorded 24,170-pentachoron closed input |
| `reference/published_fig2_S4/boundary/path_short.json.gz` | Fixed 22,304-move boundary path |
| `certificates/published_fig2_closed_path_moves.tar.xz` | All 780,847 long-path moves and checkpoints |
| `reference/published_fig2_S4/final/CERTIFICATE.json` | The 430-move terminal certificate |
| `legacy/results/` | Cap construction and its fixed reduction/layering records |

The target signature is `cPkbbbbaaaaaaaa`. The terminal source signature is
`iLAALMQkabcccfefhghhhcaaavava2aaaoacbZb6aEa8a8a`.

## Scope and validation

The mathematical identification of the marked link with the drawing and the
paper's Kirby-move comparison remain separate obligations. The diagram data and
historical comparison receipts are preserved. Building this module's fixed
triangulation and replaying its reduction does not silently certify any other
drawing or the separate 25,200-pentachoron route.

The construction includes a native Regina simplification between truncation and
boundary layering. Its exact labelled output is checked. That native operation
is not expanded into an elementary move transcript. The subsequent closed
reduction uses explicit checked Pachner moves and their stated generalizations.

Fresh preparation results are in `verification/INTEGRATION_VALIDATION.json`.
Historical receipts inside the imported reference folders remain labelled by
their original filenames; they are not substituted for fresh execution.

The original archives' confidentiality and source notices are retained. This
module is an integration artifact; creating it does not publish the diagrams
or change the live repository.

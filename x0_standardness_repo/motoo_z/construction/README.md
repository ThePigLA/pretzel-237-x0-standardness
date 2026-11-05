# Constructing the closed triangulation of Z^circ

`../build.py` reconstructs the **24,170-pentachoron closed triangulation** used
by the main proof route for the manifold denoted `Z^circ` in the paper. It
starts with the admitted 188-crossing, eight-component marked PD code, invokes
Katie, removes the ideal vertex by an explicit prism construction, builds a
standard upper piece, and glues the finite pieces along their actual
triangulated boundaries. Each numbered construction state must agree with
the published labelled triangulation and its uncompressed serialized bytes.
The archived file hashes and the actual generated file hashes are recorded
separately because different zlib builds can produce different gzip bytes.

This command produces the closed starting complex. The subsequent reduction
to the standard sphere is a separate replay described in the module README.
The earlier 25,200-pentachoron triangulation of the redraw is a different
construction and is not the output of this wrapper.

## Requirements

Use Python 3.10 or newer and **Regina 7.4.1**:

```sh
python3 -m venv regina_env
. regina_env/bin/activate
python -m pip install regina==7.4.1
```

The complete construction also requires the **Katie 7.4.1 executable**.
Installing Regina's Python package does not by itself provide that command.
Pass the executable explicitly with `--katie`. Its dependent libraries must
be available in the inherited environment. SnapPy is not needed for this
construction; the separate diagram-identification checks have their own
dependencies.

Run without `-O` or `-OO`. The wrapper refuses optimized Python because the
unchanged historical code contains assertions that must remain enabled.

## Full construction, including Katie

From the `motoo_z` module directory:

```sh
python build.py --output runs/z_from_pd_01 --katie /absolute/path/to/katie
```

`--out` is an alias of `--output`. The output directory must not exist. No
released input is overwritten. The command first checks the release manifest;
the logs and all fresh receipts go into the chosen output directory.

On success, the report has status:

```text
PASS_CONSTRUCTION_FROM_MARKED_PD_WITH_KATIE
```

The principal output is:

```text
runs/z_from_pd_01/closed/closed_candidate.rga
```

The archived reference file's SHA-256 is:

```text
e6281897c358e315757a77f851433db1527f48d71b97ec0056174e0a7eb2c467
```

The regenerated file can have a different SHA-256 solely because of gzip
compression. The wrapper requires its decompressed XML and its exact labelled
triangulation to match the archived reference. It records the actual file hash
in `BUILD_REPORT.json` and binds the fresh admission receipt to that hash.

The wrapper records the actual Katie executable hash and version. The
historical executable was compiled separately; a different compiler or the
official package need not produce the same executable bytes. The required
comparison is the generated triangulation: both its exact labelled state and
the decompressed serialized file must match. The wrapper does not infer a binary's source
provenance merely from the reported version number.

### Optional pinned Linux Katie installation

The unchanged `legacy/code/setup_katie_linux.py` fetches three packages with
pinned SHA-256 values and extracts them locally. It targets Ubuntu 24.04 on
x86-64 and requires `dpkg-deb` plus the usual system C/C++ libraries. It does
not use `apt` or `sudo` or change system packages.

To keep the released source tree separate from downloaded executables, copy
the setup script into a fresh working directory before invoking it:

```sh
mkdir -p runs/katie_install/code
cp legacy/code/setup_katie_linux.py runs/katie_install/code/
python runs/katie_install/code/setup_katie_linux.py
```

Then, from the module directory:

```sh
KATIE_PREFIX="$PWD/runs/katie_install/scratch/katie/prefix"
LD_LIBRARY_PATH="$KATIE_PREFIX/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  python build.py --output runs/z_from_pd_01 \
  --katie "$KATIE_PREFIX/usr/bin/katie"
```

The setup receipt records the download URLs, package hashes, and executable
hash. Other systems can use their own Katie 7.4.1 installation. The wrapper
always passes arguments directly to the executable and inherits its library
environment.

## Reconstruction beginning with the supplied raw Katie output

When Katie is unavailable, all subsequent stages can still be reconstructed:

```sh
python build.py --output runs/z_from_raw_01 --from-raw
```

This mode explicitly copies the pinned 1,528-pentachoron raw Katie file into
the new run directory. It recomputes the marked-PD orientation checks and all
later construction stages, but **does not rerun Katie's conversion**. Its
distinct success status is:

```text
PASS_RECONSTRUCTION_FROM_SAVED_RAW_KATIE
```

A failed Katie invocation never falls back automatically to this mode.
Neither mode substitutes the saved lower piece or saved upper cap for a
failed regeneration.

## What is reconstructed and checked

| Stage | Pentachora | How the state is produced and checked |
|---|---:|---|
| Raw Katie output | 1,528 | Recompute dotted-component orientation, invoke Katie on the marked PD, and require the exact published output; `--from-raw` explicitly starts at this state. |
| Finite truncation | 6,112 | Replace each truncated old simplex by four staircase cells; separately check every cell and facet pairing. The ideal link is recognized as `#2(S^1 x S^2)` and the other links as 3-spheres. |
| Simplified lower piece | 1,392 | Call Regina's native `simplify()` once on the regenerated finite piece; require the exact published state and its 826-tetrahedron boundary. |
| Explicit upper handlebody | 24 | Construct two copies of a staircase triangulation of `Delta^3 x S^1` and join them along the recorded embedded boundary tetrahedron. |
| Small upper handlebody | 9 | Replay all 15 guarded edge collapses and compare every labelled state; the independent cap checker reconstructs the product separately. |
| Upper piece with matched boundary | 474 | Recompute all 465 recorded boundary Pachner attachments and their affine facet maps; undo them by guarded boundary shellings to recover the generated nine-cell piece. |
| Lower piece with matched boundary | 23,696 | Replay 22,304 recorded 3D Pachner moves, attaching one 4-simplex per move. Independently check the attaching-ball and collar conditions. |
| Closed triangulation | 24,170 | Recompute the six-tetrahedron boundary isomorphism and glue the two generated pieces. Independently audit the gluing and all vertex links. |

The final closed f-vector is `(19, 12136, 48374, 60425, 24170)` and its Euler
characteristic is 2. The final admission reconstructs the 19 vertex links
from raw face pairings and checks their group-elimination certificates. It
does not recognize the closed 4-manifold as the standard sphere; that
conclusion uses the separate reduction certificate.

The lower boundary path is fixed input from
`reference/published_fig2_S4/boundary/path_short.json.gz`. The upper path is
recovered from the full recorded events in
`legacy/results/boundary/standard_cap_six_boundary_v1/layering_events.jsonl.gz`.
The wrapper writes the derived upper path into the output directory and
recomputes every recorded attachment. It does not rerun the searches that
originally found either path.

## Construction boundaries and provenance

The starting PD code and integer framings are bound by the published digest
in `reference/published_fig2_S4/diagram/fig2_published_katie.json`. This command
does not re-extract curves from a PDF or prove the paper's geometric
identification between different drawings. Those belong to the diagram
checks and the geometric argument in the paper.

The historical Katie runner first made a small adapter JSON containing the
PD code, decorations, and component names. The original bytes of that adapter
were not included in the uploaded archive. The wrapper reconstructs those
fields from the admitted PD, records their fresh digest, and repeats the
orientation computation. The PD code, decorations, component names, and all
geometric orientation checks must agree exactly with the supplied oriented
input; only the old adapter-file digest is replaced with the actual fresh
digest. The generated Katie triangulation and its decompressed serialized
file must still match the published reference exactly.

The source for Katie is preserved at `legacy/source/katie_7_4_1.cpp`. The
legacy Python files and the published reference files are preserved
unchanged. Their original receipts can contain historical absolute paths;
fresh construction receipts contain paths to the newly generated pieces, so
the active workflow does not depend on those historical locations.

Native Regina `simplify()` in the 6,112-to-1,392 step remains an explicitly
identified software trust boundary. Its individual operations were not
recorded in the construction archive. The wrapper reruns the call and checks
the exact result. By contrast, the cap's 24-to-9 reduction is an individually
guarded move replay and the two boundary-layering sequences are checked move
by move.

## Logs, failure handling, and further reduction

Each stage's complete output is retained in `logs/`. `STAGES.json` records
commands, exit codes, log hashes, and timings. `BUILD_REPORT.json` records the
overall result, exact-state checkpoints, and whether Katie was executed.
An interrupted or failed run does not receive a passing construction status.

The two adjustable time budgets are `--katie-timeout` (default 600 seconds)
and `--group-seconds` (default 120 seconds for each vertex-link group search).
There is no resume option for construction. If a stage fails, keep the
directory for diagnosis and use a new output directory for a corrected run.

For the standardness proof, keep the generated closed file together with the
build report and run the module's certificate replay. The supplied reduction
is tied to the archived closed-file hash above; the construction receipt
establishes exact labelled-state and decompressed-file identity between that
archived input and the regenerated complex. Construction and reduction have
separate success records so that beginning from a saved intermediate cannot
be mistaken for regenerating the full pipeline.

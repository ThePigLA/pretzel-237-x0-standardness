# Runtime and Katie setup

## Python and Regina

`requirements.txt` pins the Python distribution `regina==7.4.1`. Regina's
`versionString()` reports `7.4` for this distribution; the wrappers check the
installed distribution version separately.

The terminal verifier and small-search wrapper use ordinary Regina Python APIs.
The provided input-file wrapper accepts either a root triangulation packet or a
single triangulation nested in a packet tree.

## Katie

`../legacy/code/setup_katie_linux.py` fetches three recorded Ubuntu 24.04 x86-64
packages, checks their exact SHA-256 values, and extracts them below
`../legacy/scratch/katie/`. It does not use apt or sudo. The host must provide
the compatible system libraries described in that script.

```sh
python legacy/code/setup_katie_linux.py
KATIE_PREFIX="$PWD/legacy/scratch/katie/prefix"
LD_LIBRARY_PATH="$KATIE_PREFIX/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  python build.py --katie "$KATIE_PREFIX/usr/bin/katie" --output runs/build
```

For another installation, pass `--katie /absolute/path/to/katie` to `build.py`
and supply any library search path required by that installation. Katie must
report version 7.4.1. The builder checks the actual resulting labelled
triangulation, so a different binary does not bypass the input/output bindings.

The exact Katie 7.4.1 source used by the earlier package is retained at
`legacy/source/katie_7_4_1.cpp`, with its original notices. The historical
certificate used a locally compiled executable; the documented installer uses
the pinned upstream Ubuntu package. Their executable hashes differ; the
reproduced mathematical output must agree.

## Optional recognition accelerator

The unchanged launcher `legacy/code/run_with_primal_recognition.py` pins the
Linux x86-64 CPython 3.12 Regina 7.4.1 engine to:

```text
42e0c5a630d71f98c6e5ff3c70c520d1c3d136151aad4ac37fe807004fbaca1d
```

The public wrapper compiles its C++ bridge in the output directory and runs the
supplied positive and negative controls before using it. It intercepts only
three-dimensional sphere/ball recognition and leaves four-dimensional move
guards unchanged. Every shortcut is accepted only after independent replay of
its group certificate; inconclusive cases fall back to original recognition.

Do not edit the engine pin to make a different binary pass. Use original Regina
recognition without the helper, or perform a separate validation of the new
runtime. The archived Python 3.13 helper in
`reference/published_fig2_S4/code/helper/` is retained as provenance; the public
wrapper uses the unchanged CPython 3.12 launcher from the earlier package.

## Apple laptops

The small terminal replay and small search can be run in a Regina 7.4.1 Python
environment. The Ubuntu Katie installer and the ELF recognition accelerator
target Linux. A Mac build needs a separately installed Katie 7.4.1 executable;
the helper-assisted full replay should be run in the matching Linux environment.
This preparation records its actual tested platform in
`verification/INTEGRATION_VALIDATION.json`.

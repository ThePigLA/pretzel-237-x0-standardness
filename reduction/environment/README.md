# Reduction runtime

The complete reduction uses CPython 3.12 and Regina 7.4.1 on Linux x86_64.
The native macro segment checks a numbered gluing hash after every recorded
block. Matching the recorded native implementation is therefore part of
reproducing that segment.

The runtime pins are:

| Component | Pin |
| --- | --- |
| Container platform | `linux/amd64` |
| Python base | `python:3.12.14-slim-bookworm` |
| Base image manifest SHA256 | `1aaa65a85fda306ffb8b910824d4e93bdce61e212c7e87168123ea3073b41a1a` |
| Regina wheel | `regina-7.4.1-cp312-cp312-manylinux_2_28_x86_64.whl` |
| Wheel SHA256 | `d33ccc69697ec680a8dcf8f0c663428853da44f3578898ad9b908e4669e83c71` |
| Regina native engine SHA256 | `42e0c5a630d71f98c6e5ff3c70c520d1c3d136151aad4ac37fe807004fbaca1d` |

`runtime_pins.json` records the public metadata sources, resolved image and
configuration digests, and the Python version read from the image configuration.
Both registry digests were checked against the downloaded manifest/configuration
bytes. The image is pinned by digest: the tag is a descriptive label rather than
a request to use whichever image it points to later. See Docker's
[base-image pinning documentation](https://docs.docker.com/build/building/best-practices/#pin-base-image-versions).
The wheel hash comes from the
[Regina 7.4.1 PyPI metadata](https://pypi.org/pypi/regina/7.4.1/json).

**Validation boundary:** the native Linux preparation environment has the
recorded Regina engine. The full replay receipts identify the environment used
for their checks. Docker is not installed in that preparation environment, so
this container recipe has not itself been built or run there. Its build includes
an explicit interpreter, platform, package-version, and engine-hash check. A
successful image build does not replace the complete reduction replay.

## macOS: use the Linux runtime

Run the following commands from the directory containing `construction/` and
`reduction/`, with Docker available. The construction candidate is assumed to be
at `runs/x0-build/candidate`.

```sh
docker build --platform linux/amd64 \
  -f reduction/environment/Dockerfile \
  -t x0-reduction-runtime:1 .
```

```sh
docker run --rm --init --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$(pwd),dst=/work" \
  --workdir /work \
  --env PYTHONDONTWRITEBYTECODE=1 \
  x0-reduction-runtime:1 \
  python -E -B reduction/verify.py \
    --candidate runs/x0-build/candidate \
    --output runs/x0-reduction
```

The output directory must be new. The bind mount keeps the replay receipts and
endpoint files on your machine, and `--user` runs the verifier with your numeric
user/group IDs so it does not create root-owned output. Docker documents these
[bind-mount options](https://docs.docker.com/engine/storage/bind-mounts/) and
the [`--user`, `--platform`, and `--workdir` flags](https://docs.docker.com/reference/cli/docker/container/run/).
The explicit `linux/amd64` platform selects the recorded architecture on both
Intel and Apple Silicon Macs. Performance under architecture emulation can
differ from the recorded Linux timing.

If the construction candidate is outside the repository, mount its directory
separately. Replace the example path below with its absolute host path:

```sh
docker run --rm --init --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$(pwd),dst=/work" \
  --mount "type=bind,src=/absolute/path/to/candidate,dst=/candidate,readonly" \
  --workdir /work \
  --env PYTHONDONTWRITEBYTECODE=1 \
  x0-reduction-runtime:1 \
  python -E -B reduction/verify.py \
    --candidate /candidate \
    --output runs/x0-reduction
```

The Dockerfile installs only the runtime. Repository sources and certificates
are supplied through the mount, so rebuilding the image is unnecessary when
changing documentation or choosing a different new output directory.

## Linux x86_64: use a virtual environment

With CPython 3.12 installed, run from the repository root:

```sh
python3.12 -m venv .venv-reduction
.venv-reduction/bin/python -m pip install \
  --require-hashes --only-binary=:all: --no-deps \
  -r reduction/environment/requirements.txt
.venv-reduction/bin/python -E -B reduction/verify.py \
  --candidate runs/x0-build/candidate \
  --output runs/x0-reduction
```

The requirements file deliberately permits only the recorded CPython 3.12 Linux
x86_64 wheel. It is not a native macOS installation file. pip's
[hash-checking documentation](https://pip.pypa.io/en/stable/topics/secure-installs/)
explains the required hash and binary-only options.

## Partial verification scopes

The elementary link, contraction, and simplicial-move replay uses only Python's
standard library on a supported POSIX system:

```sh
python3 -E -B reduction/verify.py --scope elementary \
  --candidate runs/x0-build/candidate \
  --output runs/x0-elementary
```

The terminal explicit-move certificate can be replayed separately with a
compatible Regina installation:

```sh
python -E -B reduction/verify.py --scope terminal \
  --output runs/x0-terminal
```

The terminal scope does not use a construction candidate and does not rerun the
earlier native macro segment. Each partial scope reports its narrower result.
Use the default full scope to verify the complete supplied chain from the
construction master to the standard sphere.

Run these entry points without `-O` or `-OO`; some unchanged elementary checkers
use Python assertions. The public driver rejects optimized execution.

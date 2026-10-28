# Runtime only. Build from the repository root with --platform linux/amd64.
# Base manifest and Python version were verified against the public registry.
# The container build itself has not been executed in the preparation environment.
FROM python:3.12.14-slim-bookworm@sha256:1aaa65a85fda306ffb8b910824d4e93bdce61e212c7e87168123ea3073b41a1a

ENV PYTHONDONTWRITEBYTECODE=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# A URL requirement carries its locally pinned SHA256 into pip's hash-checking mode.
# No repository content is copied into this image; mount the repository at /work.
RUN python -m pip install --no-cache-dir --no-deps --only-binary=:all: --require-hashes \
    "https://files.pythonhosted.org/packages/68/42/b4a994e393298c08c361342ee19bf240b7b16f155c492ce0631082c14e66/regina-7.4.1-cp312-cp312-manylinux_2_28_x86_64.whl#sha256=d33ccc69697ec680a8dcf8f0c663428853da44f3578898ad9b908e4669e83c71"

RUN ["python", "-c", "import hashlib, importlib.metadata, pathlib, platform, sys; import regina; actual = hashlib.sha256(pathlib.Path(regina.engine.__file__).read_bytes()).hexdigest(); expected = '42e0c5a630d71f98c6e5ff3c70c520d1c3d136151aad4ac37fe807004fbaca1d'; ok = sys.version_info[:3] == (3, 12, 14) and platform.system() == 'Linux' and platform.machine() == 'x86_64' and importlib.metadata.version('regina') == '7.4.1' and actual == expected; sys.exit(0 if ok else 'Unsupported reduction runtime or Regina engine hash mismatch')"]

WORKDIR /work
CMD ["python", "-E", "-B", "reduction/verify.py", "--help"]

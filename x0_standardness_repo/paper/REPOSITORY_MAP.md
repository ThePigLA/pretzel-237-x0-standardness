# Appendix A and the repository layout

The manuscript retains the original C1–C4 archive designations so that its
accepted statements, file identities, and commands remain unchanged. The
repository organizes their construction and proof material into modules
with coordinated entry points.

| Paper designation | Repository location | Current entry point |
| --- | --- | --- |
| C1: marked construction dossier | [`construction/`](../construction/README.md) | `construction/build.py` and `construction/audit.py` |
| C2: reduction of the master to the closed six-pentachoron interface | [`reduction/elementary/`](../reduction/elementary/) and [`reduction/source/native/`](../reduction/source/native/) | Elementary and native stages of `reduction/verify.py` |
| C3: terminal certificate and supplied verifier | [`reduction/certificates/standardization.json`](../reduction/certificates/standardization.json) and [`reduction/source/terminal/verify_x0_standard.py`](../reduction/source/terminal/verify_x0_standard.py) | Included in full or terminal-scope reduction replay |
| C4: independent terminal checks | [`reduction/source/terminal/`](../reduction/source/terminal/) and [`reduction/terminal.py`](../reduction/terminal.py) | Independent replay, separate two-zero reconstruction, and fresh composition in `reduction/verify.py` |

The two original terminal archives are also included intact:

- [x0_is_standard_v1.1.zip](supplements/x0_is_standard_v1.1.zip)
- [x0_standardness_audit_v1.zip](supplements/x0_standardness_audit_v1.zip)

The construction and reduction directories are curated repository components,
not byte-for-byte copies of the old outer C1 and C2 ZIP files. Their own
provenance inventories identify the exact imported source members and bytes.
Their new drivers add fresh computation and explicit checks between stages.
The hashes of the original archives printed in the paper continue to identify
those original archives; they are not hashes of the new repository snapshot.

## Recommended repository sequence

Run these commands from the repository root. Each output directory must be new.

```sh
python3 -E -B construction/build.py --output runs/x0-build
```

Run the construction audit in the environment described in
[`construction/README.md`](../construction/README.md):

```sh
python -E -B construction/audit.py \
  --candidate runs/x0-build/candidate \
  --native --links \
  --output runs/x0-construction-audit-full
```

Run the complete reduction in its recorded Linux environment:

```sh
python -E -B reduction/verify.py \
  --candidate runs/x0-build/candidate \
  --output runs/x0-reduction
```

The required setup and Mac container commands are in
[`reduction/environment/README.md`](../reduction/environment/README.md).
The public entry point rejects optimized Python because the inherited
elementary checks rely on assertions. Appendix A's `python -O` examples
apply only to its named assertion-free terminal scripts; use the commands
above for the coordinated repository replay.

## Check the standalone terminal theorem

With Regina 7.4.1 installed, on a Mac or Linux:

```sh
python -E -B reduction/verify.py \
  --scope terminal --output runs/x0-terminal
```

This executes the supplied terminal checker, the independent full replay,
the independent two-zero checker, and the terminal composition check. It
establishes the terminal theorem's certificate without rebuilding the
construction or replaying the earlier reduction. Its result is explicitly
labelled `PASS_TERMINAL_ONLY`.

The repository's full driver additionally connects the freshly generated
capped endpoint to those terminal paths. Its new composition check uses
the receipts produced during the current run; the historical binding
script inside C4 remains preserved in the original archive.

## Shared exact identities

The terminal certificate has the same bytes in C3, C4, and
`reduction/certificates/standardization.json`:

```text
e9ba9f14831bb0da32043d7d93eb0148f0683e2a6e4f67d5156e8e9fc94adcdf
```

The marked master consumed by reduction has file SHA-256:

```text
251d7ad38235869fa268e9ee6ec9da2f397cf1ca71db267d0bf0ca9b434c5e4c
```

These identities are checked by the module entry points and their fresh
validation records. Compilation of the paper alone does not replay them.

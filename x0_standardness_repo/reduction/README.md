# Reduction of the marked X0 master to the standard sphere

This repository component replays the accepted reduction of the exact
**103,744-pentachoron, 4,310-vertex master** produced by `construction/`.
It carries the proof through to the explicit standard sphere obtained by
doubling one four-simplex. The module includes the readable checkers, the
complete fixed move certificates, and fresh reproduction receipts.

The main command performs an actual fresh replay. It checks the interfaces
between stages using the files produced in that run. The geometric
identification of the marked master with the intended double branched cover
is the responsibility of the separate `construction/` module.

## Full reproduction

The complete replay requires **CPython 3.12 on Linux x86-64 with Regina
7.4.1**, using the recorded native engine. This requirement comes from the
punctured macro stage. Installation instructions and a digest-pinned
container recipe for a Mac are in [environment/README.md](environment/README.md).

From the repository root, with that environment active:

```sh
python -E -B reduction/verify.py \
  --candidate runs/x0-build/candidate \
  --output runs/reduction
```

Here `runs/x0-build/candidate` is the output of:

```sh
python3 -E -B construction/build.py --output runs/x0-build
```

The candidate must contain the exact released `marked_candidate.json.gz`
and `cover_regina.json.gz`. Both file hashes are checked before any replay
output is created. A modified or different triangulation requires a new
certificate; this command does not search for one. Compressed-byte identity
is intentional because the inherited certificates bind those exact files.

The module is also self-contained. Omitting `--candidate` replays the supplied
copy of that exact master:

```sh
python -E -B reduction/verify.py --output runs/reference-reduction
```

The report distinguishes an explicitly selected construction output from the
packaged reference. Both are subjected to the same mathematical replay.

Use a new output directory for each run. The driver refuses to overwrite an
existing one, protects the source and input directories, and preserves the
logs from a failed or interrupted attempt. Successful runs remove their
temporary input copies; `--keep-work` retains them if needed.

**Do not use `python -O` or `python -OO`.** Some imported elementary checkers
use assertions for mathematical preconditions. The public driver rejects
optimization and launches its checkers with assertions enabled.

## What the command proves

| Stage | Result | Evidence checked in the fresh run |
| --- | --- | --- |
| Manifold premise | All 4,310 vertex links are PL three-spheres | 529,460 recorded three-dimensional bistellar flips |
| Prefix contractions | 103,744 to 44,620 pentachora | 3,308 contractions, each satisfying the full link condition |
| Simplicial reduction | 44,620 to 1,332 facets, with 45 vertices | 629,874 moves: 957 contractions and 628,917 bistellar flips |
| Native handoff | The same closed 1,332-facet triangulation | Exact ordered facets and all face gluings; every half-facet accounted for |
| Puncture and native reduction | 1,331 to 5 pentachora, with boundary | An embedded simplex puncture, 54,563 guarded 2–4 moves, and 5,173 native simplify blocks with an exact numbered-gluing hash after every block |
| Capping | A closed six-pentachoron triangulation | Explicit boundary isomorphism and embeddings, face pushout, and direct endpoint isomorphism without further simplification |
| Terminal standardization | PL equivalence to the explicit two-pentachoron standard sphere | Three connected paths totaling 8,436 checked moves, including a separate implementation checking all 594 admissible 2–0 moves |

The terminal certificate consists of the forward paths
`G6 -> W`, `S6 -> W`, and `S2 -> S6`. Their exact common states establish
PL equivalence between `G6` and `S2`. The JSON is a connected collection of
checked paths, rather than one forward list of decreasing-size moves.

`S2` is constructed directly from two distinct abstract four-simplices,
joining each corresponding pair of facets by the identity. Its signature is:

```text
cPkbbbbaaaaaaaa
```

The closed six-pentachoron interface has signature:

```text
gAMLMQaabcdeffff2a+a2avbmaoaxaOb0aUa
```

The character before `b0aUa` is capital `O`. Neither endpoint's role depends
on a classification lookup.

The top-level result for a completed full run is:

```text
PASS_MASTER_PL_HOMEOMORPHIC_TO_STANDARD_S4
```

This is a nonequivariant PL equivalence. The move certificates do not carry
the deck involution or the marked surfaces through the reductions.

## Trust boundaries

The elementary portion uses only Python's standard library. The native
portion uses Regina's move, triangulation, and signature implementations.
In particular, the 5,173 simplify blocks record their exact outputs but do
not expand each internal simplify operation into an elementary move list.

The terminal certificate is more explicit: all 8,436 moves replay separately.
Its independent 2–0 checker does not call Regina's `has20`, `move20`, or
`twoZeroMove`; it reconstructs the resulting facet pairings. It still shares
Regina's signature codec, skeleton queries, and triangulation constructor.

The complete chain is therefore **not a pure Pachner transcript** and is not
a proof-assistant formalization. Its checked computations and mathematical
dependencies are listed in [docs/CHAIN_AND_TRUST.md](docs/CHAIN_AND_TRUST.md).

## Shorter modes

The elementary part can be replayed with ordinary unoptimized Python on a
POSIX system, without Regina:

```sh
python3 -E -B reduction/verify.py \
  --scope elementary \
  --candidate runs/x0-build/candidate \
  --output runs/reduction-elementary
```

The explicit small terminal paths can be replayed independently on a Mac
or Linux with Regina 7.4.1:

```sh
python -E -B reduction/verify.py \
  --scope terminal \
  --output runs/reduction-terminal
```

These modes produce `PASS_ELEMENTARY_ONLY` and `PASS_TERMINAL_ONLY`,
respectively. Neither reports that the complete master-to-sphere chain
was replayed. If Regina is installed in a separate interpreter, select it
with `--native-python /absolute/path/to/python`.

For an integrity check with no mathematical replay:

```sh
python3 -E -B reduction/verify.py --integrity-only
```

## Outputs

| Path inside the run directory | Contents |
| --- | --- |
| `REDUCTION_REPLAY.json` | Overall scope and result, exact input and code identities, fresh stage results, interface bindings, and timing |
| `elementary/` | Link/contraction/flip logs and resource records, plus the actual replayed 1,332-facet export |
| `native/` | Replayed punctured endpoint, capped endpoint, all block sizes, capping maps, and native receipt |
| `terminal/independent_replay.json` | All three move-path results and the explicit standard anchor |
| `terminal/independent_two_zero.json` | Separate replay of all 594 2–0 operations |
| `terminal/bindings.json` | Exact composition of the fresh cap and the terminal paths |
| `terminal/standard_sphere.json` | The explicit two-pentachoron standard sphere's facet-pairing table |
| `logs/` | Complete output from every invoked checker |

The driver hashes both fresh receipts and fresh endpoint files. It checks
that all selected inputs and code are unchanged when the run finishes.
An old successful receipt is never substituted for executing a stage.

## Directory guide

| Path | Purpose |
| --- | --- |
| `verify.py` | Coordinated replay, fresh-output requirements, and interface checks |
| `common.py` | Exact input identities, provenance checks, and path protection |
| `terminal.py` | Independent terminal composition and explicit endpoint checks |
| `source/elementary/`, `source/native/`, `source/terminal/` | Unchanged imported audit implementations |
| `elementary/` | Elementary move checkers, their fixed input files, complete move logs, and native macro certificate |
| `certificates/standardization.json` | The unchanged complete terminal certificate |
| `reference/original/` | Original-release input copies and manifest, for exact source binding |
| `provenance/` | Per-file source archive, member name, size, and SHA-256 mappings |
| `environment/` | Pinned dependencies, container recipe, and environment records |
| `verification/` | This release's fresh reproduction and negative-check results |

Some source inputs appear twice because the accepted elementary audit checks
the replay input against the original release and its manifest. The duplicated
files are byte-identical. The small historical control master is retained for
an ancillary cap-difference check; its standardization is not a premise of
this reduction. Imported source bytes and existing notices are preserved.

## Release validation

The complete six-stage replay passed from a relocated copy on 28 October
2026 in **306.7 seconds**, using a master freshly regenerated by
`construction/build.py`. The fresh geometry build took 128.7 seconds and
reproduced the released master and exports byte-for-byte. All three rejection
checks passed: optimized Python, a changed master, and a changed terminal
certificate were rejected before any replay output was created.

The validation details and measured timing are in
`verification/release_validation.json`. They identify the freshly generated
construction output, the relocated replay directory, every code hash,
the exact native engine, and the resulting standard-sphere export.

This component is prepared for repository review. Project authorship,
acknowledgements, and licensing remain to be agreed before public publication.
No new license is assigned on behalf of the original contributors.

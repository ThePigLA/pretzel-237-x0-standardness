# Paper

This directory contains **Revision 5, dated 28 October 2026**, of
*Standardness of the double branched cover of the zero-twist, one-roll spin
of P(-2,3,7): a computer-assisted proof with explicit move certificates*.

**Authors:** Vance Faber and Michael Murphy, in that order.

Revision 5 adds the author line and an acknowledgement of substantial
assistance from AI tools developed by OpenAI and Anthropic. The authors
remain responsible for the mathematics, computational verification, and
manuscript. The mathematical statements and proof data are unchanged from
Revision 4; the bibliography and two terminal supplement archives retain
their accepted bytes.

## Read and edit

| File | Role |
| --- | --- |
| [x0_standardness_proof_note.pdf](x0_standardness_proof_note.pdf) | Current manuscript reading copy |
| [x0_standardness_proof_note.tex](x0_standardness_proof_note.tex) | Primary editable manuscript |
| [x0_standardness_proof_note.txt](x0_standardness_proof_note.txt) | Searchable manuscript text |
| [bibliography.bib](bibliography.bib) | The fourteen scholarly references as editable bibliographic data |
| [REPOSITORY_MAP.md](REPOSITORY_MAP.md) | Correspondence between Appendix A's original archives and the current modules |

Edit the `.tex` file for the next manuscript revision. Its bibliography is
embedded in a `thebibliography` environment, so BibTeX is not needed for the
current build. The `.bib` file is retained for reference management and a
future journal-format conversion; edits to it do not automatically change
the embedded bibliography. Keep both consistent when revising references.

The text export is an accompanying reading format, not a second source from
which the `.tex` is regenerated. The new build command compiles the `.tex`
directly and leaves the distributed source and PDF untouched.

## Build

From the repository root:

```sh
python3 -E -B paper/build.py --output runs/paper-build
```

Use a new output directory. The command copies the manuscript into a temporary
working directory and runs pdfLaTeX twice with shell escape disabled. It
retains compiler logs, records warnings and source identities, and writes:

- `x0_standardness_proof_note.pdf`: the newly compiled reading copy.
- `x0_standardness_proof_note.txt`: text extracted from that PDF, when
  Poppler's `pdftotext` is available.
- `BUILD_REPORT.json`: compiler version, source/output hashes, pass results,
  warnings, timing, and reference comparison.
- `logs/`: the compiler logs and available PDF metadata.

Required tools are Python 3 and a TeX installation providing `pdflatex`,
Latin Modern, `amsmath`, `amssymb`, `booktabs`, `array`, `geometry`,
`microtype`, `fancyhdr`, `fancyvrb`, `needspace`, `hyperref`, and
`glyphtounicode`. No topology software is needed to compile the paper.

For a reproduction check against the supplied reading copy, install
Poppler's `pdfinfo` and `pdftotext` and run:

```sh
python3 -E -B paper/build.py \
  --output runs/paper-reference-check \
  --verify-reference
```

That option requires matching page counts and normalized extracted text.
It is intended for reproducing this revision. The ordinary build permits
intentional edits and reports whether its output still matches the supplied
reference. PDF timestamps and identifiers may differ even when the pages
render identically.

The compilation check does not execute a mathematical certificate. The
[construction module](../construction/README.md) and
[reduction module](../reduction/README.md) provide those entry points.

## Reproduction and provenance

Appendix A retains the original C1–C4 designations and archive hashes.
[REPOSITORY_MAP.md](REPOSITORY_MAP.md) gives their current locations and
the repository commands. The original C3 and C4 archives remain in
`supplements/` with their accepted bytes intact.

`provenance/imported_files.json` records the historical import from the
Revision 4 archive, including the original member names and hashes.
`provenance/revision5_changes.json` records the subsequent authorship and
acknowledgement changes. The original archive manifest and its earlier
typography-validation receipt are retained under `provenance/`.

The current reading copy is rebuilt from the Revision 5 LaTeX source.
Build results, compiler logs, and layout-review records are kept separately
under `verification/`, including `verification/build/BUILD_REPORT.json`
and `verification/release_validation.json`. These checks concern the
current manuscript; they do not assert page-pixel identity with Revision 4.

## Publication decisions

The author attribution and acknowledgement of AI assistance are set for
this draft. Affiliations, any further acknowledgements, licensing, public
repository identity, and archival citation metadata remain open. The
scholarly bibliography is complete for this draft; repository citation
metadata and DOI links will be added after the publication decisions are
settled.

Before a public manuscript release, Appendix A's distribution paragraph can
be updated to the final repository and archive locations. The mapping
document supports use of the present draft without inventing those locations.

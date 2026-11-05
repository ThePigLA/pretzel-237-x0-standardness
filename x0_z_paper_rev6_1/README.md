# Combined paper, Revision 6.1

**5 November 2026 — confidential working draft for coauthor review**

Authors: **Vance Faber, Michael Murphy, Motoo Tange**, in the order retained from the supplied combined Revision 5.

Open `x0_z_standardness_note_rev6_1.pdf` for the paper. The complete editable source is `x0_z_standardness_note_rev6_1.tex`, with four vector PDF figures in this directory. `x0_z_rev6_1_new_material_marked.pdf` marks everything that is not in the accepted X0-only Revision 6 of 28 October, which is the useful comparison for a reader who knows only the X0 paper.

Revision 6.1 is a small follow-up to Revision 6; see the top of `REVISION_NOTES.md`.

This revision merges the new Motoo Tange / Z material in combined Revision 5 with the accepted X0-only Revision 6 of 28 October 2026. The two revision-number sequences refer to different documents. See `REVISION_NOTES.md` for the exact restoration, the new control discussion, and remaining author items.

## Build

From this directory, with a standard TeX Live installation:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error x0_z_standardness_note_rev6_1.tex
```

Alternatively, run `pdflatex -interaction=nonstopmode -halt-on-error x0_z_standardness_note_rev6_1.tex` twice. No shell escape, external bibliography processor, or network access is needed. The distributed PDF was built with pdfTeX / TeX Live 2023, has 23 pages, and embeds the mathematical fonts. The four figures remain byte-for-byte identical to the Revision 5 figures.

## Included in this manuscript bundle

- Current PDF and LaTeX source; four vector figures.
- `corrected_figures_check/`: the original four labelled PD records and historical comparison/extraction scripts, plus the new exact local diagram audit. Its own README explains the scope and status of each script.
- `history/`: combined Revisions 5 and 6 (PDF and source), the two Revision 5 marked-change PDFs, and the recovered X0-only Revision 6 source used for restoration.
- `provenance/`: input hashes, the restoration diff, the complete Revision 5-to-6 source diff, the control archive metadata/count check, the historical Revision 6 validation/build records, and the separate fresh Revision 6.1 review record.
- `SHA256SUMS`: hashes of the distributed files other than the manifest itself.

The manuscript ZIP is not the complete computational proof distribution. In particular, the separate supplements C1–C6 and the large Z proof archives are **not** included. This is stated explicitly in Appendix A.

## Separate proof data

For X0, retain the marked construction dossier, `x0_to_G6_verified_v1.zip`, `x0_is_standard_v1.1.zip`, `x0_standardness_audit_v1.zip`, and the repository's `coordinate_audit/` directory. The swapped-coordinate control is `maggie_twist_verified_v1.zip`; its strict launcher and exact archive/input hashes are specified in Appendix A.

For the main Z proof, retain `published_fig2_S4.zip` and `published_fig2_closed_path_moves.tar.xz`, including the original and standard-position PD inputs, peripheral comparison records, triangulations and complete 780,847-move log. Appendix A retains the supplied bindings to these archives. The separate direct redraw construction and paths to Z6 also require their original construction/replay archive; this small paper bundle does not replace it.

## Local diagram check

The following check runs using standard Python only:

```sh
python3 corrected_figures_check/audit_corrected_pd.py   corrected_figures_check --output LOCAL_PD_AUDIT_fresh.json
```

It checks exactly the two changed crossings, a planar empty quadrilateral supporting the slide band, and linking/framing consistency in all four supplied PD records. It does not re-extract those records from the figures, generate an automated Kirby-move certificate, or replay a four-dimensional standardization.

With SnapPy installed, `corrected_figures_check/cover_obstruction.py` reproduces the cover comparison of Section 7.2, and `corrected_figures_check/cert_candidate.py` the marked-link comparisons; see `corrected_figures_check/README.md`.

## Remaining author items

Section 7.6 records a direct reduction to the exact six-pentachoron triangulation Z6. To call that route an independently complete second proof, attach its separate standardness certificate/reference or a path to the explicit standard sphere. The main proof of Theorem 1.4 already uses the Z-circle path to the explicit standard sphere and does not depend on this item.

Remark 7.5 leaves the derivation of Figure 1 from the roll-spun construction for author completion. Figure 1 is not a premise of the independent construction proof of Theorem 1.1. Author affiliations and final repository/archival citation decisions remain to be supplied.

## Verification records

`provenance/VALIDATION.json` and `provenance/latex_build.log` are retained **historical Revision 6** records. Their hashes and test flags describe that earlier revision. Current follow-up checks are recorded separately in `provenance/VALIDATION_REV6_1_REVIEW.json` and `provenance/review_rev6_1/`.

### Fresh Revision 6.1 review

With SnapPy 3.3.2, both cover-obstruction calculations and all four positive/negative marked-link comparisons were rerun and reproduced the supplied results. The cusp-to-source-component correspondence was also checked on all six comparison inputs. The restored X0 Sections 2–6, nine isomorphism signatures, four figures and original four PD records are unchanged. The supplied file manifest was checked in full. The newly affected PDF passages were visually reviewed. No change to the paper source, PDFs, figures, PD inputs or checker code was needed.

The reviewed bundle adds these fresh receipts and clarifies historical documentation only. It retains the submitted Revision 6.1 paper unchanged.

### Historical Revision 6 review

The revision review compared the recovered and supplied sources, checked the actual control archive bytes and move counts, and freshly executed the exact local PD audit. It checked arithmetic, signature strings, bibliography claims, compilation and rendering. It did **not** rerun the large X0 or Z construction/reduction certificates, the control's 161,328 moves, or the vertex-link recognition. The paper's assertions of earlier successful replays describe the retained computational record, not fresh runs performed during this edit.

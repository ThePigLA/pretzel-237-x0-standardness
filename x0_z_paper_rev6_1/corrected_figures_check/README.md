# Corrected-diagram checks

This directory contains the original four labelled PD records from combined Revision 5, unchanged. The old/new pairs differ by the two crossing switches used in Lemma 7.1. It also contains the original `CHECK_corrected_figures.txt` summary and original `cert.py` and `orig_retrace.py` files.

## Fresh exact local check

From the manuscript root:

```sh
python3 corrected_figures_check/audit_corrected_pd.py   corrected_figures_check --output LOCAL_PD_AUDIT_fresh.json
```

`LOCAL_PD_AUDIT.json` is the fresh receipt shipped with Revision 6. The script needs standard Python only, uses explicit failures, and has also been executed with `python -O`. It checks the supplied records' planar rotation systems, exact changed crossings, the empty quadrilateral for the band, and the linking/framing arithmetic. It does not extract curves from PDFs, prove the complete local Kirby isotopy by code, or replay any 4D standardization.

## Original SnapPy comparison script

`cert.py` is retained unmodified as historical computational provenance. It prints `CERTIFIED` or `NOT_CERTIFIED`; **its exit status alone is not a pass/fail certificate**. In particular, it can finish with exit status zero after printing `NOT_CERTIFIED`. It requires the recorded SnapPy/Spherogram environment and both comparison inputs. These four PD records alone do not supply all redraw/standard-position comparison inputs used in Proposition 7.3.

`cert_candidate.py` is a fail-closed version of the same comparison that records the accepted cusp matrices and returns a nonzero status on a rejected comparison or an error. In Revision 6.1 its positive branch was executed with SnapPy 3.3.2:

| Comparison | Exit | Status | Receipt |
|---|---|---|---|
| corrected Figure 1 vs first drawing with the two spot crossings switched | 0 | `CERTIFIED_SAME_MARKED_LINK` | `CERT_new_fig1_vs_old_fig1_spot_switched.json` |
| corrected Figure 2 vs first drawing with the two spot crossings switched | 0 | `CERTIFIED_SAME_MARKED_LINK` | `CERT_new_fig2_vs_old_fig2_spot_switched.json` |
| corrected Figure 1 vs first drawing as drawn | 1 | `NOT_CERTIFIED` | `CERT_new_fig1_vs_old_fig1.json` |
| corrected Figure 2 vs first drawing as drawn | 1 | `NOT_CERTIFIED` | `CERT_new_fig2_vs_old_fig2.json` |

`old_fig{1,2}_spot_switched.json` are the earlier archive-order records of the first drawings with the two spot crossings switched; these are the records previously certified equal to the third author's freehand redrawings. The comparison uses the identity cusp correspondence induced by the archive component order, as stated in each receipt. A `NOT_CERTIFIED` result is not by itself a proof of non-equivalence; the cover obstruction below is.

## Cover obstruction (Section 7.2)

```sh
python3 corrected_figures_check/cover_obstruction.py corrected_figures_check/old_fig2.json corrected_figures_check/new_fig2.json --degrees 2,3 --output cover_fig2.json
python3 corrected_figures_check/cover_obstruction.py corrected_figures_check/old_fig1.json corrected_figures_check/new_fig1.json --degrees 2,3 --output cover_fig1.json
```

Requires SnapPy. Each dotted cusp is filled along its zero slope, giving the exterior of the five framed attaching circles in #3(S^1 x S^2). Connected covers are enumerated combinatorially and H_1 is computed exactly. Recorded results: `COVER_OBSTRUCTION_fig2.json` (31 double covers each, same multiset; 572 triple covers each, multisets differ in 63 homology classes) and `COVER_OBSTRUCTION_fig1.json` (31 double covers each, multisets differ in 9 classes). Each run takes a few seconds. This supports the remark in Section 7.2 and is not used in the proof of Theorem 1.4.

## Original extraction script

`orig_retrace.py` is preserved for provenance. It imports historical helper modules through an absolute `/home/claude/motoo_followup/v2` path. Those modules are not included here, so this file is not a standalone PDF-to-PD reproduction tool. A complete extraction distribution should include the helpers and input bindings. No such completion is claimed by this manuscript edit.

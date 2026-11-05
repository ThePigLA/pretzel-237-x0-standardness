# Revision 6.1 (5 November 2026)

Small follow-up to Revision 6. No proof payload, figure or PD record was changed.

- **Typography restored from the accepted X0-only Revision 6.** The two occurrences of "X0-to-G6" in Appendix A are again $X_0$-to-$G_6$, and the list of exact isomorphism signatures after Table 1 again uses the aligned `\makebox` layout (widened for the new labels). All nine signature strings are unchanged.
- **Cover obstruction restored in Section 7.2, now with its record.** Revision 6 had omitted the statement that no homeomorphism of #3(S^1 x S^2) carries one family of attaching circles to the other, because no record was in the bundle. The computation had been made on 3 November; it has now been rerun on the four PD records in this bundle with the new `corrected_figures_check/cover_obstruction.py` (SnapPy 3.3.2). Figure 2: 572 connected triple covers each, H_1 multisets differ. Figure 1: 31 connected double covers each, multisets differ. Receipts: `COVER_OBSTRUCTION_fig{1,2}.json`. Appendix A gives the command. The statement is a remark and is not used in the proof of Theorem 1.4.
- **Positive branch of `cert_candidate.py` executed.** With SnapPy 3.3.2 it certifies corrected Figures 1 and 2 against the first drawings with the two spot crossings switched (exit 0), and rejects the comparison with the first drawings as drawn (exit 1). Receipts and the two comparison inputs are in `corrected_figures_check/`.
- **SnapPy citation checked.** The `is_isometric_to` docstring of SnapPy 3.3.2 states: "The answer True is rigorous, but the answer False may not be as there could be numerical errors resulting in finding an incorrect canonical triangulation." This supports the sentence in the proof of Proposition 7.3.
- Revision labels updated to 6.1. `history/` now also holds the Revision 6 source and PDF; `provenance/revision6_to_revision6_1.diff` is the source diff; `provenance/changes_rev6_to_rev6_1_marked.pdf` marks the changes.

The remaining author items listed below for Revision 6 are unchanged.

---

# Revision 6 review and change log

**Combined X0 / Z manuscript — 5 November 2026**

## Finding

The supplied combined Revision 5 had incorporated the new Z argument into an older version of the X0 text. Comparison with the saved X0-only Revision 6, dated 28 October 2026, confirmed specific lost passages and already accepted presentation fixes. The current revision restores that material while preserving the corrected Tange figures, the new handle-slide argument, the new Z proof, and the current three-author list.

`provenance/MERGE_INPUTS.json` identifies the actual inputs by SHA-256. The attached PDF and the PDF inside the supplied source ZIP were identical. No proof payload or original figure was modified.

## Restored from the accepted X0-only Revision 6

| Location | Restoration |
|---|---|
| Section 1 | The precise statement that the residue class m = 1 mod 4 is not covered by Corollary 1.3; this is not a claim about all current literature. |
| Section 2.1 | Explicit parametrization of the preferred longitude and meridian; u is position along the longitude and v along the meridian. |
| Section 2.1 | The transverse-arc comparison loop has peripheral winding (1,0). |
| Section 2.1 | The covering map q(u,v-tilde) = (u,2 v-tilde), deck translation by 1/2, and the **half-turn** lift of a one-meridian shear, with the exterior deck seam. |
| Section 2.3 | Integral peripheral cocycles, the actual core and normal meridian checks, and all 516 upstairs-edge covering equations. |
| Section 2.3 | The actual 11,136-pentachoron collar audit and transverse-arc/image paths through 6 and 14 pentachora. |
| Sections 3–6 | Accepted presentation fixes and the existing precise trust boundaries for contractions, puncturing, recapping, terminal moves, and smoothing. |
| Acknowledgments | The accepted acknowledgment of substantial OpenAI and Anthropic assistance, with responsibility retained by the authors. |

After removing only the new swapped-coordinate control block, **Sections 2–6 are byte-for-byte identical to the recovered accepted X0-only Revision 6**. This was checked programmatically, not inferred from a visual comparison.

## Added completed swapped-coordinate control

Section 2.4 now records the successful meridional-shear control, and Appendix A identifies its strict verifier and hashes. The recovered `maggie_twist_verified_v1.zip` contains the input, move log and strict receipt.

| Quantity | Value |
|---|---:|
| Initial pentachora | 101,056 |
| Vertex links recognized as three-spheres | 4,226 |
| Full-link-condition edge contractions | 4,203 |
| 3–3 moves | 134,673 |
| 4–2 moves | 18,138 |
| 2–4 moves | 4,297 |
| 5–1 moves | 17 |
| Total moves | 161,328 |
| Endpoint | Boundary of a five-simplex |

The log was recounted and the published hashes were checked against the archive bytes. The retained strict receipt records successful replay and vertex-link recognition. Those large computations were not rerun for this revision.

The wording explicitly separates two conclusions: the peripheral and arc checks identify roll versus twist, while standardness of the swapped cover is a control. Standardness of both covers alone would not determine which coordinate was used. The paper does not infer smooth unknotting from the branch-complement group alone.

The strict reproduction command is `sh verify_certificate_strict.sh`, from the extracted control directory after installing Regina 7.4.1. The original launcher can conceal child failures through a shell pipeline; the strict launcher checks statuses directly and protects assertions. The control must not be run with Python optimization substituted for that launcher. This differs deliberately from the terminal X0 verifier, whose explicit guards remain active under `python -O`.

## Preserved and clarified new Z argument

The manuscript now distinguishes the first drawing's manifold Z-circle from the corrected Figure 2 manifold Z. The new single handle slide transfers the first drawing's standardization to Z. This supersedes the earlier unresolved original-versus-redraw status and has been retained.

The main chain remains:

K-circle -> Z-circle-8 -> W4 <- S2.

It has 780,847 + 423 + 7 = 781,277 recorded moves and reaches the explicit double of a four-simplex. The manuscript's reported replay results and all original hash bindings are retained. These large proof archives were not present in the supplied manuscript ZIP and were not rerun during the revision.

### Local handle-slide audit

The new `audit_corrected_pd.py` was executed on the four supplied records, including with Python optimization enabled. Its guards are explicit failures, not assertions. It confirms:

- Exactly crossing rows 13 and 29 change (zero-based), in both pairs.
- The records give connected planar rotation systems.
- The slide region is an empty quadrilateral bounded by C0, C1, C3 and C2.
- In the supplied orientations, [C0] = (-1,1,0), [C3 old] = (0,1,0), and [C3 new] = (1,0,0), so the slide subtracts C0.
- C0 has zero linking with every other framed component. Its zero framing gives the slid handle framing -1 + 0 - 2(0) = -1.

Lemma 7.1 now explicitly chooses an untwisted band in that region. The local Kirby isotopy remains a mathematical argument. The code certifies its planar-diagram premises, not an automated full handle-move sequence.

### Marked exterior comparison

Proposition 7.3 now leads with the positive marked exterior isometry and requires cusp matrices **+I or -I**, with the same sign on each meridian/longitude pair. That makes the preservation of the nonzero framing explicit. Equality of signatures and cover homology is supplementary; it does not replace the positive map. The official SnapPy documentation supports the distinction between a rigorous positive result and a potentially inconclusive negative numerical search.

Katie is included in Regina 7.4.1, so the supplied release-source build statement was retained. The framing sentence now says that the supplied integers are measured relative to Seifert longitudes, avoiding an unsupported identification with the as-drawn self-writhe. The native-simplification claim is limited to the recorded invocation reproducing the same labelled endpoint in its specified environment.

### Stronger non-equivalence claim omitted

The original Section 7.2 additionally asserted that no homeomorphism of the boundary of the one-handlebody carries the attaching-circle families to each other. That stronger claim needs the comparison of all connected three-fold covers of the **five-cusped attaching-link exteriors after zero-filling the dotted components**. Those inputs and records are not in this source bundle. It is unnecessary for the proof, so the revision retains the directly verified fixed-basis homology obstruction to isotopy and omits the stronger assertion. This is a documentation/scope correction, not a finding that the stronger assertion is false.

## Direct Z6 route: documentation still required

The supplied source asserted that

`gALLQQaacdeffeffqbgaJbCa5ajbwbQbwbwb`

is standard, but immediately contained a TODO to supply its standardness reference or certificate. The revised abstract, introduction, Section 7.6 and Appendix B consistently describe the established direct reduction to this endpoint and the extra evidence needed before calling it an independently complete second proof.

**This does not add a dependency to the main proof of Theorem 1.4**, which goes through Z-circle and the explicit S2 reference. No claim that the Z6 standardness report is false is made.

Remark 7.5 also makes clear that deriving Figure 1 from the roll-spun construction and standardizing its triangulation would be needed for a second proof of Theorem 1.1. The existing X0 construction proof is independent of that diagram.

## Presentation and reproducibility

- Retained Vance Faber, Michael Murphy, Motoo Tange, in that order, and added matching PDF author metadata.
- Restored thin spacing, mathematical subscripts, monospace code identifiers, and uniform en-dashes for move names.
- Preserved every exact isomorphism signature, including the capital O in the G6 signature.
- Removed the duplicate Regina implementation sentence through the accepted text restoration.
- Corrected the Appendix B reference from Proposition 7.2 to 7.3.
- Renumbered the Z move table as Table 4 and the supplement inventory as Table 5, in order of appearance; fixed the extra unused table column.
- Added C5 for coordinate audits and C6 for the completed control.
- Corrected the claim that the manuscript ZIP included the terminal proof archives. It contains manuscript sources, figures and local diagram checks; the proof archives are separate.
- Preserved the original scripts and records. The proposed status-checking `cert_candidate.py` is clearly marked as only partially tested and is not substituted for the original positive comparison record.
- Kept Q and the census out of the paper's exposition.

## Remaining author decisions before final distribution

1. Supply a bound reference/certificate for the separate Z6 standardness claim, or omit the independent-second-proof claim permanently while retaining the direct reduction.
2. Add the short geometric derivation of Figure 1 requested of Motoo, if that diagram is to be identified as X0 independently of the construction proof.
3. Supply affiliations and final repository/DOI/citation details.
4. Ensure the separately distributed construction and move archives accompany a public proof release. In particular retain full PD inputs and peripheral comparison records for Proposition 7.3, not just the small diagnostic summary here.

## Review limits

This is a source merge, mathematical/expository review, exact local diagram audit, archive/count/hash inspection, and PDF quality check. It is not a fresh end-to-end computational verification. Earlier recorded successes in the manuscript remain assertions about the existing proof record. The large X0/Z reductions, the swapped-control replay and vertex-link recognition were not rerun during this edit.

# Published Figure 2 is the standard S⁴ (computed 2026-11-05)

**Confidential: built from Motoo Tange's unpublished diagrams. Do not circulate.**

## Result

The closed 4-manifold presented by **Motoo's published Figure 2** is PL-homeomorphic, and so diffeomorphic, to
the standard 4-sphere. This is the version in which the −1-framed handle passes under the two dotted circles
near the "−1" label.

Together with the earlier chain for the redrawn version (25,200 → 16 → 8 → 6 → known standard), **both
readings of the disputed spot give the standard S⁴**.

## The chain

| # | Step | Object | Check |
|---|---|---|---|
| 1 | Published Figure 2 → Katie-position diagram | 188-crossing PD (`diagram/`) | identical decorated isometry signature; isometry matching framings; identical double-cover homology (`diagram/CERTIFICATION.txt`) |
| 2 | Katie (Regina 7.4.1 `utils/katie.cpp`, compiled from the release source) | ideal, 1,528 pentachora (`katie/`) | 13 vertex links: 12 S³ and one #2(S¹×S²); H₁ = 0 |
| 3 | Staircase truncation of the ideal vertex | finite, 6,112 (`truncation/`) | `PASS_EXPLICIT_FINITE_TRUNCATION`; the boundary is exactly the ideal link (1,528 tetrahedra) |
| 4 | Regina `simplify()` | lower piece, 1,392, boundary 826 tetrahedra (`lower/`) | deterministic: a rerun gives the identical labelled triangulation |
| 5 | Boundary walk 826 → 6 tetrahedra by pure 3D Pachner moves, each realised by attaching one pentachoron | 22,304 layers, giving 23,696 (`boundary/`) | `PASS_LABELLED_3D_REPLAY_AND_EXPLICIT_4D_LAYERING`; collar audit `PASS_ALL_LOCAL_COLLAR_ATTACHMENT_HYPOTHESES` |
| 6 | Glue the archive's ♮2(S¹×B³) cap (474 pentachora, boundary `gLALQccbefefdpcghfc`) | closed, 24,170 (`closed/`) | `PASS_INDEPENDENT_RAW_CLOSED_GLUING`; `PASS_CLOSED_ORIENTABLE_PL_4_MANIFOLD` (19 vertex links all PL S³; f = (19, 12136, 48374, 60425, 24170), χ = 2) |
| 7 | Logged guarded moves (`reduce_closed_sampled.py`, seed 11) | 24,170 → 8, 780,847 moves (`closed_path/`, moves in `published_fig2_closed_path_moves.tar.xz`) | independent replay `PASS_PARTIAL_CLOSED_PATH`; endpoint tight sha 43fdb328… = `iLAALMQkabcccfefhghhhcaaavava2aaaoacbZb6aEa8a8a` |
| 8 | `simplify_penta.py` search | 8 → 6 → 4, then joined to `cPkbbbbaaaaaaaa` = `Example4.sphere()`, 430 moves (`final/`) | `PASS_STANDARD_S4` |

Closing up uses the 474-pentachoron cap, the archive's boundary connected sum of two triangulated S¹×B³.
Any gluing gives the same closed manifold, by Laudenbach–Poénaru.

## Replaying it

```bash
# step 8 (seconds)
python simplify_penta.py verify final/CERTIFICATE.json --expect-start iLAALMQkabcccfefhghhhcaaavava2aaaoacbZb6aEa8a8a
# step 7 (about 23 min with the recognition helper; much slower without it)
tar -xJf published_fig2_closed_path_moves.tar.xz          # gives best8_prefix/
python -O code/helper/run_with_primal_recognition.py --proofs proofs -- \
    <archive>/code/verify_closed_bursts.py best8_prefix --source closed/closed_candidate.rga --out REPLAY.json
# steps 5-6: archive scripts audit_layering_collar.py, audit_fast_closed_gluing.py, audit_raw_algebra4.py --method primal
```

`code/ARCHIVE_SCRIPTS_USED.sha256` lists the scripts from `motoo_figure2_simplification_v1` that were used,
with their hashes.

## Trust boundary

- **Transcription of the published PDF.** Two independent transcriptions agree. One tight clasp hook was read
  the natural way.
- **Regina 7.4.1.** Its checked moves and Katie (compiled from the 7.4.1 source tag).
- **Steps 4 and 6 (native simplify).** Step 4 uses Regina's native `simplify()`, which is deterministic and
  reproduced exactly. An earlier native simplification of the closed triangulation is *not* part of the chain;
  step 7 starts from the glued triangulation itself.
- **Vertex-link recognition.** Steps 6–7 use the archive's recognition helper. It is rebuilt here from
  `primal_recognition_bridge.cpp`, with its engine pin changed to the local Python 3.13 wheel of Regina 7.4.1.
  It certifies 3-sphere and 3-ball links by a replayed trivial-π₁ elimination, using the Poincaré conjecture,
  and falls back to Regina's own recognition otherwise.

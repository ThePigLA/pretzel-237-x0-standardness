# Fresh preparation records

`INTEGRATION_VALIDATION.json` is the concise inventory of the fresh checks
performed on 5 November 2026. `FRESH_COMPOSED_REPLAY.json` binds the regenerated
closed construction to the complete move replay and the explicit standard
anchor.

- `construction_run/`: the full 343.10-second Katie-to-closed build, with
  generated intermediate states, cap and boundary data, logs, checks, and
  the 19-link manifold admission.
- `closed_path_run/`: the fresh 1,256.30-second independent replay of all
  780,847 long-path moves, including the helper controls and all 19 positive
  recognition proofs with their separate replay records.
- `terminal_run/`: the new public terminal command's fresh guarded replay.
- `large_search_run/`: a short search producing a 10,799-move journal and
  an independent replay of that journal. Its endpoint is a partial result;
  no new standardization discovery is claimed.
- `small_search_run/`: fresh small-search and resume smoke checks. The saved
  best chains remained at the initial triangulation, with zero moves.
- `focused_checks/`: structural archive/interface audit, positive marked-PD
  comparison, and deliberately malformed terminal-certificate rejection cases.
- `preflight/`: source-identity, recompression, safe extraction, and altered-XML
  rejection records. Deliberately failing fixtures are identified by their
  receipt labels.

The complete long replay ran once through the unchanged independent historical
entry point. The new full wrapper's extraction, input binding, terminal replay,
and composition functions were separately exercised against those fresh
results. A second duplicate long replay through the new top-level command was
not performed. That distinction is recorded explicitly.

The run receipts preserve the actual paths and hashes used in preparation.
Their manifest hashes describe the frozen input sets at the time of those
runs; the final package manifest additionally includes these receipts and
documentation. The computational source hashes are recorded in the summary.
Historical receipts under `reference/` and `legacy/` are separate from these
fresh execution records.

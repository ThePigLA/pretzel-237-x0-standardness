# Explicit guarded proof for the standard upper cap

The saved nine-pentachoron reference cap is reached from the explicit
24-pentachoron boundary connected sum of two copies of `Delta3 x S1` by
**15 guarded edge collapses**. The endpoint is exactly the existing
labelled reference, including an identical RGA file hash.

This replaces the unlogged `simplify()` step in the historical reference
construction. It does not change that saved reference or its boundary.
It proves a fact about the upper handlebody; it does not establish
standardness of the closed four-manifold.

## Verification

With `regina==7.4.1` installed, run from the project root:

```sh
python -O code/verify_upper_cap_guarded.py \
  --certificate results/figure2_working_20261102/upper_cap_guarded_v1 \
  --target results/boundary/upper_cap_selftest_v1/reference_cap.rga \
  --out upper_cap_replay_fresh.json
```

The independent verifier imports neither the source constructor nor the
search driver. It reconstructs the twelve staircase simplices of each
`Delta3 x S1` product directly from their coordinates, checks that the
connecting boundary tetrahedron is embedded, forms the boundary sum,
and matches the exact saved source. It then checks every collapse guard,
labelled-state hash and size change, and requires exact equality with
the saved target. All input hashes are checked before loading, after
loading and at the end. The completed replay took approximately
0.04 seconds, excluding Python startup.

## Files and bindings

- `constructed_upper24.rga`: explicit geometric source.
- `SOURCE_CONSTRUCTION.json`: product coordinates and connecting facet.
- `moves.jsonl`: all 15 guarded edge collapses.
- `guarded_upper9.rga`: exact endpoint.
- `RESULT.json`: discovery receipt and identity map to the saved target.
- `INDEPENDENT_REPLAY.json`: independent reconstruction and replay.

The target signature is `jafbGbafqajaadiLaNa`. Its file SHA-256 is:

```text
6ff5b588fb2cf655e1a5a40d630312c5fe82da9a2327ea22f2ca5176ab1c0cb1
```

The complete move log SHA-256 is:

```text
154094594ee460da0cc459659a38d09d1f990ea86bab1e76653246b35ea1a66f
```

The remaining software trust is Regina 7.4.1's checked edge-collapse
implementation and public triangulation operations. No unlogged
four-dimensional simplification is used in this proof.

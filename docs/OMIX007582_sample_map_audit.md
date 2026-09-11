# OMIX007582 Sample-Map Audit

## Purpose

This audit records what can be mapped defensibly between the public `OMIX007582`
Mammal40 methylation files shipped in this repository and what remains
unrecoverable without inventing assumptions.

It exists because `OMIX007582` is scientifically useful as an orthogonal
validation asset, but the public release does not currently provide a clean
biological sample key for the beta matrix.

## Files checked

- `data/RAW/data/OMIX007582_beta_matrix.csv`
- `data/RAW/data/OMIX007582-02.csv`
- `data/RAW/data/OMIX007582-03.zip`
- `data/RAW/data/OMIX007582_idat`

## Defensible findings

1. The beta matrix exposes `643` unique technical IDs such as
   `207925070004_R01C01`.
2. The unpacked IDAT directory and the raw `OMIX007582-03.zip` archive contain
   exactly the same `643` technical prefixes represented in the beta matrix.
3. The public metadata file contains `620` biological rows with group and
   tissue labels, not `643`.
4. The candidate identifier columns in the public metadata
   (`OriginalSampleName`, `OriginalSampleName.1`, `sample`) have `0` exact
   overlaps with the beta-matrix technical IDs.
5. The raw archive contains no sample-sheet-like sidecar file that could bridge
   technical Sentrix IDs to biological sample names.

## Conclusion

The only defensible `OMIX007582` map recoverable from the public files bundled
here is:

- `technical ID <-> beta matrix column <-> IDAT prefix`

A defensible biological map of:

- `technical ID <-> metadata sample / tissue / group`

cannot be recovered from the current repository contents alone.

Any order-based assignment of the `620` metadata rows onto the `643` technical
IDs remains exploratory and should not be promoted into the main pipeline,
validation summaries, or interpretation-facing outputs.

## Operational rule for this repo

Use `OMIX007582` in one of two ways only:

1. As a documented non-estimable validation asset when no external sample sheet
   is available.
2. As a biological validation layer only after a separately sourced sample key
   is added and audited.

Do not silently rename beta-matrix columns by metadata order.

## Reproducible audit command

```bash
conda run -n srsc python -m src.omix007582_audit
```

This writes:

- `results/omix007582_sample_map_summary.csv`
- `results/omix007582_candidate_overlap_audit.csv`
- `results/omix007582_technical_id_inventory.csv`
- `results/omix007582_metadata_inventory.csv`

## Current repo implication

The methylation block can stay in the repository as a guarded validation layer,
but its current public-file implementation must keep writing explicit stubs when
it requires a biological sample map that is not defensible from the available
data.

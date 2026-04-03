# SRSC: Reproducible macaque rejuvenation pipeline for tissue vs plasma-exosome signals

SRSC resolves one core problem: **run a transparent, reproducible analysis of macaque multi-omics data to test whether observed rejuvenation signals look more cell-intrinsic or plasma/exosome-associated, with explicit estimability gates when causal linkage is not possible.**

## At a Glance

### What it does
- Builds a **cross-validated transcriptomic aging clock** from primate bulk RNA-seq and derives sample-level rejuvenation proxies (`rejuvenation_score`, `delta_age`).
- Quantifies rejuvenation signal by **group and tissue** with bootstrap summaries and covariate-adjusted tissue expression effects (age/sex/batch when available).
- Integrates plasma proteomics and enforces **strict estimability gates** for mediation/decomposition, plus uncertainty-aware exosome-fraction summaries.

### Phase 0.5 updates
- Adds a conservative plasma-to-animal linkage layer for OMIX007581 sample IDs, with explicit provenance and confidence labels.
- Writes linkage diagnostics as first-class outputs (`plasma_to_animal_map.csv`, `linkage_qc_report.csv`) before estimability gating.
- Runs mediation on **animal-level aggregated rows** (not tissue-replicated rows) to avoid pseudo-replication.

### Run in 3 commands
```bash
conda env create -f environment.yml
conda activate srsc
py -m src.run_pipeline
```

### What you get (main outputs)

| Output | Purpose | Interpretable when unlinked? |
|---|---|---|
| `results/clock_metrics_primates.csv` | CV clock performance (MAE, RMSE, Pearson/Spearman, calibration, CV strategy) | Yes |
| `figures/primates_age_scatter.png` | Chronological vs predicted age sanity check | Yes |
| `results/rejuvenation_by_tissue.csv` | Tissue-level rejuvenation summary (treated vs controls) | Yes |
| `results/tissue_expression_effects.csv` | Covariate-adjusted tissue treatment effects | Yes |
| `results/plasma_biomarkers.csv` | Plasma proteins ranked by association metric | Yes |
| `results/plasma_to_animal_map.csv` | Plasma sample mapping audit (`sample_id -> animal_id`, confidence, rule, validity) | Yes |
| `results/linkage_qc_report.csv` | Mapping coverage/collision diagnostics used before linkage gate | Yes |
| `results/linkage_audit.csv` | Animal-level overlap quality diagnostics across bulk/plasma | Yes |
| `results/estimability_report.csv` | Hard gate (`unlinked`, `partially_linked`, `fully_linked`) for mediation/decomposition | Yes |
| `results/mediation_summary.csv` | Real mediation only when strictly allowed; otherwise explicit structured stub | Conditionally |
| `results/exosome_fraction_summary.csv` | Exosome fraction with bootstrap CI + permutation-based falsification (or explicit non-estimable stub) | Conditionally |
| `results/sensitivity_summary.csv` | Control-set and top-feature-threshold sensitivity runs with direction/stability flags | Yes |

## Scientific objective and claim discipline

### Central question
Can the rejuvenation-associated signal in the macaque study be partitioned into:
- **cell-intrinsic/tissue expression shifts**, vs
- **plasma/exosome-associated signatures**?

### What SRSC currently supports
- Strong, reproducible evidence on **clock behavior** and **group/tissue rejuvenation trends**.
- Covariate-adjusted tissue effects with uncertainty metrics and per-tissue FDR.
- Exosome-fraction estimation with bootstrap uncertainty and permutation falsification when estimable.
- Causal mediation/decomposition only when strict linked-data criteria are met.

### What SRSC explicitly does not over-claim
Without reliable plasma-to-bulk `animal_id` linkage, SRSC does **not** present mediation/decomposition as causal evidence. It writes structured stubs with reasons instead.

## Design Decisions and Scientific Guardrails

- **Guardrailed I/O over convenience**: OMIX loaders enforce metadata sanity checks, alignment checks, and safe fallbacks for messy public files.
- **No silent causal leakage**: clock training uses grouped CV when valid `animal_id` groups exist; otherwise it falls back to KFold and records the CV strategy in metrics.
- **Strict causal gate**: `enable_mediation=True` is not sufficient; decomposition also requires `enable_causal_decomposition=True`, `tier=fully_linked`, and minimum sample criteria.
- **Conservative linkage policy**: deterministic mapping is accepted only for directly compatible plasma codes (`V/WT/GES`) and is validated against observed bulk `orig.ident`.
- **Animal-level mediation inputs**: linked mediation uses one aggregated row per animal to reduce pseudo-replication risk.
- **Machine-consistent summaries**: key CSVs follow a common schema (`available`, `estimable`, `reason`, `n_used`, `method`, `ci_low`, `ci_high`) and include `evidence_level` to reduce over-interpretation.

## Pipeline structure

`src/run_pipeline.py` orchestrates:
1. Load and harmonize bulk RNA-seq metadata/matrix.
2. Train transcriptomic clock and derive rejuvenation proxies.
3. Summarize rejuvenation globally and by tissue.
4. Estimate covariate-adjusted tissue expression effects (with p-values and FDR).
5. Load and clean plasma proteomics, rank plasma biomarkers.
6. Build conservative plasma-to-animal mapping and emit linkage QC diagnostics.
7. Compute linkage audit + estimability tier.
8. Run mediation/decomposition only when strict gate is satisfied; otherwise write structured stubs.
9. Estimate exosome fraction with bootstrap CI + permutation null when estimable.
10. Run sensitivity analyses across control definitions and top-feature thresholds.
11. Write standardized summaries and plots.

## Data scope

SRSC is wired for public OMIX macaque datasets and optional exploratory modules:
- `OMIX007580` (bulk transcriptomics)
- `OMIX007581` (plasma proteomics)
- `OMIX007582` (Mammal40 methylation; optional)
- optional mouse exosome block (translation scaffolding)

Raw full datasets are not redistributed in this repository. Example-compatible execution is supported through the repository layout and config defaults.

## Full-data vs clean/demo usage

- **SRSC-work**: full-data development and debugging environment.
- **SRSC (clean)**: reproducible public-facing run profile for portfolio and collaboration onboarding (typically with reduced example files).

Both use the same pipeline entrypoint and produce the same output schema.
Depending on data reduction in the clean profile, linkage tiers may differ from full-data runs.

## Interpreting current limitation correctly

Public OMIX007581 plasma columns are mixed in linkage quality:
- some are directly compatible with bulk animal IDs (`FV_2`, `MWT_3`, `FGES_1`),
- others remain unresolved without external key metadata (for example, `FY_*`, `MY_*`).

Implication:
- a **high-confidence linked subset** may be analyzable when overlap/quality thresholds are satisfied,
- unresolved subsets stay excluded from linked causal analyses by design,
- if thresholds fail, `estimability_report.csv` and downstream causal outputs remain structured non-estimable stubs.

This is a data linkage constraint, not a software bug.

## Maturity roadmap (next steps)

1. Add robust external mapping to enable true plasma-bulk `animal_id` linkage.
2. Extend exosome-fraction falsification with tissue-aware/per-group constrained null models.
3. Expand cross-species translation module once mouse pathways are fully standardized.
4. Add broader CI coverage for full output schema and sensitivity stability checks.

## Who This Repository Is For

SRSC demonstrates end-to-end ownership across:
- scientific modeling and statistical prudence,
- production-minded data guardrails under messy public omics metadata,
- reproducibility and transparent failure modes,
- interpretable outputs suitable for peer review and collaboration handoff.

## Minimal repository map

- `src/run_pipeline.py`: main entrypoint.
- `src/omix_io.py`: robust matrix/metadata loading and alignment guardrails.
- `src/clocks.py`: transcriptomic clock training and prediction.
- `src/rejuvenation.py`: rejuvenation and tissue-level summaries.
- `src/linkage_audit.py`: bulk-plasma linkage diagnostics and estimability gating.
- `src/viz.py`: plotting utilities.
- `results/`, `figures/`: generated artifacts.

## Reproducibility notes

- Environment: `environment.yml` (conda).
- Determinism: fixed seeds in config and bootstrap paths.
- Main command: `py -m src.run_pipeline`.
- Guardrail tests: `pytest -q tests/test_scientific_guardrails.py`.

For laptop-conservative runs:
```bash
py -m src.run_pipeline --safe
```

# SRSC: Reproducible macaque rejuvenation pipeline for tissue vs plasma-exosome signals

SRSC resolves one core problem: **run a transparent, reproducible analysis of macaque multi-omics data to test whether observed rejuvenation signals look more cell-intrinsic or plasma/exosome-associated, with explicit estimability gates when causal linkage is not possible.**

## At a Glance

### What it does
- Builds a **cross-validated transcriptomic aging clock** from primate bulk RNA-seq and derives sample-level rejuvenation proxies (`rejuvenation_score`, `delta_age`).
- Quantifies rejuvenation signal by **group and tissue** with bootstrap summaries and covariate-adjusted tissue expression effects (age/sex/batch when available).
- Integrates plasma proteomics and enforces **strict estimability gates** for mediation/decomposition, plus uncertainty-aware exosome-fraction summaries.

### Current phase updates
- Adds a conservative plasma-to-animal linkage layer for OMIX007581 sample IDs, with explicit provenance and confidence labels.
- Integrates `OMIX009283` as a mouse exosome mechanism-support block and writes cross-species alignment summaries instead of relying only on a placeholder exosome fraction.
- Promotes methylation to an orthogonal validation layer with explicit multimodal concordance output.
- Keeps mediation on **animal-level aggregated rows** (not tissue-replicated rows) to avoid pseudo-replication.

### Run in 3 commands
```bash
conda env create -f environment.yml
conda activate srsc
python -m src.run_pipeline
```

### Data profiles
- `python -m src.run_pipeline --profile full`: use the full-data working layout under `data/RAW/data`.
- `python -m src.run_pipeline --profile demo`: use the example/demo layout under `data/PROCESSED`.
- `python -m src.run_pipeline --profile auto`: prefer full data when present, otherwise fall back to demo data.
- `python -m src.run_pipeline --data-root <PATH>`: point the selected profile at a different data directory without editing source code.

### What you get (main outputs)

| Output | Purpose | Interpretable when unlinked? |
|---|---|---|
| `results/clock_metrics_primates.csv` | CV clock performance (MAE, RMSE, Pearson/Spearman, calibration, CV strategy) | Yes |
| `figures/primates_age_scatter.png` | Chronological vs predicted age sanity check | Yes |
| `results/rejuvenation_by_tissue.csv` | Tissue-level rejuvenation summary (treated vs controls) | Yes |
| `results/tissue_expression_effects.csv` | Covariate-adjusted tissue treatment effects | Yes |
| `results/mouse_exosome_effects.csv` | Mouse exosome tissue effects from `OMIX009283` using tissue-specific age-clock contrasts | Yes |
| `results/mouse_exosome_signature_summary.csv` | Contrast-level summary of the mouse exosome mechanism-support block | Yes |
| `results/exosome_alignment_by_tissue.csv` | Tissue-wise macaque vs mouse exosome alignment metrics | Yes |
| `results/exosome_alignment_summary.csv` | Primary exosome-aligned contribution summary with explicit evidence tiering | Yes |
| `results/methylation_rejuvenation_by_tissue.csv` | Tissue-level methylation rejuvenation validation summary | Yes |
| `results/multimodal_concordance_summary.csv` | Transcriptome vs methylation concordance summary on overlapping tissues | Yes |
| `results/ovary_subset_validation.csv` | Targeted ovary subset validation audit | Yes |
| `results/hippocampus_subset_validation.csv` | Targeted hippocampus subset validation audit | Yes |
| `results/plasma_biomarkers.csv` | Plasma proteins ranked by association metric | Yes |
| `results/plasma_to_animal_map.csv` | Plasma sample mapping audit (`sample_id -> animal_id`, confidence, rule, validity) | Yes |
| `results/linkage_qc_report.csv` | Mapping coverage/collision diagnostics used before linkage gate | Yes |
| `results/linkage_audit.csv` | Animal-level overlap quality diagnostics across bulk/plasma | Yes |
| `results/estimability_report.csv` | Hard gate (`unlinked`, `partially_linked`, `fully_linked`) for mediation/decomposition | Yes |
| `results/mediation_summary.csv` | Real mediation only when strictly allowed; otherwise explicit structured stub | Conditionally |
| `results/exosome_fraction_summary.csv` | Compatibility exosome-fraction artifact retained while alignment summaries are primary | Conditionally |
| `results/sensitivity_summary.csv` | Control-set and top-feature-threshold sensitivity runs with direction/stability flags | Yes |

## Scientific objective and claim discipline

### Central question
Can the rejuvenation-associated signal in the macaque study be partitioned into:
- **cell-intrinsic/tissue expression shifts**, vs
- **plasma/exosome-associated signatures**?

### What SRSC currently supports
- Strong, reproducible evidence on **clock behavior** and **group/tissue rejuvenation trends**.
- Covariate-adjusted tissue effects with uncertainty metrics and per-tissue FDR.
- Cross-species **exosome-aligned contribution** summaries against `OMIX009283`, with the legacy exosome-fraction CSV kept only as a compatibility layer.
- Methylation-based orthogonal validation and multimodal concordance summaries.
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
- **Canonical group legend**: use `Documents/group_label_crosswalk.md` when translating article labels (`A4-Ctrl`, `A4-WTC`, `A4-SRC`) into repo/OMIX labels (`O_V`, `O_WT`, `O_GES`).

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
9. Build mouse exosome mechanism-support effects and exosome-alignment summaries.
10. Run methylation validation and multimodal concordance checks.
11. Write targeted subset-validation audits for ovary and hippocampus.
12. Estimate the legacy compatibility exosome fraction and run sensitivity analyses.
13. Write standardized summaries and plots.

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

Group naming should also be interpreted through the repository crosswalk:
- `A4-WTC -> O_WT` is strongly supported by the article terminology plus OMIX/BioProject naming.
- `A4-SRC -> O_GES` is strongly supported by the article `SRC` arm plus public BioProject sample names such as `GESMSC-F-1`.
- `A4-Ctrl -> O_V` is strongly supported by the article saline control arm plus public BioProject sample names such as `O-V-F-1`.

Use `Documents/group_label_crosswalk.md` as the operational source of truth for these translations.

Implication:
- a **high-confidence linked subset** may be analyzable when overlap/quality thresholds are satisfied,
- unresolved subsets stay excluded from linked causal analyses by design,
- if thresholds fail, `estimability_report.csv` and downstream causal outputs remain structured non-estimable stubs.

This is a data linkage constraint, not a software bug.

## Maturity roadmap (next steps)

1. Add robust external mapping to enable true plasma-bulk `animal_id` linkage.
2. Extend exosome-fraction falsification with tissue-aware/per-group constrained null models.
3. Extend subset validation from sample-sheet audits to expression-level pseudobulk checks when runtime permits.
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
- Main command: `python -m src.run_pipeline`.
- Guardrail tests: `pytest -q tests/test_scientific_guardrails.py`.
- Methylation sample-map audit: `python -m src.omix007582_audit`.

For laptop-conservative runs:
```bash
python -m src.run_pipeline --safe
```

For the current `OMIX007582` mapping status:
- see `Documents/OMIX007582_sample_map_audit.md`
- regenerate machine-readable audit tables with `python -m src.omix007582_audit`

For low-friction promotion from the full-data work repo into the clean demo repo:
```bash
python -m src.repo_promotion
python -m src.repo_promotion --apply
```
The promotion plan is controlled by `promotion_manifest.json` and copies only approved files into the sibling `SRSC` repo.

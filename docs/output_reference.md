# Output Reference

The pipeline writes machine-readable analysis tables to `results/`, diagnostic
plots to `figures/`, and interpretation-facing plots to `figures/report/`.

## Primary Scientific Outputs

| Output | Purpose | Interpretation boundary |
|---|---|---|
| `clock_metrics_primates.csv` | Cross-validated transcriptomic clock performance | Model validation, not treatment evidence |
| `rejuvenation_by_tissue.csv` | Tissue-level treated-versus-control `delta_age` summaries | Use confidence intervals; ranking is not confirmation |
| `tissue_expression_effects.csv` | Covariate-adjusted tissue expression effects | Depends on released covariate quality |
| `plasma_biomarkers.csv` | Protein ranking with direction and stability fields | Small-cohort association only |
| `plasma_age_axis_summary.csv` | Oriented plasma PC1 group summary | Age-state axis, not a biological clock or causal mediator |
| `plasma_age_axis_delta_age_correlation.csv` | Linked plasma-axis versus tissue `delta_age` association | Conditional on animal linkage and sample size |
| `exosome_alignment_by_tissue.csv` | Tissue-wise macaque/mouse alignment | Cross-species mechanism support only |
| `exosome_alignment_summary.csv` | Contrast-level alignment summary | Primary exosome-alignment interface |
| `mediation_summary.csv` | Animal-level mediation estimates or structured stub | Causal interpretation requires stable estimates and assumptions |
| `multimodal_concordance_summary.csv` | Transcriptome/methylation validation status | Blocked when biological sample mapping is unavailable |

`clock_metrics_primates.csv` records the fitted feature count, numerical input
dtype, Ridge alpha, and solver. The pipeline uses float64 model inputs to avoid
single-precision conditioning warnings; this is a numerical-stability control,
not a change in feature scaling or model class.

For estimable mediation, `Total`, `ADE`, and `ACME` are the canonical effect
fields and their `*_CI` columns contain effect-specific bootstrap intervals.
The common `ci_low` and `ci_high` fields repeat the `Total_CI` bounds so the
table retains the repository-wide result contract. ADE and ACME must not be
relabeled as cellular and exosome-causal effects without stronger identifying
evidence.

## Linkage and Estimability Diagnostics

| Output | Purpose |
|---|---|
| `plasma_to_animal_map.csv` | Mapping provenance, confidence, and validity per plasma sample |
| `linkage_qc_report.csv` | Coverage and collision diagnostics |
| `linkage_audit.csv` | Bulk/plasma animal-overlap diagnostics |
| `estimability_report.csv` | Gate for mediation and linked decomposition |

These are required scientific diagnostics. A mechanism-facing result should not
be interpreted without its linkage and estimability context.

## Supporting and Compatibility Outputs

| Output | Purpose |
|---|---|
| `mouse_exosome_effects.csv` | Tissue-level `OMIX009283` effects |
| `mouse_exosome_signature_summary.csv` | Mouse contrast summaries |
| `methylation_rejuvenation_by_tissue.csv` | Methylation validation or structured limitation |
| `ovary_subset_validation.csv` | Targeted ovary subset audit |
| `hippocampus_subset_validation.csv` | Targeted hippocampus subset audit |
| `sensitivity_summary.csv` | Control-set and feature-threshold sensitivity |
| `exosome_fraction_summary.csv` | Legacy compatibility artifact; not the primary attribution result |
| `omix007582_sample_map_summary.csv` | Mammal40 mapping audit |
| `omix009284_audit_summary.csv` | PBMC single-cell structural audit |

## Standard Result Contract

Interpretation-facing tables use the following fields where applicable:

- `available`: whether the required source data were available.
- `estimable`: whether the configured analysis could be evaluated.
- `reason` and `reason_code`: human- and machine-readable status.
- `missing_author_key`: exact external metadata required when blocked.
- `n_used`: observations used by the analysis.
- `method`: compact method identifier.
- `ci_low` and `ci_high`: uncertainty interval where defined.
- `evidence_level`: design/evidence rung, not a statistical significance score.

## Report Figures

Run `python -m src.report_figures` after the pipeline. The generated
`figures/report/report_figure_manifest.csv` records each figure's source tables,
status, and interpretation note. The report set includes:

- tissue rejuvenation and tissue-priority views;
- exosome-alignment summary and tissue drivers;
- estimability, evidence-ladder, and public-data-ceiling views;
- mediation uncertainty;
- plasma biomarker stability and heuristic categories;
- the oriented plasma aging axis;
- analysis sensitivity summaries.

Three report figures are curated separately for the public README:

| Asset | Source tables | Interpretation boundary |
|---|---|---|
| `docs/assets/aging_rejuvenation_signal.png` | `clock_metrics_primates.csv`; `rejuvenation_by_tissue.csv` | Grouped clock validation plus nominal tissue prioritization; not confirmed tissue effects |
| `docs/assets/multimodal_evidence_architecture.png` | Rejuvenation, linkage, alignment, and multimodal status tables | Status-aware evidence map; not a pooled multimodal effect |
| `docs/assets/estimability_guardrail.png` | `linkage_qc_report.csv`; `estimability_report.csv`; `mediation_summary.csv` | A passed gate means estimable, not causally established |

These tracked images are release snapshots. Regenerate them from reviewed
full-profile results with:

```bash
python -m src.report_figures --portfolio-assets-dir docs/assets
```

Protein categories are heuristic symbol groupings unless a dedicated annotation
field is supplied; they are not formal pathway-enrichment results.

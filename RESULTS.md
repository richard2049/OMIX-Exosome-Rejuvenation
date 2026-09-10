# Results: What the Public OMIX Data Currently Support

This document is the scientific reading guide for the current full-data output
tables. It separates computational estimability from statistical stability and
biological interpretation. Values below reflect the latest local run and should
be regenerated before formal citation.

## Executive Interpretation

The analysis recovers a strong transcriptomic age signal and a usable linked
plasma subset. It does not yet provide stable evidence that macaque
rejuvenation is mediated primarily by plasma or exosomes. Tissue effects are
uncertain, linked mediation confidence intervals cross zero, cross-species
alignment is based on only four shared tissues, and methylation validation is
blocked by a missing biological sample map.

The current contribution is therefore methodological and hypothesis-generating:
the workflow identifies which claims are estimable, quantifies their
uncertainty, and names the metadata required to move beyond the public-data
ceiling.

## 1. Transcriptomic Aging Signal

The Ridge transcriptomic clock was evaluated with five-fold
`GroupKFold(animal_id)` cross-validation across 2,058 tissue samples from 61
animals.

| Metric | Current value | Interpretation |
|---|---:|---|
| MAE | 2.65 years | Mean absolute cross-validated prediction error |
| RMSE | 3.24 years | Error with greater weight on large deviations |
| Pearson correlation | 0.86 | Strong linear age association |
| Spearman correlation | 0.89 | Strong rank-order age association |
| Calibration slope | 0.65 | Predictions are compressed toward the mean |

This supports use of the clock as an age-associated signal model. It does not,
by itself, prove biological rejuvenation or establish causal treatment effects.

Source: `results/clock_metrics_primates.csv`.

## 2. Tissue Rejuvenation

The tissue analysis estimates treated-versus-control differences in
cross-validated `delta_age`. Thirty-nine tissues are estimable, but every
bootstrap confidence interval crosses zero. The apparent younger- and
older-shifted tissues should therefore be treated as candidates for follow-up,
not as confirmed tissue-specific responses.

Differences in tissue sample size are represented through uncertainty and the
reported `n_ctrl` and `n_trt`; effects are not rescaled by sample size. The
`signal_to_uncertainty` field is a prioritization aid, not a second effect size.

Sources: `results/rejuvenation_by_tissue.csv` and
`figures/report/report_tissue_rejuvenation_forest.png`.

## 3. Plasma State and Linkage

The deterministic linkage audit maps 24 of 32 plasma samples to animals present
in the bulk transcriptomic cohort, with no detected mapping collisions. Eight
plasma samples remain unresolved and are excluded from linked analyses.

The oriented plasma PC1 explains 30.8% of the selected protein variance. It is
oriented so that higher values mean older-like only because the old-control
median exceeds the young median. The confidence interval for that old-young
gap crosses zero, so the orientation denominator and any fraction-of-gap
summary are unstable.

In the current run, the GES group lies above the old-control median on this
axis rather than shifting toward the young group. This is a result worth
falsifying, not evidence that the treatment accelerates aging: PC1 is an
unsupervised plasma state axis, the cohort contains 32 samples, and its linked
correlation with tissue `delta_age` is weak and uncertain (`rho = 0.17`, 95%
CI `-0.24` to `0.56`, permutation `p = 0.448`, `n = 24`).

Sources: `results/linkage_qc_report.csv`,
`results/plasma_age_axis_summary.csv`, and
`results/plasma_age_axis_delta_age_correlation.csv`.

## 4. Linked Mediation

The estimability gate permits animal-level mediation on 24 linked animals: 8
treated and 16 controls. This is an important design achievement, but the
estimates are not stable enough for causal partitioning.

| Quantity | Estimate | Bootstrap confidence interval |
|---|---:|---:|
| Indirect effect (ACME) | -0.59 | -2.35 to 1.31 |
| Direct effect (ADE) | 0.50 | -2.04 to 2.61 |
| Total effect | -0.09 | -1.21 to 1.03 |
| Proportion mediated | 6.44 | -30.22 to 20.67 |

All intervals cross zero, and the proportion-mediated estimate is unstable
because the total effect is near zero. Level 4 here means that linked mediation
was computationally estimable under the stated assumptions and configured gate;
it does not mean the resulting causal decomposition is supported.

Source: `results/mediation_summary.csv`.

## 5. Cross-Species Exosome Alignment

The mouse `OMIX009283` support block shares four tissues with the macaque
analysis. Both GES-versus-vehicle and WT-versus-vehicle contrasts are signed
concordant in three of four tissues. However, tissue rank correlations are
negative (`rho = -0.80` and `-0.40`) and permutation tests are not significant
(`p = 0.528` and `0.661`). Hippocampus is directionally discordant in both
contrasts.

The scientifically defensible reading is limited directional overlap with weak
rank/magnitude agreement. Cross-species alignment can support or challenge a
mechanistic hypothesis, but it cannot substitute for exosome-linked
measurements in the treated macaques.

Sources: `results/exosome_alignment_summary.csv` and
`results/exosome_alignment_by_tissue.csv`.

## 6. Orthogonal Validation

Biological methylation validation is not currently estimable. The released
Mammal40 beta matrix uses technical Sentrix identifiers, while the public files
available to this project do not provide the required mapping to biological
sample, tissue, group, sex, and age. The workflow records this as
`OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING` instead of applying an order-based
join.

`OMIX009284` is auditable as PBMC single-cell RNA data, but its PBMC-only scope
does not resolve cross-tissue exosome attribution. It remains a candidate for a
separate immune-cell-state analysis.

Sources: `results/multimodal_concordance_summary.csv`,
`Documents/OMIX007582_sample_map_audit.md`, and
`Documents/OMIX009284_audit.md`.

## Bottom Line

The public evidence does not distinguish a predominantly exosome-mediated
effect from a predominantly tissue-intrinsic effect with reliable precision.
The strongest next scientific upgrade is not a more complex decomposition
model. It is recovery of the missing sample maps and direct exosome
cargo/donor/recipient linkage described in
`Documents/public_data_ceiling.md`.

Until those data are available, tissue rankings, plasma associations, and
cross-species alignment should remain uncertainty-aware prioritization tools.

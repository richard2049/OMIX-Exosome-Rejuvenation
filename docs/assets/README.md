# Curated README Figures

This directory contains exactly three interpretation-facing figure snapshots
for the public README. They are generated from reviewed full-profile result
tables, not assembled manually and not selected from the diagnostic plot set.

| Figure | Source tables | Scientific boundary |
|---|---|---|
| `aging_rejuvenation_signal.png` | `clock_metrics_primates.csv`; `rejuvenation_by_tissue.csv` | Shows grouped clock validation and uncertainty-aware tissue prioritization; it does not establish tissue-specific rejuvenation |
| `multimodal_evidence_architecture.png` | `rejuvenation_by_tissue.csv`; `linkage_qc_report.csv`; `exosome_alignment_summary.csv`; `multimodal_concordance_summary.csv` | Shows which evidence streams are available or blocked; it is not a pooled concordance estimate |
| `estimability_guardrail.png` | `linkage_qc_report.csv`; `estimability_report.csv`; `mediation_summary.csv` | Shows whether mediation may be estimated and how failure is reported; gate passage does not establish causality |

Regenerate both the complete report layer and these public snapshots with:

```bash
python -m src.report_figures --portfolio-assets-dir docs/assets
```

Only run the publication step after reviewing full-profile result tables and
the generated `figures/report/report_figure_manifest.csv`. The remaining report
and diagnostic figures stay generated and are not part of the README image set.

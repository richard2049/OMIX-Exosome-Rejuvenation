# OMIX Exosome Rejuvenation: Primate-First Attribution Plan

## Summary

- Replace the current narrative roadmap with a phased implementation plan whose primary endpoint remains reproduced macaque rejuvenation, while attribution outputs are framed as evidence-weighted exosome alignment rather than causal partition.
- Use a hybrid roadmap: Phase 1 is implementation-ready and should be executed next; later phases are explicitly deferred and gated by earlier validation.
- Preserve the current `results/exosome_fraction_summary.csv` contract during transition, but make new alignment outputs the primary attribution interface.

## Strategic Objective

- The repository is not intended only to reproduce the SRSC paper. Its purpose is to test whether the public data can support new, defensible insight about therapeutic promise, limitations, and mechanism for the SRSC-exosome and FOXO3A-enrichment strategy.
- The preferred endpoint remains a separation of exosome-aligned and cell-intrinsic rejuvenation signals, but that endpoint must be downgraded to evidence-weighted alignment when public data do not support causal decomposition.
- Treat the public-data ceiling as an explicit decision point. If all reproducible public-data analyses converge on linkage, sample-map, or design limitations, the correct next step is to request missing sample maps or design details from the authors rather than extending speculative analyses.
- Use `Documents/scientific_objectives_and_decision_framework.md` to classify future work as reproduction, mechanism support, therapeutic insight, hypothesis generation, or unsupported speculation.
- Use `Documents/public_data_ceiling_and_author_request.md` as the operational boundary between defensible public-data analyses, hypothesis-generating extensions, and author-side metadata requests.

## Implementation Changes

### Phase 1: Mouse exosome mechanism support and alignment

- Integrate `OMIX009283` as a separate mouse exosome mechanism module driven by parsed sample names, with deterministic fields `age`, `sex`, `tissue`, and `arm`.
- Normalize the mouse intervention arms to `Veh`, `WT`, and `GES`; compute per-tissue effects for `GES vs Veh`, `WT vs Veh`, and `GES vs WT`.
- Emit `results/mouse_exosome_effects.csv` and `results/mouse_exosome_signature_summary.csv` using the repo's standard result schema where applicable.
- Replace the current attribution target with two primary outputs: `results/exosome_alignment_by_tissue.csv` and `results/exosome_alignment_summary.csv`.
- Per tissue, compute signed concordance, Spearman rank correlation, standardized effect similarity, bootstrap confidence intervals, permutation p-values, and explicit `available`, `estimable`, `reason`, `reason_code`, `missing_author_key`, `n_used`, `method`, `ci_low`, `ci_high`, and `evidence_level` fields.
- Keep `results/exosome_fraction_summary.csv` in Phase 1 as a compatibility artifact only. It should continue to be written, but downstream interpretation and documentation should point to the alignment summaries as the primary attribution outputs.
- Do not merge mouse exosome evidence into macaque mediation outputs; cross-species evidence can strengthen support tiering, but it cannot be presented as direct primate causal mediation.

### Phase 2: Orthogonal validation

- Promote `OMIX007582` from optional side module to validation-only support for macaque rejuvenation, not mechanism discovery.
- Emit `results/methylation_rejuvenation_by_tissue.csv` and `results/multimodal_concordance_summary.csv`, focused on DNAmAge rescue and transcriptome-methylation directional concordance by tissue.
- Use `OMIX007583` and `OMIX007586` as targeted validation sets only, producing `results/ovary_subset_validation.csv` and `results/hippocampus_subset_validation.csv`.
- Require these subset outputs to preserve intervention-era controls separately from stage controls; never collapse `O_C` and `O_V`.

### Phase 3: Evidence system and deferred audit

- Expand repo-wide `evidence_level` semantics to a fixed ladder:
  - `0`: not estimable
  - `1`: macaque rejuvenation reproduced
  - `2`: macaque rejuvenation plus linkage-supported plasma association
  - `3`: macaque rejuvenation plus orthogonal exosome-alignment support
  - `4`: linked mediation estimable under stated assumptions
- Apply that ladder consistently to all interpretation-facing CSVs and report text.
- Treat the `OMIX009284` read-only audit as complete. It identifies a PBMC single-cell RNA export with recoverable suffix-to-sample mapping, but PBMC-only scope means it should not be integrated into the current cross-tissue attribution path unless a dedicated PBMC cell-state question is added.

## Important Interfaces

- Keep `PipelineConfig` as the single control surface. Add explicit toggles and thresholds for the mouse exosome alignment block, subset validation block, and minimum shared-tissue requirements instead of hardcoding them inside analysis functions.
- Keep the standard machine-readable CSV contract consistent across new outputs: `available`, `estimable`, `reason`, `reason_code`, `missing_author_key`, `n_used`, `method`, `ci_low`, `ci_high`, and `evidence_level`.
- Update public-facing docs and future report language to prefer "exosome-aligned contribution" and "mechanism support"; avoid "causal decomposition" or "proves exosomes are the main cause" unless Level 4 criteria are actually met.

## Test Plan

- Parser and metadata tests: `OMIX009283` sample parsing must recover `age`, `sex`, `tissue`, and `arm` deterministically, and arm normalization must preserve `Veh`, `WT`, and `GES` as distinct groups.
- Guardrail tests: unlinked primate runs must still emit structured stubs for mediation; cross-species alignment outputs must downgrade cleanly when tissue overlap is insufficient; mouse evidence must never be promoted to direct macaque causal evidence.
- Validation tests: methylation outputs must still write explicit summaries under partial overlap; ovary and hippocampus subset validations must preserve `O_C` versus `O_V` distinctions and report disagreement explicitly.
- Smoke test: `python -m src.run_pipeline` must regenerate the new result files alongside the existing core outputs, and `pytest -q tests/test_scientific_guardrails.py` must cover the evidence-tier and alignment downgrade behavior.

## Assumptions and Defaults

- The plan is internal and implementation-facing, not public narrative documentation.
- Primary scientific question: how much macaque rejuvenation is consistent with an exosome-derived mechanism, and how much remains non-exosome-aligned.
- `OMIX009283` is the next implementation priority and is treated as mechanistic support only.
- `OMIX007582`, `OMIX007583`, and `OMIX007586` are validation assets, not new primary cohorts.
- `OMIX009284` is deferred after audit; use it only for a future PBMC-focused mechanistic module.

# Multimodal Validation

## Purpose

The project uses multiple modalities to ask whether a result is reproducible,
orthogonally consistent, or compatible with a proposed mechanism. It does not
pool heterogeneous assays into a single evidence score or treat agreement as
proof of causality.

## Evidence Layers

| Layer | Validation question | Role in this repository |
|---|---|---|
| Macaque transcriptomics | Are age relationships reproducible, and do treated tissues shift toward a younger-like transcriptomic state? | Primary analysis |
| Macaque plasma proteomics | Is a circulating state associated with age or intervention? | Exploratory association; linked analyses require animal identity |
| Macaque methylation | Is transcriptomic rejuvenation directionally supported by an independent epigenetic modality? | Orthogonal validation, currently blocked when the sample map is unavailable |
| Ovary and hippocampus subsets | Are selected tissue findings consistent in targeted released subsets? | Validation-only support; controls remain design-specific |
| Mouse exosome perturbation | Are macaque tissue effects directionally aligned with a controlled exosome-related perturbation? | Cross-species mechanism support or contradiction, not direct macaque mediation |

## Integration Rules

1. Preserve assay-specific effect directions, uncertainty, sample counts, and
   reasons for non-estimability.
2. Compare compatible tissues and explicitly defined contrasts rather than
   merging all measurements into one pooled endpoint.
3. Treat discordance as informative; it may narrow a mechanism rather than
   represent a technical failure.
4. Keep cross-species and unlinked evidence below direct causal attribution.
5. Report blocked modalities explicitly so that missing validation is not
   mistaken for agreement.

The principal machine-readable summaries are
`multimodal_concordance_summary.csv`, `exosome_alignment_summary.csv`, and the
dataset-specific validation tables described in the
[output reference](outputs.md). Current findings and their uncertainty
are summarized in [RESULTS](../RESULTS.md).

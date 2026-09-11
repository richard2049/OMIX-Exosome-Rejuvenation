# Public Data Limitations and Evidence Ceiling

## Purpose

This page defines which conclusions can be supported by the public OMIX
data used in this repository and which require additional metadata or new
experimental evidence. Its purpose is to prevent additional model complexity
from being mistaken for stronger biological identification.

## Defensible Uses of the Public Data

The guarded workflow can support:

- evaluation of a macaque transcriptomic age model and uncertainty-aware
  tissue treatment effects;
- covariate-aware tissue expression analyses, subject to released metadata;
- conservative plasma biomarker ranking and plasma-to-animal linkage audits;
- animal-level mediation as an exploratory analysis when its linkage and
  sample-support gate passes;
- mouse `OMIX009283` comparisons as cross-species, non-causal mechanism
  support or contradiction;
- technical audits of `OMIX007582` and `OMIX009284` with explicit structured
  limitations.

These outputs support reproduction, falsification, and hypothesis
prioritization. They do not establish that exosomes caused the macaque tissue
response.

## Claims Not Established

The current public study design does not establish:

- a numeric causal partition into exosome-derived and tissue-intrinsic effects;
- exosomes as the principal cause of macaque rejuvenation;
- causal FOXO3A necessity or sufficiency;
- biological methylation concordance from `OMIX007582`, because the Sentrix
  technical identifiers lack a defensible biological sample map;
- a cross-tissue mechanism from the PBMC-only `OMIX009284` release.

## Metadata Required for Stronger Inference

The following author-side records would materially increase evidentiary
strength:

- a complete `OMIX007581` plasma proteomics
  `sample_id -> animal_id -> group -> sex -> age` mapping;
- confirmation of the `OMIX007580` bulk RNA-seq
  `sample_id -> animal_id -> tissue -> group -> sex -> age` mapping;
- exosome cargo, preparation, donor, batch, dose, treatment, recipient, and
  outcome linkage;
- the `OMIX007582` Sentrix barcode and position to biological sample, tissue,
  group, sex, and age mapping;
- assay-specific QC and exclusion tables;
- direct exosome uptake, biodistribution, or cargo-transfer measurements linked
  to recipient animals and tissues.

## Evidence Upgrade Rules

- Level 1 records an estimable macaque rejuvenation analysis.
- Level 2 adds linkage-supported plasma association.
- Level 3 adds an estimable orthogonal exosome-alignment analysis.
- Level 4 records linked mediation estimable under stated animal-level
  assumptions.

These levels describe design requirements reached, not effect significance or
causal certainty. In particular, Level 4 does not establish a stable causal
partition when mediation confidence intervals cross zero or identifying
assumptions remain unsupported.

## Operational Rule

When an analysis depends on unavailable information, the pipeline must emit a
structured non-estimable row with a stable `reason_code` and the exact external
requirement in `missing_author_key`, or label the result exploratory and state
the limitation. Missing biological identity must not be replaced by positional
joins, row order, group matching, or cross-species substitution.

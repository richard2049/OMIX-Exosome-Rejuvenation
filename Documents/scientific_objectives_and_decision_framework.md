# Scientific Objectives and Decision Framework

## Purpose

This repository is not only a reproduction of the SRSC paper. Its primary
purpose is to test whether the public data can generate new, defensible insight
about the therapeutic promise, limitations, and likely mechanism of the
SRSC-exosome and FOXO3A-enrichment strategy.

Replication remains necessary because it establishes that the local pipeline,
metadata harmonization, and effect summaries behave consistently. It is not the
final scientific endpoint.

## Central Scientific Question

The preferred question is:

> How much of the observed macaque rejuvenation signal is consistent with an
> exosome-associated mechanism, and how much appears non-exosome-aligned or
> tissue-intrinsic?

With the current public data, this should be implemented as evidence-weighted
alignment and falsification rather than direct causal decomposition. A numeric
exosome-versus-cellular partition requires stronger linked sample metadata than
the public files currently provide.

## Objective Ladder

### Objective 1: Reproduce the macaque rejuvenation signal

This is the baseline requirement. The pipeline should continue to reproduce
aging-clock behavior, tissue-level rejuvenation summaries, covariate-aware
tissue expression effects, and sensitivity diagnostics.

Interpretation: reproduced study signal.

Not sufficient for: mechanism claims, exosome causality, or therapeutic
prioritization.

### Objective 2: Identify robust tissue and pathway patterns

The next useful contribution is to determine which tissues, expression programs,
and biomarkers show the most reproducible treatment-associated changes and
which results are fragile to model choices, control definitions, or sample
linkage.

Interpretation: therapeutic signal prioritization.

Useful outputs include:

- tissue ranking by rejuvenation magnitude and uncertainty;
- tissue-specific discordance or weak-response flags;
- plasma biomarker ranking with linkage confidence;
- sensitivity summaries for control set and feature-threshold choices.

### Objective 3: Estimate exosome-aligned mechanism support

The repository should compare macaque rejuvenation-associated effects with
independent exosome-relevant evidence, especially `OMIX009283`, without
claiming direct primate mediation.

Interpretation: mechanism support or mechanism inconsistency.

Preferred language:

- "exosome-aligned contribution";
- "mechanism support";
- "consistent with an exosome-associated mechanism";
- "non-exosome-aligned residual signal".

Avoid unless Level 4 evidence is available:

- "causal decomposition";
- "the exosome fraction is X percent";
- "exosomes are the main cause";
- "purely cellular effect".

### Objective 4: Test orthogonal validation and falsification

Public methylation, subset, and single-cell resources should be used to test
whether the transcriptomic conclusions are directionally supported or
contradicted. These modules should not be promoted into primary mechanism
claims unless their sample mapping and design support that use.

Interpretation: validation, contradiction, or hypothesis generation.

Current boundaries:

- `OMIX007582` supports technical-ID audit and bounded Mammal40 rebuild tests
  until a biological sample map is available.
- `OMIX007583` and `OMIX007586` should remain targeted subset-validation
  assets.
- `OMIX009284` is PBMC-only single-cell RNA data and is best suited for a
  future immune-cell-state question, not cross-tissue exosome attribution.

### Objective 5: Assess therapeutic promise and risk

The most useful downstream contribution is not a single rejuvenation score. It
is a balanced assessment of where the strategy looks promising, where it is
weak, and which mechanistic or safety questions remain unresolved.

Interpretation: therapeutic hypothesis prioritization.

Useful analyses include:

- tissue-specific response strength and uncertainty;
- evidence for broad versus tissue-restricted benefit;
- disagreement between transcriptomic, methylation, plasma, and single-cell
  signals;
- pathway-level evidence for beneficial repair versus stress, inflammation, or
  compensatory responses;
- identification of the smallest author metadata set needed to validate or
  reject the strongest hypotheses.

## Claim Classes

Future outputs and documentation should classify claims into one of these
classes:

- Reproduced result: directly supported by public data and guarded pipeline
  outputs.
- Mechanism support: consistent with exosome or FOXO3A-related biology, but
  not direct causality.
- Therapeutic insight: useful for prioritization, risk assessment, or
  hypothesis selection.
- Hypothesis-generating result: plausible but dependent on missing linkage,
  incomplete metadata, or cross-species assumptions.
- Not supported: would require private metadata, new experiments, or direct
  intervention-linked measurements.

## Public-Data Ceiling

The repository has reached the practical public-data ceiling for causal
attribution when all available public datasets have been audited or integrated
and the remaining question still depends on one of these missing keys:

- plasma proteomics sample-to-animal mapping;
- bulk RNA-seq sample-to-animal confirmation;
- exosome cargo, donor, preparation, batch, dose, and recipient linkage;
- `OMIX007582` Sentrix technical ID to biological sample mapping;
- direct exosome uptake, biodistribution, or cargo-transfer measurements linked
  to recipient animals and tissues.

At that point, adding more model complexity would mainly increase apparent
precision without increasing scientific validity. The correct next step is to
request author metadata or design a new validation experiment.

## Operational Decision Rule

For each new analysis, ask:

1. Does it strengthen reproduction, mechanism support, therapeutic
   prioritization, falsification, or metadata triage?
2. Can it be run from public data without positional or order-based assumptions?
3. Does it produce a machine-readable output with uncertainty, sample counts,
   an explicit evidence level, a stable `reason_code`, and a
   `missing_author_key` when the result depends on author-side metadata?
4. Would a negative or discordant result be interpretable?
5. Does it avoid upgrading cross-species or unlinked evidence into direct
   primate causal mediation?

If the answer to any of these is no, keep the analysis exploratory or defer it.

## Current Evidence Boundary

The remaining evidence upgrades and their required metadata are documented in
`Documents/public_data_ceiling.md`. Until those requirements are met, new
mechanism-facing analyses should remain explicitly exploratory.

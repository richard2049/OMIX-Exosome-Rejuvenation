# OMIX Scientific Inference and Estimability Contract

## Purpose

This document defines the evidence requirements that govern linkage-dependent,
mediation, exosome-related, and cross-species claims in the OMIX primate
rejuvenation project.

Its purpose is to distinguish:

1. technically computable quantities;
2. statistically estimable quantities;
3. biologically supported interpretations;
4. mechanistic or causal claims.

A result must satisfy the requirements of each level before it is promoted to
the next.

---

## Evidence States

Use the following project-wide evidence vocabulary:

- **Observed** — directly estimated from the relevant dataset under the
  documented design.
- **Supported** — consistent with additional orthogonal evidence but not
  independently causal.
- **Exploratory** — useful for hypothesis generation but not strong enough for
  primary interpretation.
- **Not estimable** — the required linkage, sample size, design, or information
  is unavailable.
- **Not established** — the analysis can be performed, but the requested
  mechanistic or causal conclusion does not follow from it.

---

## 1. Linkage Contract

### Question
Can measurements from different modalities be assigned to the same biological
animal with sufficient confidence for individual-level inference?

### Required evidence

A valid link must have:

- explicit sample identifiers;
- documented mapping provenance;
- mapping rule;
- confidence level;
- collision/duplication checks;
- overlap with the relevant modality;
- biological plausibility;
- no unresolved many-to-one or one-to-many ambiguity unless explicitly modeled.

### Forbidden shortcuts

Do not infer linkage solely from:

- column order;
- similar sample counts;
- treatment group;
- sex;
- tissue;
- approximate naming patterns;
- expected experimental design.

### Outputs

Canonical artifacts:

- `plasma_to_animal_map.csv`
- `linkage_qc_report.csv`
- `linkage_audit.csv`

Unresolved samples remain unresolved.

---

## 2. Estimability Contract

### Principle

Technical computability does not imply scientific estimability.

An analysis is estimable only when its required biological unit, linkage,
sample support, and model assumptions are satisfied.

### States

- `unlinked`
- `partially_linked`
- `fully_linked`

The exact state must be produced from documented criteria rather than analyst
judgment after observing results.

### Failure behavior

When requirements fail:

- do not calculate an interpretable causal estimate;
- do not substitute group averages for individual linkage;
- do not silently weaken the model;
- emit a structured non-estimable result;
- document the failed criterion.

Canonical artifact:

- `estimability_report.csv`

---

## 3. Mediation Contract

### Required causal structure

A mediation analysis conceptually requires:

Treatment → Mediator → Outcome

with the mediator and outcome measured at compatible biological units.

### Minimum requirements

- valid treatment assignment;
- valid animal-level linkage;
- mediator measured for the relevant animal;
- outcome measured for the relevant animal;
- sufficient linked sample support;
- no tissue-level pseudo-replication of one animal-level mediator;
- model assumptions documented;
- configured mediation and causal-decomposition gates both satisfied.

### Interpretation

Allowed:

> "The linked subset is compatible with an indirect association through X."

Not automatically allowed:

> "X mediates the rejuvenation effect."

Forbidden without stronger identification:

> "Exosomes caused X% of the rejuvenation."

A statistically non-zero mediation estimate does not independently establish
the causal assumptions needed for mechanistic attribution.

---

## 4. Exosome-Causality Contract

Distinguish four levels:

### Level A — Exosome-associated
A feature is observed in plasma/exosome-related measurements.

### Level B — Exosome-aligned
The direction or molecular pattern is concordant with an exosome experimental
dataset.

### Level C — Mechanism-supportive
Independent experimental evidence increases the plausibility of an
exosome-related mechanism.

### Level D — Causal exosome mediation
Requires a design capable of identifying an exosome-mediated causal pathway.

The repository currently must not collapse Levels A-C into Level D.

The following are not, by themselves, evidence of exosome causality:

- plasma biomarker associations;
- rejuvenation-clock shifts;
- mouse exosome effects;
- cross-species concordance;
- inferred exosome fractions;
- pathway overlap.

The legacy `exosome_fraction_summary.csv` must therefore remain explicitly
described as a compatibility/exploratory artifact unless a defensible causal
decomposition becomes estimable.

---

## 5. Cross-Species Evidence Contract

### Purpose

The mouse exosome dataset is an orthogonal mechanism-support layer.

### Required checks

Before claiming concordance:

- define ortholog mapping;
- define comparable tissues or biological contexts;
- define compatible treatment contrasts;
- preserve effect direction;
- report available uncertainty;
- record missing or non-comparable features.

### Interpretation hierarchy

Allowed:

> "The macaque tissue signal is directionally aligned with an independent
> mouse exosome perturbation."

Stronger but still non-causal:

> "The cross-species concordance provides mechanism-supportive evidence
> consistent with an exosome-related process."

Not allowed from cross-species concordance alone:

> "Exosomes caused the macaque rejuvenation response."

Not allowed:

> "The same mechanism is established in humans."

Cross-species evidence changes biological plausibility, not the underlying
identifiability of the macaque causal effect.

---

## Claim Promotion Rule

A claim may only be promoted when:

1. the underlying data are valid;
2. required linkage is established;
3. the requested quantity is estimable;
4. the statistical result passes its predefined checks;
5. its evidence tier is correctly assigned;
6. human scientific review approves the wording.

Failure at an earlier stage cannot be repaired by stronger downstream
presentation.

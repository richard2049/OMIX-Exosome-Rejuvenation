# Data Linkage

## Why Linkage Matters

The repository combines tissue transcriptomics, plasma proteomics, methylation,
and supporting datasets. Samples from the same study group are not necessarily
measurements from the same animal. Group membership, file order, and similar
sample names therefore cannot substitute for demonstrated biological identity.

## Biological Units

| Analysis | Required unit |
|---|---|
| Transcriptomic age model | Animal, with repeated tissues kept within the same validation fold |
| Tissue treatment effect | Animal within tissue |
| Plasma biomarker ranking | Plasma sample; animal identity is retained when available |
| Plasma-to-tissue association | Linked animal |
| Mediation | One valid analytical row per linked animal unless repeated structure is modeled explicitly |
| Cross-species alignment | Tissue and contrast after explicit species-aware matching; not individual linkage |

## Required Evidence Artifacts

The pipeline writes four complementary linkage outputs:

| Output | Question answered |
|---|---|
| `plasma_to_animal_map.csv` | What mapping was proposed, by which rule, and with what confidence? |
| `linkage_qc_report.csv` | Are coverage, uniqueness, and collisions acceptable? |
| `linkage_audit.csv` | Which animal identifiers overlap across the relevant modalities? |
| `estimability_report.csv` | May a linkage-dependent analysis proceed under configured requirements? |

Actual coverage and gate status are profile-specific and should be read from
these outputs. A failed gate must produce a structured non-estimable result,
not a group-level surrogate.

## Permitted and Blocked Uses

- Unlinked plasma data may support descriptive or group-level exploratory
  summaries when labeled accordingly.
- Linked analyses must exclude unresolved, low-confidence, or colliding maps
  according to the configured rules.
- Tissue rows must not be replicated to manufacture multiple observations for
  one animal-level plasma measurement.
- Missing identifiers must not be filled from row order, treatment balance, or
  the nearest-looking sample name.

The [inference framework](inference_framework.md) is authoritative if this
reader-oriented guide and an implementation detail appear to conflict.

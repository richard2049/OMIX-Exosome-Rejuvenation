# Dataset Guide

## Purpose

This page summarizes how each public OMIX dataset contributes to the project.
Availability of a modality does not imply that cross-modal or causal analysis
is estimable; those decisions follow the
[inference framework](inference_framework.md).

## Dataset Roles

| Dataset | Organism and modality | Repository role | Current interpretation boundary |
|---|---|---|---|
| `OMIX007580` | Macaque bulk transcriptomics across tissues | Primary age-model and tissue treatment-effect analysis | Tissue estimates require uncertainty-aware interpretation and are not mechanism attribution |
| `OMIX007581` | Macaque plasma proteomics | Exploratory circulating biomarker analysis and linkage-aware plasma aging axis | Animal-level integration depends on explicit plasma-to-animal linkage |
| `OMIX007582` | Macaque Mammal40 DNA methylation | Intended orthogonal validation layer | Biological validation is blocked without a technical-to-biological sample map |
| `OMIX007583` | Macaque ovary expression subset | Targeted validation and naming support | Intervention-era `O_V` and stage control `O_C` must remain distinct |
| `OMIX007586` | Macaque hippocampus expression subset | Targeted validation and strongest article-to-OMIX naming bridge | Supports subset validation, not a new independent cohort |
| `OMIX009283` | Mouse exosome perturbation transcriptomics | Cross-species mechanism-supportive alignment | Cannot establish macaque or human exosome causality |
| `OMIX009284` | Macaque PBMC single-cell RNA data | Read-only structural audit and possible future PBMC analysis | PBMC-only scope does not resolve cross-tissue attribution |

## Data Profiles

- `demo` uses reduced, tracked examples from `data/PROCESSED` and is intended
  for reproducibility and smoke testing.
- `full` uses locally downloaded source data under `data/RAW/data`.
- `auto` selects the full profile when those inputs are present and otherwise
  falls back to the demo profile.

Raw source data are not distributed through Git and must remain immutable.
Generated tables and figures belong in `results/` and `figures/`.

## Reading Dataset Status

Dataset-specific sample counts and estimability can differ by profile. Consult
the generated linkage, audit, and estimability tables rather than treating this
guide as a frozen numerical results report. See the
[output reference](outputs.md) for the relevant files.

# AGENTS.md

## Scope
Repository-level operating notes for the clean `SRSC` repo. Keep edits consistent with the public/demo pipeline profile and with the current synchronization workflow from `SRSC - work`.

## Setup and execution
- Environment file: `environment.yml`
- Conda env name in file: `srsc`
- Main entrypoint: `src/run_pipeline.py`
- Public/demo run:
  ```bash
  conda env create -f environment.yml
  conda activate srsc
  python -m src.run_pipeline --profile demo
  ```
- Auto-detect profile:
  ```bash
  python -m src.run_pipeline --profile auto
  ```
- Conservative laptop run:
  ```bash
  python -m src.run_pipeline --profile demo --safe
  ```

## Validation and smoke tests
- Canonical guardrail test file: `tests/test_scientific_guardrails.py`
- Fast validation path:
  ```bash
  pytest -q tests/test_scientific_guardrails.py
  ```
- Public smoke test:
  ```bash
  python -m src.run_pipeline --profile demo --safe
  ```
- Confirm key outputs are regenerated in `results` and `figures`.

## Data, outputs, and layout
- Default clean-repo data root: `data/PROCESSED`
- Generated outputs go to:
  - `results`
  - `figures`
- This repo is the clean/demo distribution. Do not wire canonical behavior here to private full-data paths from `SRSC - work`.
- Raw full datasets are not redistributed here. Keep the clean repo runnable with example/sample-compatible files.

## Promotion and sync workflow
- Development with full data happens in `SRSC - work`.
- Promotion into this clean repo is curated and manifest-driven from the work repo through:
  - `src/repo_promotion.py`
  - `promotion_manifest.json`
- Do not manually rewrite OMIX paths in `src/run_pipeline.py` when preparing this repo for GitHub. Use profile-aware config instead.

## Canonical sources of truth
Before modifying analysis logic, metadata mappings, workflow rules, figure generation, or interpretation-related code, read:

- `Documents/study_design.md`
- `Documents/group_label_crosswalk.md`
- `Documents/OMIX007582_sample_map_audit.md`
- `README.md`
- `src/run_pipeline.py`
- `src/config.py`
- `src/omix_io.py`
- `src/linkage_audit.py`
- `tests/test_scientific_guardrails.py`

Use `Documents/study_design.md` as the working summary of study design.
Use `Documents/group_label_crosswalk.md` as the canonical article-to-repo naming legend.
If there is a conflict between docs, code, notebooks, metadata, or the paper, prefer the paper plus released OMIX/BioProject records and document the discrepancy explicitly instead of normalizing it silently.

## Study-specific guardrails
- Do not overwrite raw or downloaded source data.
- Treat sample/group mappings and metadata joins as high-risk operations.
- Do not silently change reported group definitions, endpoints, thresholds, or reference mappings.
- Clearly separate reproduced findings from exploratory extensions.
- Label any conclusion based on inferred or missing study details as tentative.

## High-risk areas
- Metadata joins and column auto-detection are conservative but fragile. Check joins around `sample_id`, `animal_id`, `group`, `tissue`, `age`, and `sex` before changing candidate-column logic.
- `OMIX007582` still does not have a defensible public biological sample map inside this repo. Do not promote positional or order-based methylation column renaming into canonical behavior.
- Plasma sample-to-animal mapping is explicit and confidence-scored. Treat `results/plasma_to_animal_map.csv`, `results/linkage_qc_report.csv`, `results/linkage_audit.csv`, and `results/estimability_report.csv` as required diagnostics, not optional side outputs.
- Keep causal claims gated by the estimability flow in `src/run_pipeline.py` and `src/linkage_audit.py`.
- Keep the clean repo focused on reproducible public/demo execution. Private full-data assumptions belong in `SRSC - work`, not here.

## Missing canonical study context
- There is still no dedicated machine-readable study-spec file defining all sample naming semantics, cohort design, or how every unused raw dataset should be integrated.
- If a study-specific assumption is not explicitly documented in `README.md`, `src/config.py`, or the current pipeline code, do not invent it.

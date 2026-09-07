# Setup and Reproducible Execution

## Environment

Create the main environment once:

```bash
conda env create -f environment.yml
conda activate srsc
```

If an existing environment must be updated after the specification changes:

```bash
conda env update -n srsc -f environment.yml --prune
```

## Data Profiles

The pipeline has one entry point and three data-selection profiles:

```bash
python -m src.run_pipeline --profile demo --safe
python -m src.run_pipeline --profile full
python -m src.run_pipeline --profile auto
```

- `demo` reads reduced public examples from `data/PROCESSED`.
- `full` reads locally downloaded inputs from `data/RAW/data`.
- `auto` prefers the full layout when present and otherwise uses the demo
  layout.
- `--safe` applies conservative laptop settings.
- `--data-root <PATH>` overrides the selected profile's data directory without
  changing source code.

Raw inputs are immutable. Generated tables and figures are written to
`results/` and `figures/`.

## Validation

Run the focused scientific guardrail suite:

```bash
pytest -q tests/test_scientific_guardrails.py
```

Run the pipeline and then regenerate the interpretation-facing report layer:

```bash
python -m src.run_pipeline --profile demo --safe
python -m src.report_figures
```

The report command reads existing `results/*.csv`; it does not recompute the
analysis or replace diagnostic plots.

When preparing a public release from reviewed full-profile results, refresh the
three tracked README figures explicitly:

```bash
python -m src.report_figures --portfolio-assets-dir docs/assets
```

Do not publish these snapshots from a reduced demo run unless they are clearly
labelled as demo outputs.

## Optional Audit Workflows

```bash
python -m src.omix007582_audit
python -m src.omix009283_metadata
python -m src.omix009284_audit
```

The optional Mammal40 IDAT rebuild uses R and Bioconductor. The main `srsc`
environment includes R and `BiocManager` as a bootstrap:

```bash
conda run -n srsc Rscript -e "BiocManager::install(c('sesame', 'sesameData', 'BiocParallel'), ask = FALSE, update = FALSE)"
conda run -n srsc Rscript src/scripts/process_OMIX007582_Mammal40.R --max-prefixes 1 --prep-candidates default --output-dir results/omix007582_rebuild_smoke
```

For isolation, create the optional R-only environment instead:

```bash
conda env create -f environment-omix007582-r.yml
conda run -n srsc-omix007582-r Rscript src/scripts/process_OMIX007582_Mammal40.R --help
```

The first rebuild may cache sesame reference resources through ExperimentHub.
Use `--max-prefixes` or `--prefixes` before attempting the complete archive.

## Full-Data Development and Public Promotion

The local full-data workspace and clean public repository use the same source
code and profile interface. Promotion is manifest-driven:

```bash
python -m src.repo_promotion
python -m src.repo_promotion --apply
python -m src.repo_promotion --include-assets
python -m src.repo_promotion --include-assets --apply
```

`promotion_manifest.json` controls code, documentation, tests, and environment
files, including explicit file-only removal of obsolete public artifacts.
`promotion_assets_manifest.json` separately controls reduced demo data. Run the
command without `--apply` first and inspect the proposed changes.

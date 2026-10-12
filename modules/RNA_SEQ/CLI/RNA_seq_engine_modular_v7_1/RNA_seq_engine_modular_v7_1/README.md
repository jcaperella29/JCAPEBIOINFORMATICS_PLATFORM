# RNA-seq Engine v7.1 — Modularized Bundle

This bundle is a mechanical modularization of the statistically cleared v7.1 candidate. The statistical methods were intentionally not redesigned during this split.

## Layout

- `rnaseq_cli.R` — CLI entry point and orchestration
- `R/validation.R` — CLI parsing/help
- `R/io.R` — matrix IO, sample alignment, QC, count preparation
- `R/design_de.R` — group/contrast resolution, design matrix, limma/voom DE
- `R/annotation.R` — biomaRt annotation and one-row-per-gene join invariant
- `R/enrichment.R` — Enrichr/g:Profiler enrichment, tested-gene background, network edges
- `R/classifier.R` — leakage-safe classifier preprocessing, train/test and CV paths
- `R/plots.R` — PCA, UMAP, volcano, heatmap, enrichment plots, HTML saving
- `R/power.R` — generic effect-size power reference
- `R/reporting.R` — logging, JSON, manifest, run summary/status, ZIP output
- `tests/statistical_release_regression_tests.R` — 7-part statistical regression suite adapted to the module tree

## First check

From the project directory:

```bash
Rscript tests/statistical_release_regression_tests.R .
```

Expected final line:

```text
PASS: targeted statistical release regression tests for items 1-5 plus positive-class reporting semantics.
```

## Airway CV example

If `sample_data(airways)` is beside `rnaseq_cli.R`:

```bash
Rscript rnaseq_cli.R \
  --counts "sample_data(airways)"/counts_data* \
  --phenotype "sample_data(airways)"/phenotype_data* \
  --phenotype-column dex \
  --reference-group untrt \
  --contrast trt-untrt \
  --species hsapiens \
  --classifier-validation cv \
  --seed 42 \
  --outdir airway_modular_cv_test
```

## Refactor rule

Treat the first modularization pass as behavior-preserving engineering work. Run the regression suite after structural changes, especially around `design_de.R`, `annotation.R`, `enrichment.R`, and `classifier.R`.

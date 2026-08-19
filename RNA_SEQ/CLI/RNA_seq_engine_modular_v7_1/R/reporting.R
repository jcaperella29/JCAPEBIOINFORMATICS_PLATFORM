# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

log_msg <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = ""), "\n")
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
}

write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

rel_path <- function(path) {
  tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) path)
}

file_info_record <- function(path, label = NULL, type = NULL) {
  if (!file.exists(path)) return(NULL)
  list(
    name = tools::file_path_sans_ext(basename(path)),
    path = sub(paste0("^", normalizePath(outdir, mustWork = FALSE), "/?"), "", normalizePath(path, mustWork = FALSE)),
    label = if (is.null(label)) basename(path) else label,
    type = if (is.null(type)) tools::file_ext(path) else type,
    bytes = unname(file.info(path)$size)
  )
}

statistical_notes <- function() {
  c(
    "Bulk RNA-seq DE uses edgeR::filterByExpr, TMM normalization, limma::voom precision weights, explicit contrasts, and limma empirical Bayes moderation.",
    "DE_results contains genes passing the configured FDR and absolute logFC thresholds; limma_all_results contains the full tested table.",
    "PCA/UMAP/heatmap use TMM-normalized logCPM values and are exploratory visualizations, not formal hypothesis tests.",
    "Enrichment uses the significant DE gene list with the DE-tested gene universe (genes surviving edgeR::filterByExpr) as the explicit background when the backend supports custom background input.",
    "Random Forest/logistic classifier outputs are exploratory; feature selection is performed inside each training split/fold and classifier logCPM preprocessing is sample-local so held-out samples cannot alter training-sample normalization.",
    "Generic effect-size power reference is an approximation for a generic t-test/ANOVA model and is not RNA-seq-specific, gene-wise, or FDR-aware power."
  )
}

make_manifest <- function(date_tag, enrichment_background = NULL) {
  table_files <- list.files(tables_dir, full.names = TRUE, recursive = FALSE)
  plot_files <- list.files(plots_dir, full.names = TRUE, recursive = FALSE)
  object_files <- list.files(objects_dir, full.names = TRUE, recursive = FALSE)
  list(
    module = "RNA_SEQ",
    version = "v7.1_stat_clearance_candidate",
    generated = as.character(Sys.time()),
    outdir = normalizePath(outdir, mustWork = FALSE),
    parameters = list(
      counts = counts_path,
      phenotype = phenotype_path,
      phenotype_column = phenotype_column,
      species = species_code,
      contrast = contrast_arg,
      reference_group = reference_group,
      covariates = covariates,
      fdr = fdr_cutoff,
      logfc_cutoff = logfc_cutoff,
      enrich_backend = enrich_backend,
      gprof_sources = gprof_sources,
      enrichr_db = enrichr_db,
      classifier_validation = classifier_validation,
      classifier_model = classifier_model,
      rf_top_n = rf_top_n,
      enrichment_background = enrichment_background
    ),
    tables = Filter(Negate(is.null), lapply(table_files, file_info_record, type = "table")),
    plots = Filter(Negate(is.null), lapply(plot_files, file_info_record, type = "html")),
    objects = Filter(Negate(is.null), lapply(object_files, file_info_record)),
    statistical_notes = statistical_notes()
  )
}

write_run_status <- function(status, started, finished = NULL, error = NULL, zip_path = NULL) {
  obj <- list(
    module = "RNA_SEQ",
    version = "v7.1_stat_clearance_candidate",
    started = as.character(started),
    finished = if (is.null(finished)) NULL else as.character(finished),
    status = status,
    outdir = normalizePath(outdir, mustWork = FALSE),
    error = error,
    zip_path = zip_path,
    parameters = list(
      counts = counts_path,
      phenotype = phenotype_path,
      phenotype_column = phenotype_column,
      species = species_code,
      fdr = fdr_cutoff,
      logfc_cutoff = logfc_cutoff,
      classifier_validation = classifier_validation,
      enrich_backend = enrich_backend
    )
  )
  write_json_file(obj, file.path(outdir, "run_status.json"))
}

write_run_summary <- function(de_all = NULL, final_results = NULL, enrichment_counts = list(), rf_metrics = NULL, matched_samples = NA_integer_, enrichment_background = NULL) {
  summary <- list(
    module = "RNA_SEQ",
    version = "v7.1_stat_clearance_candidate",
    generated = as.character(Sys.time()),
    matched_samples = matched_samples,
    tested_genes = if (is.null(de_all)) NA_integer_ else nrow(de_all),
    significant_genes = if (is.null(final_results)) NA_integer_ else nrow(final_results),
    enrichment_rows = enrichment_counts,
    enrichment_background = enrichment_background,
    classifier_metrics = rf_metrics,
    statistical_notes = statistical_notes()
  )
  write_json_file(summary, file.path(outdir, "run_summary.json"))
}

zip_results <- function() {
  zip_path <- paste0(outdir, ".zip")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  parent <- dirname(normalizePath(outdir, mustWork = FALSE))
  base <- basename(normalizePath(outdir, mustWork = FALSE))
  setwd(parent)
  if (file.exists(zip_path)) unlink(zip_path)
  utils::zip(zipfile = zip_path, files = base, flags = "-r9Xq")
  zip_path
}


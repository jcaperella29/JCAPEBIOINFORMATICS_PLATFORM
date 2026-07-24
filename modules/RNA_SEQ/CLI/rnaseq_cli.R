#!/usr/bin/env Rscript
# RNA_SEQ_CLI_V6_API_MANIFEST_CV_ZIP_2026_07_02
# GPROFILER_GENE_EDGE_FIX_V5_2026_07_02
# NUMERIC_CONVERSION_FIX_V4_ENRICHMENT_ID_FALLBACK_2026_07_02

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(biomaRt)
  library(enrichR)
  library(gprofiler2)
  library(randomForest)
  library(caret)
  library(pROC)
  library(pwr)
  library(doParallel)
  library(foreach)
  library(umap)
  library(ggplot2)
  library(dplyr)
  library(plotly)
  library(pheatmap)
  library(jsonlite)
})

parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!grepl("^--", key)) stop("Unexpected positional argument: ", key)
    key <- sub("^--", "", key)
    if (key %in% c("help", "skip-enrichment", "skip-annotation", "skip-classifier", "skip-plots", "drop-library-outliers", "zip")) {
      args[[key]] <- TRUE
      i <- i + 1
    } else {
      if (i == length(argv)) stop("Missing value for --", key)
      args[[key]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  args
}

show_help <- function() {
  cat("
JCAP RNA-SEQ Analyzer CLI - production statistics edition

Required:
  --counts PATH
  --phenotype PATH
  --phenotype-column NAME

Strongly recommended for production:
  --reference-group NAME                            Reference/control group for 2-group DE
  --contrast GROUP_B-GROUP_A                        Explicit comparison. Example: treated-control

Optional:
  --covariates COL1,COL2                            Optional design covariates, e.g. batch,sex,RIN
  --species hsapiens|mmusculus|drerio|dmelanogaster  [default: hsapiens]
  --outdir PATH                                      [default: rnaseq_cli_output]
  --enrich-db NAME                                   Backward-compatible enrichment DB/source
  --enrich-backend gprof|enrichr|auto                [default: auto]
  --gprof-sources GO:BP,GO:MF,GO:CC,KEGG,REAC        [default: GO:BP]
  --enrichr-db NAME                                  [default: GO_Biological_Process_2023]
  --fdr NUMBER                                       [default: 0.05]
  --logfc-cutoff NUMBER                              [default: 1]
  --effect-size NUMBER                               [default: 0.8]
  --power-test-type ttest|anova                      [default: ttest]
  --curve-n-min INTEGER                              [default: 2]
  --curve-n-max INTEGER                              [default: 30]
  --rf-top-n INTEGER                                 [default: 50]
  --classifier-validation cv|train_test|none          [default: cv]
  --classifier-model auto|rf|logistic                 [default: auto]
  --seed INTEGER                                     [default: 42]
  --drop-library-outliers                            Opt-in sample dropping based on library-size IQR QC
  --skip-enrichment
  --skip-annotation
  --skip-classifier
  --skip-plots
  --zip                                             Zip output directory at end
  --help

Example:
  Rscript rnaseq_cli_production_fixed.R \
    --counts counts.csv \
    --phenotype phenotype.csv \
    --phenotype-column condition \
    --reference-group control \
    --species hsapiens \
    --outdir results

Outputs:
  run_status.json / run_summary.json / manifest.json
  tables/DE_results_<date>.csv                       Significant DE genes only
  tables/limma_all_results_<date>.csv                Full limma/voom contrast table
  tables/normalized_logCPM_<date>.csv                TMM-normalized logCPM matrix
  tables/sample_qc_<date>.csv                        Library-size QC table
  tables/classifier_feature_stability_<date>.csv     CV feature stability
  tables/rf_feature_importance_<date>.csv            Exploratory RF feature importance
  tables/rf_predictions_<date>.csv
  tables/rf_metrics_<date>.csv
  tables/exploratory_power_summary_<date>.csv
  tables/enrichment_all.csv / upregulated / downregulated
  tables/enrichment_*_network_edges.csv
  plots/pca_plot.html / umap_plot.html / volcano_plot.html / heatmap_plot.html
  plots/rf_roc_plot.html
  objects/dge.rds / voom_object.rds / limma_fit.rds / contrast_matrix.rds
  sessionInfo.txt
  run_log.txt
")
}

args <- parse_args()
if (isTRUE(args$help)) {
  show_help()
  quit(status = 0)
}

required_args <- c("counts", "phenotype", "phenotype-column")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args) > 0) {
  show_help()
  stop("Missing required argument(s): ", paste0("--", missing_args, collapse = ", "))
}

counts_path <- args$counts
phenotype_path <- args$phenotype
phenotype_column <- args[["phenotype-column"]]
species_code <- if (!is.null(args$species)) args$species else "hsapiens"
outdir <- if (!is.null(args$outdir)) args$outdir else "rnaseq_cli_output"
effect_size <- if (!is.null(args[["effect-size"]])) as.numeric(args[["effect-size"]]) else 0.8
power_test_type <- if (!is.null(args[["power-test-type"]])) args[["power-test-type"]] else "ttest"
curve_n_min <- if (!is.null(args[["curve-n-min"]])) as.integer(args[["curve-n-min"]]) else 2L
curve_n_max <- if (!is.null(args[["curve-n-max"]])) as.integer(args[["curve-n-max"]]) else 30L
seed <- if (!is.null(args$seed)) as.integer(args$seed) else 42L
fdr_cutoff <- if (!is.null(args$fdr)) as.numeric(args$fdr) else 0.05
logfc_cutoff <- if (!is.null(args[["logfc-cutoff"]])) as.numeric(args[["logfc-cutoff"]]) else 1
rf_top_n <- if (!is.null(args[["rf-top-n"]])) as.integer(args[["rf-top-n"]]) else 50L
enrich_backend <- if (!is.null(args[["enrich-backend"]])) tolower(args[["enrich-backend"]]) else "auto"
gprof_sources <- if (!is.null(args[["gprof-sources"]])) trimws(strsplit(args[["gprof-sources"]], ",")[[1]]) else c("GO:BP")
enrichr_db <- if (!is.null(args[["enrichr-db"]])) args[["enrichr-db"]] else "GO_Biological_Process_2023"
classifier_validation <- if (!is.null(args[["classifier-validation"]])) tolower(args[["classifier-validation"]]) else "cv"
classifier_model <- if (!is.null(args[["classifier-model"]])) tolower(args[["classifier-model"]]) else "auto"
zip_output <- isTRUE(args[["zip"]])
reference_group <- args[["reference-group"]]
contrast_arg <- args[["contrast"]]
covariates <- if (!is.null(args$covariates) && nzchar(args$covariates)) trimws(strsplit(args$covariates, ",")[[1]]) else character(0)

skip_enrichment <- isTRUE(args[["skip-enrichment"]])
skip_annotation <- isTRUE(args[["skip-annotation"]])
skip_classifier <- isTRUE(args[["skip-classifier"]])
skip_plots <- isTRUE(args[["skip-plots"]])
drop_library_outliers <- isTRUE(args[["drop-library-outliers"]])

valid_species <- c("hsapiens", "mmusculus", "drerio", "dmelanogaster")
if (!species_code %in% valid_species) stop("--species must be one of: ", paste(valid_species, collapse = ", "))
if (!power_test_type %in% c("ttest", "anova")) stop("--power-test-type must be ttest or anova")
if (curve_n_min < 2 || curve_n_max < curve_n_min) stop("Invalid power curve range.")
if (is.na(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff >= 1) stop("--fdr must be between 0 and 1.")
if (is.na(logfc_cutoff) || logfc_cutoff < 0) stop("--logfc-cutoff must be >= 0.")
if (!enrich_backend %in% c("auto", "gprof", "enrichr")) stop("--enrich-backend must be auto, gprof, or enrichr")
if (!classifier_validation %in% c("cv", "train_test", "none")) stop("--classifier-validation must be cv, train_test, or none")
if (!classifier_model %in% c("auto", "rf", "logistic")) stop("--classifier-model must be auto, rf, or logistic")

if (drop_library_outliers) {
  warning("--drop-library-outliers is enabled. This can remove biologically valid samples. Use only after QC review.")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(outdir, "tables"); dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(outdir, "plots"); dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
objects_dir <- file.path(outdir, "objects"); dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(outdir, "run_log.txt")

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
    "Enrichment uses the significant DE gene list and the selected public enrichment backend/background; if biomaRt annotation fails, Ensembl IDs are passed to g:Profiler.",
    "Random Forest/logistic classifier outputs are exploratory; CV feature selection is performed inside each training fold when --classifier-validation cv is used.",
    "Generic power analysis is an effect-size approximation and is not RNA-seq-specific FDR-aware/gene-wise power."
  )
}

make_manifest <- function(date_tag) {
  table_files <- list.files(tables_dir, full.names = TRUE, recursive = FALSE)
  plot_files <- list.files(plots_dir, full.names = TRUE, recursive = FALSE)
  object_files <- list.files(objects_dir, full.names = TRUE, recursive = FALSE)
  list(
    module = "RNA_SEQ",
    version = "v6_api_manifest_cv_zip",
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
      rf_top_n = rf_top_n
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
    version = "v6_api_manifest_cv_zip",
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

write_run_summary <- function(de_all = NULL, final_results = NULL, enrichment_counts = list(), rf_metrics = NULL, matched_samples = NA_integer_) {
  summary <- list(
    module = "RNA_SEQ",
    version = "v6_api_manifest_cv_zip",
    generated = as.character(Sys.time()),
    matched_samples = matched_samples,
    tested_genes = if (is.null(de_all)) NA_integer_ else nrow(de_all),
    significant_genes = if (is.null(final_results)) NA_integer_ else nrow(final_results),
    enrichment_rows = enrichment_counts,
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

save_plotly_html <- function(plot_obj, path) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Package htmlwidgets is required.")
  htmlwidgets::saveWidget(plot_obj, file = path, selfcontained = TRUE)
}

read_matrix_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    read.csv(path, row.names = 1, check.names = FALSE)
  } else {
    read.table(path, header = TRUE, row.names = 1, check.names = FALSE, sep = "", quote = "\"", comment.char = "")
  }
}

strip_ens_version <- function(x) sub("\\..*", "", x)

bm_dataset_for <- function(sp) {
  switch(
    sp,
    "hsapiens"      = "hsapiens_gene_ensembl",
    "mmusculus"     = "mmusculus_gene_ensembl",
    "drerio"        = "drerio_gene_ensembl",
    "dmelanogaster" = "dmelanogaster_gene_ensembl",
    "hsapiens_gene_ensembl"
  )
}

gp_org_for <- function(sp) {
  switch(sp,
         "hsapiens" = "hsapiens",
         "mmusculus" = "mmusculus",
         "drerio" = "drerio",
         "dmelanogaster" = "dmelanogaster",
         "hsapiens")
}

get_species_enrich_choices <- function(species_code) {
  if (species_code == "hsapiens") {
    c("GO_Biological_Process_2021", "KEGG_2021_Human", "WikiPathway_2021_Human", "Reactome_2022")
  } else {
    c("GO: Biological Process" = "GO:BP", "KEGG" = "KEGG", "Reactome" = "REAC", "WikiPathways" = "WP")
  }
}

default_enrich_db <- function(species_code) unname(get_species_enrich_choices(species_code)[[1]])

annotate_genes <- function(ensembl_ids, species_code) {
  ens_ids <- unique(strip_ens_version(ensembl_ids))
  dataset <- bm_dataset_for(species_code)
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  ann <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = ens_ids,
    mart = mart
  )
  names(ann) <- c("ensembl_gene_id", "symbol", "description")
  ann
}

align_counts_pheno <- function(counts, pheno) {
  common <- intersect(colnames(counts), rownames(pheno))
  if (length(common) < 2) stop("Fewer than two matching sample IDs between counts columns and phenotype row names.")
  list(counts = counts[, common, drop = FALSE], pheno = pheno[common, , drop = FALSE])
}

make_sample_qc <- function(counts, pheno, phenotype_column) {
  library_size <- colSums(counts, na.rm = TRUE)
  detected_genes <- colSums(counts > 0, na.rm = TRUE)
  q1 <- stats::quantile(library_size, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(library_size, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  data.frame(
    sample_id = colnames(counts),
    phenotype = as.character(pheno[colnames(counts), phenotype_column]),
    library_size = as.numeric(library_size),
    detected_genes = as.numeric(detected_genes),
    library_size_iqr_lower = as.numeric(lower),
    library_size_iqr_upper = as.numeric(upper),
    suggested_library_outlier = as.logical(library_size < lower | library_size > upper),
    stringsAsFactors = FALSE
  )
}

prepare_counts <- function(counts, pheno, phenotype_column, date_tag) {
  # read.csv/read.table returns a data.frame. Convert each sample column explicitly.
  # Do not use storage.mode() on a data.frame, because data.frames are lists in R.
  counts_df <- as.data.frame(counts, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- rownames(counts_df)
  cn <- colnames(counts_df)

  numeric_cols <- lapply(seq_along(counts_df), function(j) {
    x <- counts_df[[j]]
    if (is.factor(x)) x <- as.character(x)
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    y <- suppressWarnings(as.numeric(x))
    bad <- is.na(y) & !is.na(x) & trimws(as.character(x)) != ""
    if (any(bad)) {
      stop("Counts matrix contains non-numeric values in sample column: ", cn[[j]])
    }
    y
  })

  counts <- do.call(cbind, numeric_cols)
  counts <- as.matrix(counts)
  rownames(counts) <- rn
  colnames(counts) <- cn
  storage.mode(counts) <- "double"
  if (anyNA(counts)) stop("Counts matrix contains NA values.")
  if (any(counts < 0)) stop("Counts matrix contains negative values.")
  if (any(abs(counts - round(counts)) > .Machine$double.eps^0.5)) {
    log_msg("WARNING: Counts matrix contains non-integer values. edgeR/voom expects raw integer-like counts.")
  }

  qc <- make_sample_qc(counts, pheno, phenotype_column)
  qc_path <- file.path(tables_dir, paste0("sample_qc_", date_tag, ".csv"))
  write.csv(qc, qc_path, row.names = FALSE)
  log_msg("Wrote sample QC: ", qc_path)

  if (drop_library_outliers && ncol(counts) >= 10) {
    keep_samples <- qc$sample_id[!qc$suggested_library_outlier]
    dropped <- setdiff(colnames(counts), keep_samples)
    if (length(keep_samples) < 2) stop("Library-size outlier removal would leave fewer than two samples.")
    if (length(dropped) > 0) log_msg("Dropped library-size QC outliers by user request: ", paste(dropped, collapse = ", "))
    counts <- counts[, keep_samples, drop = FALSE]
    pheno <- pheno[keep_samples, , drop = FALSE]
  } else {
    n_out <- sum(qc$suggested_library_outlier, na.rm = TRUE)
    if (n_out > 0) log_msg("QC flagged ", n_out, " library-size outlier(s); no samples dropped unless --drop-library-outliers is set.")
  }

  keep_total <- rowSums(counts) > 10
  counts <- counts[keep_total, , drop = FALSE]
  if (nrow(counts) < 2) stop("Fewer than two genes remain after simple total-count filtering.")
  list(counts = counts, pheno = pheno, sample_qc = qc)
}

resolve_groups <- function(pheno, phenotype_column, reference_group = NULL, contrast_arg = NULL) {
  group <- factor(pheno[[phenotype_column]])
  names(group) <- rownames(pheno)
  if (length(levels(group)) < 2) stop("Phenotype column must contain at least two groups.")

  if (!is.null(reference_group)) {
    if (!reference_group %in% levels(group)) stop("Reference group not found in phenotype column: ", reference_group)
    group <- stats::relevel(group, ref = reference_group)
  }

  if (!is.null(contrast_arg)) {
    parts <- strsplit(contrast_arg, "-", fixed = TRUE)[[1]]
    if (length(parts) != 2 || any(!nzchar(parts))) stop("--contrast must look like GROUP_B-GROUP_A, e.g. treated-control")
    numerator <- parts[[1]]
    denominator <- parts[[2]]
    if (!all(c(numerator, denominator) %in% levels(group))) {
      stop("Contrast group(s) not found. Available groups: ", paste(levels(group), collapse = ", "))
    }
  } else if (length(levels(group)) == 2) {
    denominator <- levels(group)[[1]]
    numerator <- levels(group)[[2]]
    log_msg("No --contrast supplied; using two-group contrast: ", numerator, "-", denominator)
  } else {
    stop("More than two phenotype groups detected. Supply --contrast GROUP_B-GROUP_A.")
  }

  list(group = group, numerator = numerator, denominator = denominator)
}

make_design <- function(pheno, group, covariates) {
  missing_covars <- setdiff(covariates, colnames(pheno))
  if (length(missing_covars) > 0) stop("Covariate column(s) not found: ", paste(missing_covars, collapse = ", "))

  design_df <- pheno
  design_df$.group <- group

  for (cv in covariates) {
    if (is.character(design_df[[cv]])) design_df[[cv]] <- factor(design_df[[cv]])
  }

  rhs <- c("0 + .group", covariates)
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(form, data = design_df)

  group_cols <- grep("^\\.group", colnames(design))
  old_group_cols <- colnames(design)[group_cols]
  for (i in seq_along(old_group_cols)) {
    lvl <- sub("^\\.group", "", old_group_cols[[i]])
    colnames(design)[group_cols[[i]]] <- make.names(lvl)
  }

  if (qr(design)$rank < ncol(design)) {
    stop("Design matrix is not full rank. Check phenotype groups and covariates for confounding or empty levels.")
  }

  design
}

annotate_de_table <- function(de_table, species_code, skip_annotation = FALSE) {
  de_table$Ensembl_IDs <- strip_ens_version(de_table$Ensembl_IDs)
  de_table$.row_order <- seq_len(nrow(de_table))
  if (skip_annotation) {
    de_table$symbol <- de_table$Ensembl_IDs
    de_table$description <- NA_character_
    de_table <- de_table[order(de_table$.row_order), , drop = FALSE]
    de_table$.row_order <- NULL
    return(de_table)
  }

  log_msg("Annotating genes with biomaRt...")
  ann <- tryCatch(
    annotate_genes(de_table$Ensembl_IDs, species_code),
    error = function(e) {
      log_msg("Gene annotation failed; continuing with Ensembl IDs as symbols: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(ann) || nrow(ann) == 0) {
    de_table$symbol <- de_table$Ensembl_IDs
    de_table$description <- NA_character_
    de_table <- de_table[order(de_table$.row_order), , drop = FALSE]
    de_table$.row_order <- NULL
    return(de_table)
  }

  out <- merge(de_table, ann, by.x = "Ensembl_IDs", by.y = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$.row_order), , drop = FALSE]
  out$.row_order <- NULL
  out$symbol[is.na(out$symbol) | out$symbol == ""] <- out$Ensembl_IDs[is.na(out$symbol) | out$symbol == ""]
  out
}

run_differential_expression <- function(counts, pheno, phenotype_column, species_code, covariates = character(0), reference_group = NULL, contrast_arg = NULL, skip_annotation = FALSE) {
  if (!phenotype_column %in% colnames(pheno)) stop("Phenotype column not found: ", phenotype_column)

  resolved <- resolve_groups(pheno, phenotype_column, reference_group, contrast_arg)
  group <- resolved$group

  dge <- DGEList(counts = counts, group = group)
  design <- make_design(pheno, group, covariates)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 2) stop("Fewer than two genes remain after edgeR::filterByExpr.")
  dge <- calcNormFactors(dge)

  v <- voom(dge, design, plot = FALSE)

  numerator <- make.names(resolved$numerator)
  denominator <- make.names(resolved$denominator)
  if (!all(c(numerator, denominator) %in% colnames(design))) {
    stop("Resolved contrast columns not found in design. Design columns: ", paste(colnames(design), collapse = ", "))
  }
  contrast_string <- paste0(numerator, "-", denominator)
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)

  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  res$Ensembl_IDs <- rownames(res)
  res$contrast <- contrast_string
  res$comparison <- paste0(resolved$numerator, " vs ", resolved$denominator)

  res_annot <- annotate_de_table(res, species_code, skip_annotation)
  sig <- res_annot[!is.na(res_annot$adj.P.Val) & res_annot$adj.P.Val < fdr_cutoff & abs(res_annot$logFC) >= logfc_cutoff, , drop = FALSE]
  sig <- sig[order(sig$adj.P.Val, -abs(sig$logFC)), , drop = FALSE]

  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

  list(
    de_all = res_annot,
    de_sig = sig,
    logcpm = logcpm,
    counts_filtered = dge$counts,
    pheno_final = pheno[colnames(dge$counts), , drop = FALSE],
    group = group[colnames(dge$counts)],
    contrast = contrast_string,
    comparison = paste0(resolved$numerator, " vs ", resolved$denominator),
    dge = dge,
    voom_object = v,
    limma_fit = fit,
    contrast_matrix = contrast_matrix
  )
}

make_pca_plot <- function(logcpm, pheno, phenotype_column) {
  df <- data.frame(t(logcpm))
  df$Phenotype <- as.factor(pheno[rownames(df), phenotype_column])
  df_vars <- df[, -ncol(df), drop = FALSE]
  zero_var_cols <- apply(df_vars, 2, function(x) var(x) == 0)
  if (all(zero_var_cols)) stop("All features have zero variance for PCA.")
  df_clean <- df_vars[, !zero_var_cols, drop = FALSE]
  pca <- prcomp(df_clean, scale. = TRUE)
  pca_df <- data.frame(pca$x, Sample = rownames(df), Phenotype = df$Phenotype)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype, text = Sample)) + geom_point() + theme_minimal()
  ggplotly(p, tooltip = "text")
}

make_umap_plot <- function(logcpm, pheno, phenotype_column) {
  df <- data.frame(t(logcpm))
  df$Phenotype <- as.factor(pheno[rownames(df), phenotype_column])
  n_samples <- nrow(df)
  if (n_samples < 3) stop("UMAP requires at least 3 samples.")
  n_neighbors <- min(15, max(2, n_samples - 1))
  umap_config <- umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_res <- umap(df[, -ncol(df), drop = FALSE], config = umap_config)
  umap_df <- data.frame(umap_res$layout, Sample = rownames(df), Phenotype = df$Phenotype)
  p <- ggplot(umap_df, aes(x = X1, y = X2, color = Phenotype, text = Sample)) + geom_point() + theme_minimal()
  ggplotly(p, tooltip = "text")
}

make_volcano_plot <- function(de_all) {
  df <- de_all
  df$adj.P.Val[df$adj.P.Val == 0] <- .Machine$double.xmin
  df$log10p <- -log10(df$adj.P.Val)
  df$significant <- ifelse(df$adj.P.Val < fdr_cutoff & abs(df$logFC) >= logfc_cutoff, "Significant", "Not Significant")
  label_col <- if ("symbol" %in% names(df)) df$symbol else df$Ensembl_IDs
  df$label <- label_col
  p <- ggplot(df, aes(x = logFC, y = log10p, color = significant, text = label)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("gray70", "firebrick")) +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value") +
    theme_minimal()
  ggplotly(p, tooltip = "text")
}

make_heatmap_plot <- function(logcpm, de_sig, de_all) {
  if (!is.null(de_sig) && nrow(de_sig) >= 2) {
    final_df <- de_sig
  } else {
    log_msg("Fewer than two significant genes for heatmap; using top genes from full DE table.")
    final_df <- head(de_all[order(de_all$adj.P.Val), , drop = FALSE], 50)
  }
  logcpm2 <- logcpm
  rownames(logcpm2) <- strip_ens_version(rownames(logcpm2))
  final_df$Ensembl_IDs <- strip_ens_version(final_df$Ensembl_IDs)
  matched_genes <- intersect(rownames(logcpm2), final_df$Ensembl_IDs)
  if (length(matched_genes) < 2) stop("Not enough DE genes matched in logCPM matrix for heatmap.")
  mat <- logcpm2[matched_genes, , drop = FALSE]
  mat <- t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
  mat[is.na(mat)] <- 0
  row_order <- hclust(dist(mat))$order
  col_order <- hclust(dist(t(mat)))$order
  mat <- mat[row_order, col_order, drop = FALSE]
  plot_ly(
    z = mat, x = colnames(mat), y = rownames(mat), type = "heatmap",
    colorscale = "Viridis",
    colorbar = list(title = "Z-score"),
    hovertemplate = paste("Gene: %{y}<br>", "Sample: %{x}<br>", "Value: %{z:.2f}<extra></extra>")
  ) %>% layout(title = "Heatmap - DE Genes / Top Genes (logCPM Z-score)", xaxis = list(title = "Samples"), yaxis = list(title = "Genes"), margin = list(l = 100, b = 100))
}

is_ensembl_like <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^ENS[A-Z]*G[0-9]+", x)) >= 0.5
}

gprofiler_source_for_db <- function(db) {
  db_txt <- toupper(as.character(db))
  if (db_txt %in% c("GO:BP", "GO_BP", "GOBP") || grepl("GO", db_txt) || grepl("BIOLOGICAL", db_txt)) return("GO:BP")
  if (db_txt %in% c("KEGG") || grepl("KEGG", db_txt)) return("KEGG")
  if (db_txt %in% c("REAC", "REACTOME") || grepl("REACTOME", db_txt)) return("REAC")
  if (db_txt %in% c("WP", "WIKIPATHWAYS") || grepl("WIKI", db_txt)) return("WP")
  "GO:BP"
}

format_gprofiler_enrichment <- function(g, source_label = "g:Profiler") {
  if (is.null(g) || is.null(g$result) || nrow(g$result) == 0) return(NULL)
  df <- g$result
  comb <- -log10(df$p_value) * (df$intersection_size / df$term_size)
  genes_col <- vapply(seq_len(nrow(df)), function(i) {
    x <- df$intersection[[i]]
    if (is.null(x) || length(x) == 0 || all(is.na(x))) return("")
    paste(unique(as.character(x)), collapse = ";")
  }, character(1))
  gene_count_col <- if ("intersection_size" %in% names(df)) df$intersection_size else vapply(strsplit(genes_col, ";", fixed = TRUE), function(x) sum(x != ""), integer(1))
  out <- data.frame(
    Term = df$term_name,
    Genes = genes_col,
    Gene.Count = gene_count_col,
    Adjusted.P.value = df$p_value,
    Combined.Score = comb,
    Source = source_label,
    stringsAsFactors = FALSE
  )
  out[order(out$Adjusted.P.value), ]
}

run_gprofiler_enrichment <- function(gene_list, db, species_code) {
  source <- gprofiler_source_for_db(db)
  org <- gp_org_for(species_code)
  log_msg("Running g:Profiler enrichment fallback/source ", source, " with ", length(gene_list), " identifiers...")
  g <- tryCatch(
    gprofiler2::gost(query = gene_list, organism = org, sources = source, correction_method = "g_SCS", evcodes = TRUE),
    error = function(e) {
      log_msg("g:Profiler enrichment failed for ", source, ": ", conditionMessage(e))
      NULL
    }
  )
  format_gprofiler_enrichment(g, source_label = paste0("g:Profiler ", source))
}

perform_enrichment <- function(gene_list, db, species_code) {
  gene_list <- unique(na.omit(as.character(gene_list)))
  gene_list <- gene_list[gene_list != ""]
  if (length(gene_list) < 2) return(NULL)

  # Production behavior:
  # - Enrichr generally prefers human gene symbols.
  # - If BioMart annotation fails and we only have Ensembl IDs, use g:Profiler,
  #   which can resolve Ensembl IDs directly. This keeps enrichment functional
  #   without making annotation a hard dependency.
  if (species_code == "hsapiens" && !is_ensembl_like(gene_list)) {
    enrichr_res <- tryCatch(enrichR::enrichr(gene_list, databases = db)[[1]], error = function(e) {
      log_msg("Human Enrichr failed for ", db, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
      if (!"Source" %in% names(enrichr_res)) enrichr_res$Source <- paste0("Enrichr ", db)
      return(enrichr_res)
    }
    log_msg("Enrichr returned no usable results; trying g:Profiler fallback.")
    return(run_gprofiler_enrichment(gene_list, db, species_code))
  }

  run_gprofiler_enrichment(gene_list, db, species_code)
}


perform_enrichment_v6 <- function(gene_list, db, species_code, backend = "auto", gprof_sources = c("GO:BP"), enrichr_db = "GO_Biological_Process_2023") {
  gene_list <- unique(na.omit(as.character(gene_list)))
  gene_list <- gene_list[gene_list != ""]
  if (length(gene_list) < 2) return(NULL)
  backend <- tolower(backend)

  if (backend %in% c("auto", "enrichr") && species_code == "hsapiens" && !is_ensembl_like(gene_list)) {
    er_db <- if (!is.null(enrichr_db) && nzchar(enrichr_db)) enrichr_db else db
    enrichr_res <- tryCatch(enrichR::enrichr(gene_list, databases = er_db)[[1]], error = function(e) {
      log_msg("Enrichr failed for ", er_db, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
      if (!"Source" %in% names(enrichr_res)) enrichr_res$Source <- paste0("Enrichr ", er_db)
      return(enrichr_res)
    }
    if (backend == "enrichr") return(NULL)
    log_msg("Enrichr returned no usable results; trying g:Profiler fallback.")
  }

  org <- gp_org_for(species_code)
  sources <- unique(trimws(gprof_sources))
  sources <- sources[nzchar(sources)]
  if (length(sources) < 1) sources <- gprof_source_for_db(db)
  log_msg("Running g:Profiler enrichment source(s) ", paste(sources, collapse = ","), " with ", length(gene_list), " identifiers...")
  g <- tryCatch(
    gprofiler2::gost(query = gene_list, organism = org, sources = sources, correction_method = "g_SCS", evcodes = TRUE),
    error = function(e) {
      log_msg("g:Profiler enrichment failed: ", conditionMessage(e))
      NULL
    }
  )
  format_gprofiler_enrichment(g, source_label = paste0("g:Profiler ", paste(sources, collapse = ",")))
}

plot_enrichment_bar <- function(df, title) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  top <- head(df[order(df$Adjusted.P.value), ], 10)
  top$Term <- factor(top$Term, levels = rev(top$Term))
  plotly::plot_ly(
    data = top, x = ~-log10(Adjusted.P.value), y = ~Term,
    type = "bar", orientation = "h", hoverinfo = "text",
    text = ~paste0("P.adj: ", signif(Adjusted.P.value, 3), "<br>Score: ", round(Combined.Score, 2))
  ) %>% layout(title = list(text = title), xaxis = list(title = "-log10 Adjusted P-value"), yaxis = list(title = ""), margin = list(l = 200))
}

normalize_col_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

find_first_col <- function(df, candidates) {
  norm_names <- normalize_col_name(names(df))
  norm_candidates <- normalize_col_name(candidates)
  hit <- match(norm_candidates, norm_names, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) == 0) return(NULL)
  names(df)[hit[[1]]]
}

make_network_edges <- function(enrich_df) {
  empty <- data.frame(gene = character(), term = character(), adjusted_pvalue = numeric(), stringsAsFactors = FALSE)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(empty)
  if ("Message" %in% names(enrich_df)) return(empty)

  term_col <- find_first_col(enrich_df, c("term", "Term", "term_name", "name", "description"))
  genes_col <- find_first_col(enrich_df, c("genes", "Genes", "intersection", "overlapping genes", "overlap_genes", "core_enrichment"))
  padj_col <- find_first_col(enrich_df, c("adjusted_pvalue", "Adjusted.P.value", "Adjusted P-value", "adjusted p value", "padj", "p.adjust", "p_adjust", "qvalue", "fdr", "p_value", "pvalue"))

  if (is.null(term_col) || is.null(genes_col) || is.null(padj_col)) {
    stop("Cannot create network edges. Need term, genes/intersection, and adjusted p-value columns. Found: ", paste(names(enrich_df), collapse = ", "))
  }

  rows <- list()
  k <- 1L
  for (i in seq_len(nrow(enrich_df))) {
    term <- as.character(enrich_df[[term_col]][[i]])
    padj <- suppressWarnings(as.numeric(enrich_df[[padj_col]][[i]]))
    raw_genes <- enrich_df[[genes_col]][[i]]
    if (is.null(raw_genes) || length(raw_genes) == 0 || is.na(term) || term == "") next
    genes_text <- paste(as.character(raw_genes), collapse = ";")
    genes <- unlist(strsplit(genes_text, "[;,/|]+"))
    genes <- trimws(genes)
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) next
    for (gene in unique(genes)) {
      rows[[k]] <- data.frame(gene = gene, term = term, adjusted_pvalue = padj, stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }

  if (length(rows) == 0) return(empty)
  unique(do.call(rbind, rows))
}

run_power_analysis <- function(pheno, phenotype_column, effect_size) {
  groups <- as.factor(pheno[[phenotype_column]])
  group_sizes <- table(groups)
  k <- length(group_sizes)
  sig <- 0.05
  note <- "Generic effect-size approximation only; not RNA-seq-specific gene-wise/FDR power."
  if (k == 2) {
    n1 <- as.numeric(group_sizes[1])
    n2 <- as.numeric(group_sizes[2])
    power_res <- pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = effect_size, sig.level = sig)
    data.frame(Test = "t-test", `Effect Size (d)` = effect_size, Groups = paste(names(group_sizes), collapse = ", "), n1 = n1, n2 = n2, `Significance Level` = sig, `Estimated Power` = round(power_res$power, 3), Note = note, check.names = FALSE)
  } else if (k > 2) {
    total_n <- sum(group_sizes)
    power_res <- pwr::pwr.anova.test(k = k, n = total_n / k, f = effect_size, sig.level = sig)
    data.frame(Test = "ANOVA", `Effect Size (f)` = effect_size, Groups = paste(names(group_sizes), collapse = ", "), `Samples per Group` = round(total_n / k), `Significance Level` = sig, `Estimated Power` = round(power_res$power, 3), Note = note, check.names = FALSE)
  } else {
    data.frame(Message = "Not enough groups for power analysis", Note = note)
  }
}

make_power_curve_plot <- function(effect_size, test_type, n_min, n_max) {
  n_seq <- seq(n_min, n_max)
  sig <- 0.05
  power_vals <- sapply(n_seq, function(n) {
    if (test_type == "ttest") pwr::pwr.t.test(n = n, d = effect_size, sig.level = sig, type = "two.sample")$power
    else pwr::pwr.anova.test(k = 3, n = n, f = effect_size, sig.level = sig)$power
  })
  df <- data.frame(SampleSize = n_seq, Power = power_vals)
  plot_ly(df, x = ~SampleSize, y = ~Power, type = "scatter", mode = "lines+markers", line = list(color = "#00cc99", width = 3)) %>%
    layout(title = "Exploratory Generic Power Curve", xaxis = list(title = "Sample Size per Group"), yaxis = list(title = "Power", range = c(0, 1)),
           shapes = list(list(type = "line", x0 = min(n_seq), x1 = max(n_seq), y0 = 0.8, y1 = 0.8, line = list(dash = "dash", color = "red"))))
}

select_train_features_by_de <- function(train_counts, train_pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n) {
  resolved <- resolve_groups(train_pheno, phenotype_column, reference_group, contrast_arg)
  group <- resolved$group
  dge <- DGEList(counts = train_counts, group = group)
  design <- make_design(train_pheno, group, covariates)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 2) stop("Too few train-set genes after filterByExpr for RF feature selection.")
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  contrast_string <- paste0(make.names(resolved$numerator), "-", make.names(resolved$denominator))
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(v, design), contrast_matrix))
  res <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  res$Ensembl_IDs <- rownames(res)
  sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < fdr_cutoff & abs(res$logFC) >= logfc_cutoff, , drop = FALSE]
  if (nrow(sig) < 2) {
    log_msg("RF train-only feature selection found fewer than 2 significant genes; using top ", top_n, " train-set DE-ranked genes as exploratory features.")
    sig <- head(res[order(res$adj.P.Val), , drop = FALSE], top_n)
  }
  strip_ens_version(head(sig$Ensembl_IDs, top_n))
}

run_rf_classifier <- function(counts_raw, pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n) {
  pheno_vec <- as.factor(pheno[[phenotype_column]])
  if (nrow(pheno) < 6) stop("Too few samples for train/test classifier.")
  if (any(table(pheno_vec) < 2)) stop("Each phenotype class needs at least 2 samples for train/test classifier.")

  set.seed(seed)
  train_idx <- caret::createDataPartition(pheno_vec, p = 0.7, list = FALSE)
  train_samples <- rownames(pheno)[train_idx]
  test_samples <- setdiff(rownames(pheno), train_samples)
  if (length(test_samples) < 1) stop("Train/test split produced empty test set.")

  train_counts <- counts_raw[, train_samples, drop = FALSE]
  test_counts <- counts_raw[, test_samples, drop = FALSE]
  train_pheno <- pheno[train_samples, , drop = FALSE]
  test_pheno <- pheno[test_samples, , drop = FALSE]

  feature_ids <- select_train_features_by_de(train_counts, train_pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n)
  rownames(train_counts) <- strip_ens_version(rownames(train_counts))
  rownames(test_counts) <- strip_ens_version(rownames(test_counts))
  feature_ids <- intersect(feature_ids, intersect(rownames(train_counts), rownames(test_counts)))
  if (length(feature_ids) < 1) stop("No train-selected RF features matched both train and test count matrices.")

  combined_counts <- cbind(train_counts[feature_ids, , drop = FALSE], test_counts[feature_ids, , drop = FALSE])
  dge <- DGEList(counts = combined_counts)
  dge <- calcNormFactors(dge)
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
  ml_df <- data.frame(t(logcpm))
  ml_df$Phenotype <- as.factor(pheno[rownames(ml_df), phenotype_column])

  train_data <- ml_df[train_samples, , drop = FALSE]
  test_data <- ml_df[test_samples, , drop = FALSE]

  set.seed(seed)
  rf_model <- randomForest(Phenotype ~ ., data = train_data, ntree = 500, importance = TRUE)
  probs <- predict(rf_model, newdata = test_data, type = "prob")
  preds <- predict(rf_model, newdata = test_data)

  pred_table <- data.frame(Sample = rownames(test_data), Actual = test_data$Phenotype, Predicted = preds, Prob = apply(probs, 1, max), check.names = FALSE)
  cm <- caret::confusionMatrix(preds, test_data$Phenotype)

  pheno_levels <- levels(test_data$Phenotype)
  if (length(pheno_levels) == 2) {
    sens <- unname(cm$byClass["Sensitivity"])
    spec <- unname(cm$byClass["Specificity"])
    roc_obj <- pROC::roc(test_data$Phenotype, probs[, pheno_levels[2]], quiet = TRUE)
    auc_val <- round(as.numeric(roc_obj$auc), 3)
    roc_df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
    roc_plot <- plot_ly(data = roc_df, x = ~FPR, y = ~TPR, type = "scatter", mode = "lines", line = list(color = "#1f77b4", width = 2)) %>%
      layout(title = paste("Exploratory ROC Curve (AUC =", auc_val, ")"), xaxis = list(title = "False Positive Rate"), yaxis = list(title = "True Positive Rate"), showlegend = FALSE)
    metrics <- data.frame(Accuracy = unname(cm$overall["Accuracy"]), Sensitivity = round(sens, 3), Specificity = round(spec, 3), `AUC (ROC)` = auc_val, Note = "Exploratory classifier; feature selection performed on train set only.", check.names = FALSE)
  } else {
    by_class <- as.data.frame(cm$byClass)
    by_class$Class <- rownames(by_class)
    macro_sens <- mean(by_class$Sensitivity, na.rm = TRUE)
    macro_spec <- mean(by_class$Specificity, na.rm = TRUE)
    roc_plot <- plotly::plot_ly() %>% layout(title = "ROC only supported for 2-class problems")
    metrics <- data.frame(Accuracy = unname(cm$overall["Accuracy"]), MacroSensitivity = round(macro_sens, 3), MacroSpecificity = round(macro_spec, 3), `AUC (ROC)` = "N/A", Note = "Exploratory multi-class classifier; train-only feature selection.", check.names = FALSE)
  }

  imp <- randomForest::importance(rf_model, type = 1)
  imp_df <- data.frame(Feature = rownames(imp), MeanDecreaseAccuracy = as.numeric(imp[, 1]), stringsAsFactors = FALSE)
  imp_df <- imp_df[order(-imp_df$MeanDecreaseAccuracy), , drop = FALSE]

  list(predictions = pred_table, metrics = metrics, roc_plot = roc_plot, importance = imp_df)
}


train_classifier_model <- function(train_data, test_data, y_train, model_choice = "auto") {
  y_train <- droplevels(as.factor(y_train))
  model_choice <- tolower(model_choice)
  if (model_choice == "auto") model_choice <- "rf"

  if (model_choice == "logistic" && length(levels(y_train)) == 2 && nrow(train_data) > 3) {
    dat <- data.frame(y = as.numeric(y_train == levels(y_train)[2]), train_data, check.names = FALSE)
    fit <- tryCatch(stats::glm(y ~ ., data = dat, family = stats::binomial()), error = function(e) NULL)
    if (!is.null(fit)) {
      prob <- tryCatch(as.numeric(stats::predict(fit, newdata = test_data, type = "response")), error = function(e) NULL)
      if (!is.null(prob) && all(is.finite(prob))) {
        pred <- factor(ifelse(prob >= 0.5, levels(y_train)[2], levels(y_train)[1]), levels = levels(y_train))
        return(list(pred = pred, prob_pos = prob, model = "logistic", importance = data.frame(Feature = colnames(train_data), MeanDecreaseAccuracy = NA_real_)))
      }
    }
    log_msg("Logistic classifier failed or was unstable; falling back to random forest.")
  }

  rf_model <- randomForest(x = train_data, y = y_train, ntree = 500, importance = TRUE)
  pred <- predict(rf_model, newdata = test_data, type = "response")
  probs <- tryCatch(predict(rf_model, newdata = test_data, type = "prob"), error = function(e) NULL)
  prob_pos <- if (!is.null(probs) && ncol(probs) >= 2) as.numeric(probs[, levels(y_train)[2]]) else rep(NA_real_, nrow(test_data))
  imp <- randomForest::importance(rf_model, type = 1)
  imp_df <- data.frame(Feature = rownames(imp), MeanDecreaseAccuracy = as.numeric(imp[, 1]), stringsAsFactors = FALSE)
  list(pred = pred, prob_pos = prob_pos, model = "rf", importance = imp_df)
}

compute_classifier_outputs <- function(preds, importance_rows, pheno_levels, validation_label) {
  pred_factor <- factor(preds$Predicted, levels = pheno_levels)
  actual_factor <- factor(preds$Actual, levels = pheno_levels)
  cm <- caret::confusionMatrix(pred_factor, actual_factor)

  if (length(pheno_levels) == 2) {
    sens <- unname(cm$byClass["Sensitivity"])
    spec <- unname(cm$byClass["Specificity"])
    auc_val <- NA_real_
    roc_plot <- plotly::plot_ly() %>% layout(title = "ROC unavailable")
    if (!all(is.na(preds$Prob_Positive))) {
      roc_obj <- pROC::roc(actual_factor, preds$Prob_Positive, quiet = TRUE)
      auc_val <- round(as.numeric(roc_obj$auc), 3)
      roc_df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
      roc_plot <- plot_ly(data = roc_df, x = ~FPR, y = ~TPR, type = "scatter", mode = "lines", line = list(color = "#1f77b4", width = 2)) %>%
        layout(title = paste("Classifier ROC (AUC =", auc_val, ")"), xaxis = list(title = "False Positive Rate"), yaxis = list(title = "True Positive Rate"), showlegend = FALSE)
    }
    metrics <- data.frame(Accuracy = unname(cm$overall["Accuracy"]), Sensitivity = round(sens, 3), Specificity = round(spec, 3), `AUC (ROC)` = ifelse(is.na(auc_val), "N/A", auc_val), Validation = validation_label, Note = "Exploratory classifier; feature selection performed inside training split/fold where applicable.", check.names = FALSE)
  } else {
    by_class <- as.data.frame(cm$byClass)
    metrics <- data.frame(Accuracy = unname(cm$overall["Accuracy"]), MacroSensitivity = round(mean(by_class$Sensitivity, na.rm = TRUE), 3), MacroSpecificity = round(mean(by_class$Specificity, na.rm = TRUE), 3), `AUC (ROC)` = "N/A", Validation = validation_label, Note = "Exploratory multi-class classifier.", check.names = FALSE)
    roc_plot <- plotly::plot_ly() %>% layout(title = "ROC only supported for 2-class problems")
  }

  if (length(importance_rows) > 0) {
    imp_all <- do.call(rbind, importance_rows)
    imp_all$Rank <- ave(-imp_all$MeanDecreaseAccuracy, imp_all$Fold, FUN = rank, ties.method = "first")
    stability <- aggregate(list(selection_count = imp_all$Feature), by = list(Feature = imp_all$Feature), FUN = length)
    stability$selection_frequency <- stability$selection_count / length(unique(imp_all$Fold))
    mean_imp <- aggregate(MeanDecreaseAccuracy ~ Feature, imp_all, mean, na.rm = TRUE)
    mean_rank <- aggregate(Rank ~ Feature, imp_all, mean, na.rm = TRUE)
    stability <- merge(stability, mean_imp, by = "Feature", all.x = TRUE)
    stability <- merge(stability, mean_rank, by = "Feature", all.x = TRUE)
    stability <- stability[order(-stability$selection_frequency, stability$Rank), , drop = FALSE]
  } else {
    imp_all <- data.frame()
    stability <- data.frame()
  }

  list(predictions = preds, metrics = metrics, roc_plot = roc_plot, importance = imp_all, stability = stability)
}

run_cv_classifier <- function(counts_raw, pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n, model_choice = "auto") {
  pheno_vec <- droplevels(as.factor(pheno[[phenotype_column]]))
  if (nrow(pheno) < 4) stop("Too few samples for CV classifier.")
  if (any(table(pheno_vec) < 2)) stop("Each phenotype class needs at least 2 samples for CV classifier.")
  pheno_levels <- levels(pheno_vec)
  k <- min(5, as.integer(min(table(pheno_vec))))
  k <- max(2, k)
  set.seed(seed)
  folds <- caret::createFolds(pheno_vec, k = k, returnTrain = FALSE)

  prediction_rows <- list()
  importance_rows <- list()
  for (fold_i in seq_along(folds)) {
    test_samples <- rownames(pheno)[folds[[fold_i]]]
    train_samples <- setdiff(rownames(pheno), test_samples)
    train_counts <- counts_raw[, train_samples, drop = FALSE]
    test_counts <- counts_raw[, test_samples, drop = FALSE]
    train_pheno <- pheno[train_samples, , drop = FALSE]

    feature_ids <- tryCatch(select_train_features_by_de(train_counts, train_pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n), error = function(e) {
      log_msg("Fold ", fold_i, " DE feature selection failed: ", conditionMessage(e))
      character(0)
    })
    rownames(train_counts) <- strip_ens_version(rownames(train_counts))
    rownames(test_counts) <- strip_ens_version(rownames(test_counts))
    feature_ids <- intersect(feature_ids, intersect(rownames(train_counts), rownames(test_counts)))
    if (length(feature_ids) < 2) {
      log_msg("Fold ", fold_i, " skipped: fewer than 2 train-selected features.")
      next
    }

    combined_counts <- cbind(train_counts[feature_ids, , drop = FALSE], test_counts[feature_ids, , drop = FALSE])
    dge <- DGEList(counts = combined_counts)
    dge <- calcNormFactors(dge)
    logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
    x_train <- data.frame(t(logcpm[, train_samples, drop = FALSE]), check.names = FALSE)
    x_test <- data.frame(t(logcpm[, test_samples, drop = FALSE]), check.names = FALSE)
    y_train <- droplevels(as.factor(pheno[train_samples, phenotype_column]))
    actual <- as.factor(pheno[test_samples, phenotype_column])

    pred_obj <- train_classifier_model(x_train, x_test, y_train, model_choice = model_choice)
    prediction_rows[[length(prediction_rows) + 1]] <- data.frame(Fold = fold_i, Sample = test_samples, Actual = as.character(actual), Predicted = as.character(pred_obj$pred), Prob_Positive = pred_obj$prob_pos, Model = pred_obj$model, Validation = "cv", stringsAsFactors = FALSE)
    imp <- pred_obj$importance
    if (nrow(imp) > 0) {
      imp$Fold <- fold_i
      imp$Model <- pred_obj$model
      importance_rows[[length(importance_rows) + 1]] <- imp
    }
  }

  if (length(prediction_rows) == 0) stop("No CV folds produced classifier predictions.")
  preds <- do.call(rbind, prediction_rows)
  compute_classifier_outputs(preds, importance_rows, pheno_levels, validation_label = paste0("cv_", length(unique(preds$Fold)), "fold"))
}

main <- function() {
  set.seed(seed)
  run_started <- Sys.time()
  write_run_status("running", run_started)
  date_tag <- as.character(Sys.Date())
  log_msg("Starting JCAP RNA-SEQ Analyzer CLI - production statistics edition")
  log_msg("Counts: ", counts_path)
  log_msg("Phenotype: ", phenotype_path)
  log_msg("Phenotype column: ", phenotype_column)
  log_msg("Species: ", species_code)
  log_msg("Output directory: ", outdir)
  log_msg("FDR cutoff: ", fdr_cutoff)
  log_msg("logFC cutoff: ", logfc_cutoff)
  if (length(covariates) > 0) log_msg("Covariates: ", paste(covariates, collapse = ", "))

  if (!file.exists(counts_path)) stop("Counts file does not exist: ", counts_path)
  if (!file.exists(phenotype_path)) stop("Phenotype file does not exist: ", phenotype_path)

  log_msg("Reading input files...")
  counts <- read_matrix_file(counts_path)
  pheno <- read_matrix_file(phenotype_path)

  if (!phenotype_column %in% colnames(pheno)) stop("Phenotype column not found. Available columns: ", paste(colnames(pheno), collapse = ", "))
  aligned <- align_counts_pheno(counts, pheno)
  counts <- aligned$counts
  pheno <- aligned$pheno

  log_msg("Matched samples: ", ncol(counts))
  log_msg("Input genes: ", nrow(counts))

  prepared <- prepare_counts(counts, pheno, phenotype_column, date_tag)
  counts_final <- prepared$counts
  pheno_final <- prepared$pheno

  log_msg("Running limma/voom differential expression with explicit contrast...")
  de <- run_differential_expression(
    counts_final, pheno_final, phenotype_column, species_code,
    covariates = covariates,
    reference_group = reference_group,
    contrast_arg = contrast_arg,
    skip_annotation = skip_annotation
  )

  final_results <- de$de_sig
  de_all <- de$de_all
  logcpm <- de$logcpm

  saveRDS(de$dge, file.path(objects_dir, "dge.rds"))
  saveRDS(de$voom_object, file.path(objects_dir, "voom_object.rds"))
  saveRDS(de$limma_fit, file.path(objects_dir, "limma_fit.rds"))
  saveRDS(de$contrast_matrix, file.path(objects_dir, "contrast_matrix.rds"))
  log_msg("Wrote RDS analysis objects to: ", objects_dir)

  de_path <- file.path(tables_dir, paste0("DE_results_", date_tag, ".csv"))
  write.csv(final_results, de_path, row.names = FALSE)
  log_msg("Wrote significant DE table: ", de_path, " (", nrow(final_results), " genes)")

  limma_path <- file.path(tables_dir, paste0("limma_all_results_", date_tag, ".csv"))
  write.csv(de_all, limma_path, row.names = FALSE)
  log_msg("Wrote full limma table: ", limma_path)

  logcpm_path <- file.path(tables_dir, paste0("normalized_logCPM_", date_tag, ".csv"))
  write.csv(logcpm, logcpm_path)
  log_msg("Wrote normalized logCPM matrix: ", logcpm_path)

  if (!skip_plots) {
    log_msg("Creating PCA plot from normalized logCPM...")
    tryCatch(save_plotly_html(make_pca_plot(logcpm, pheno_final, phenotype_column), file.path(plots_dir, "pca_plot.html")), error = function(e) log_msg("PCA plot failed: ", conditionMessage(e)))

    log_msg("Creating UMAP plot from normalized logCPM...")
    tryCatch(save_plotly_html(make_umap_plot(logcpm, pheno_final, phenotype_column), file.path(plots_dir, "umap_plot.html")), error = function(e) log_msg("UMAP plot failed: ", conditionMessage(e)))

    log_msg("Creating volcano plot from full DE table...")
    tryCatch(save_plotly_html(make_volcano_plot(de_all), file.path(plots_dir, "volcano_plot.html")), error = function(e) log_msg("Volcano plot failed: ", conditionMessage(e)))

    log_msg("Creating heatmap plot from normalized logCPM...")
    tryCatch(save_plotly_html(make_heatmap_plot(logcpm, final_results, de_all), file.path(plots_dir, "heatmap_plot.html")), error = function(e) log_msg("Heatmap plot failed: ", conditionMessage(e)))
  }

  enrichment_counts <- list()
  if (!skip_enrichment) {
    enrich_db <- if (!is.null(args[["enrich-db"]])) args[["enrich-db"]] else default_enrich_db(species_code)
    log_msg("Using enrichment database/source: ", enrich_db)

    if (!"symbol" %in% names(final_results)) {
      final_results$symbol <- final_results$Ensembl_IDs
    }

    enrichment_sets <- list(
      all = list(genes = na.omit(final_results$symbol), csv = "enrichment_all.csv", plot = "enrichment_all_plot.html", title = "Significant DE Genes"),
      upregulated = list(genes = na.omit(final_results$symbol[final_results$logFC >= logfc_cutoff & final_results$adj.P.Val < fdr_cutoff]), csv = "enrichment_upregulated.csv", plot = "enrichment_upregulated_plot.html", title = "Upregulated Significant Genes"),
      downregulated = list(genes = na.omit(final_results$symbol[final_results$logFC <= -logfc_cutoff & final_results$adj.P.Val < fdr_cutoff]), csv = "enrichment_downregulated.csv", plot = "enrichment_downregulated_plot.html", title = "Downregulated Significant Genes")
    )

    for (nm in names(enrichment_sets)) {
      item <- enrichment_sets[[nm]]
      log_msg("Running enrichment: ", nm)
      res <- perform_enrichment_v6(item$genes, enrich_db, species_code, backend = enrich_backend, gprof_sources = gprof_sources, enrichr_db = enrichr_db)
      csv_path <- file.path(tables_dir, item$csv)
      network_csv <- sub("\\.csv$", "_network_edges.csv", item$csv)
      network_path <- file.path(tables_dir, network_csv)

      if (is.null(res) || nrow(res) == 0) {
        write.csv(data.frame(Message = paste("No enrichment results for", nm, "using significant DE genes")), csv_path, row.names = FALSE)
        write.csv(make_network_edges(NULL), network_path, row.names = FALSE)
        enrichment_counts[[nm]] <- 0
        log_msg("No enrichment results for ", nm)
        log_msg("Wrote empty network-ready enrichment edges: ", network_path)
      } else {
        write.csv(res, csv_path, row.names = FALSE)
        enrichment_counts[[nm]] <- nrow(res)
        log_msg("Wrote: ", csv_path)
        network_edges <- tryCatch(make_network_edges(res), error = function(e) {
          log_msg("Network-ready enrichment edge export failed for ", nm, ": ", conditionMessage(e))
          data.frame(gene = character(), term = character(), adjusted_pvalue = numeric(), stringsAsFactors = FALSE)
        })
        write.csv(network_edges, network_path, row.names = FALSE)
        log_msg("Wrote network-ready enrichment edges: ", network_path)
        if (!skip_plots) {
          plot_obj <- tryCatch(plot_enrichment_bar(res, item$title), error = function(e) {
            log_msg("Enrichment plot failed for ", nm, ": ", conditionMessage(e))
            NULL
          })
          if (!is.null(plot_obj)) save_plotly_html(plot_obj, file.path(plots_dir, item$plot))
        }
      }
    }
  }

  log_msg("Running generic exploratory power analysis...")
  power_summary <- tryCatch(run_power_analysis(pheno_final, phenotype_column, effect_size), error = function(e) data.frame(Error = conditionMessage(e)))
  power_path <- file.path(tables_dir, paste0("exploratory_power_summary_", date_tag, ".csv"))
  write.csv(power_summary, power_path, row.names = FALSE)
  log_msg("Wrote: ", power_path)

  if (!skip_plots) {
    log_msg("Creating exploratory power curve plot...")
    tryCatch(save_plotly_html(make_power_curve_plot(effect_size, power_test_type, curve_n_min, curve_n_max), file.path(plots_dir, "power_curve_plot.html")), error = function(e) log_msg("Power curve plot failed: ", conditionMessage(e)))
  }

  rf_metrics_for_summary <- NULL
  if (!skip_classifier && classifier_validation != "none") {
    log_msg("Running exploratory classifier with validation mode: ", classifier_validation)
    rf <- tryCatch({
      if (classifier_validation == "cv") {
        run_cv_classifier(counts_final, pheno_final, phenotype_column, reference_group, contrast_arg, covariates, rf_top_n, model_choice = classifier_model)
      } else {
        run_rf_classifier(counts_final, pheno_final, phenotype_column, reference_group, contrast_arg, covariates, rf_top_n)
      }
    }, error = function(e) {
      log_msg("Classifier failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(rf)) {
      pred_path <- file.path(tables_dir, paste0("rf_predictions_", date_tag, ".csv"))
      metric_path <- file.path(tables_dir, paste0("rf_metrics_", date_tag, ".csv"))
      imp_path <- file.path(tables_dir, paste0("rf_feature_importance_", date_tag, ".csv"))
      stability_path <- file.path(tables_dir, paste0("classifier_feature_stability_", date_tag, ".csv"))
      write.csv(rf$predictions, pred_path, row.names = FALSE)
      write.csv(rf$metrics, metric_path, row.names = FALSE)
      write.csv(rf$importance, imp_path, row.names = FALSE)
      if (!is.null(rf$stability)) write.csv(rf$stability, stability_path, row.names = FALSE)
      rf_metrics_for_summary <- rf$metrics
      log_msg("Wrote: ", pred_path)
      log_msg("Wrote: ", metric_path)
      log_msg("Wrote: ", imp_path)
      if (!is.null(rf$stability)) log_msg("Wrote: ", stability_path)
      if (!skip_plots) save_plotly_html(rf$roc_plot, file.path(plots_dir, "rf_roc_plot.html"))
    }
  } else if (classifier_validation == "none") {
    log_msg("Skipping classifier metrics because --classifier-validation none was selected.")
  }

  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()

  write_run_summary(de_all = de_all, final_results = final_results, enrichment_counts = enrichment_counts, rf_metrics = rf_metrics_for_summary, matched_samples = nrow(pheno_final))
  write_json_file(make_manifest(date_tag), file.path(outdir, "manifest.json"))

  zip_path <- NULL
  if (zip_output) {
    log_msg("Zipping output directory...")
    zip_path <- tryCatch(zip_results(), error = function(e) {
      log_msg("Zip failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(zip_path)) log_msg("Wrote zip bundle: ", zip_path)
  }

  write_run_status("complete", run_started, finished = Sys.time(), zip_path = zip_path)
  log_msg("Done.")
}

tryCatch(main(), error = function(e) {
  log_msg("FATAL: ", conditionMessage(e))
  write_run_status("failed", Sys.time(), finished = Sys.time(), error = conditionMessage(e))
  quit(status = 1)
})

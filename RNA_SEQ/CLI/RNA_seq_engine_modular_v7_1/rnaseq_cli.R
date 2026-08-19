#!/usr/bin/env Rscript
# RNA_SEQ_CLI_V7_STAT_CLEARANCE_CANDIDATE_2026_08_10
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

# Resolve modules relative to this script, while leaving user input paths relative to the current working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}

module_files <- c(
  "validation.R",
  "reporting.R",
  "io.R",
  "annotation.R",
  "design_de.R",
  "enrichment.R",
  "plots.R",
  "power.R",
  "classifier.R"
)
for (module_file in module_files) {
  source(file.path(script_dir, "R", module_file), local = .GlobalEnv)
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
  enrichment_background_info <- NULL
  if (!skip_enrichment) {
    enrich_db <- if (!is.null(args[["enrich-db"]])) args[["enrich-db"]] else default_enrich_db(species_code)
    log_msg("Using enrichment database/source: ", enrich_db)

    if (!"symbol" %in% names(final_results)) final_results$symbol <- final_results$Ensembl_IDs
    if (!"symbol" %in% names(de_all)) de_all$symbol <- de_all$Ensembl_IDs
    background_genes <- unique(na.omit(as.character(de_all$symbol)))
    background_genes <- background_genes[nzchar(background_genes)]
    background_path <- file.path(tables_dir, "enrichment_background_tested_genes.csv")
    write.csv(data.frame(gene = background_genes, stringsAsFactors = FALSE), background_path, row.names = FALSE)
    enrichment_background_info <- list(
      definition = "Genes surviving edgeR::filterByExpr and therefore eligible for differential-expression testing; represented in the same identifier field used for enrichment.",
      n_genes = length(background_genes),
      file = "tables/enrichment_background_tested_genes.csv"
    )
    log_msg("Enrichment background: ", length(background_genes), " DE-tested genes; wrote ", background_path)

    enrichment_sets <- list(
      all = list(genes = na.omit(final_results$symbol), csv = "enrichment_all.csv", plot = "enrichment_all_plot.html", title = "Significant DE Genes"),
      upregulated = list(genes = na.omit(final_results$symbol[final_results$logFC >= logfc_cutoff & final_results$adj.P.Val < fdr_cutoff]), csv = "enrichment_upregulated.csv", plot = "enrichment_upregulated_plot.html", title = "Upregulated Significant Genes"),
      downregulated = list(genes = na.omit(final_results$symbol[final_results$logFC <= -logfc_cutoff & final_results$adj.P.Val < fdr_cutoff]), csv = "enrichment_downregulated.csv", plot = "enrichment_downregulated_plot.html", title = "Downregulated Significant Genes")
    )

    for (nm in names(enrichment_sets)) {
      item <- enrichment_sets[[nm]]
      log_msg("Running enrichment: ", nm)
      res <- perform_enrichment_v6(item$genes, enrich_db, species_code, backend = enrich_backend, gprof_sources = gprof_sources, enrichr_db = enrichr_db, background_genes = background_genes)
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

  log_msg("Running generic effect-size power reference (not RNA-seq-specific)...")
  power_summary <- tryCatch(run_power_analysis(pheno_final, phenotype_column, effect_size), error = function(e) data.frame(Error = conditionMessage(e)))
  power_path <- file.path(tables_dir, paste0("generic_effect_size_power_reference_", date_tag, ".csv"))
  write.csv(power_summary, power_path, row.names = FALSE)
  log_msg("Wrote: ", power_path)

  if (!skip_plots) {
    log_msg("Creating generic effect-size power reference curve...")
    tryCatch(save_plotly_html(make_power_curve_plot(effect_size, power_test_type, curve_n_min, curve_n_max), file.path(plots_dir, "generic_effect_size_power_reference_curve.html")), error = function(e) log_msg("Power curve plot failed: ", conditionMessage(e)))
  }

  rf_metrics_for_summary <- NULL
  if (!skip_classifier && classifier_validation != "none") {
    log_msg("Running exploratory classifier with validation mode: ", classifier_validation)
    rf <- tryCatch({
      if (classifier_validation == "cv") {
        run_cv_classifier(counts_final, pheno_final, phenotype_column, reference_group, contrast_arg, covariates, rf_top_n, model_choice = classifier_model)
      } else {
        run_rf_classifier(counts_final, pheno_final, phenotype_column, reference_group, contrast_arg, covariates, rf_top_n, model_choice = classifier_model)
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

  write_run_summary(de_all = de_all, final_results = final_results, enrichment_counts = enrichment_counts, rf_metrics = rf_metrics_for_summary, matched_samples = nrow(pheno_final), enrichment_background = enrichment_background_info)
  write_json_file(make_manifest(date_tag, enrichment_background = enrichment_background_info), file.path(outdir, "manifest.json"))

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

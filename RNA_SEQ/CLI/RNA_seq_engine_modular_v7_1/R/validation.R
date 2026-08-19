# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

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
  tables/generic_effect_size_power_reference_<date>.csv
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


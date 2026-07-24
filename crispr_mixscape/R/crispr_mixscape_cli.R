#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

# Resolve source path whether the CLI is run from the module root or elsewhere.
this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
core_candidates <- c(
  file.path(getwd(), "R", "crispr_mixscape_core.R"),
  file.path(dirname(this_file), "crispr_mixscape_core.R")
)
core_file <- core_candidates[file.exists(core_candidates)][1]
if (is.na(core_file)) stop("Could not find R/crispr_mixscape_core.R")
source(core_file)

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(paste("Unexpected argument:", key))
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

a <- parse_args(args)

required <- c("counts", "metadata", "outdir")
missing <- setdiff(required, names(a))
if (length(missing)) stop(paste("Missing required args:", paste(missing, collapse = ", ")))

dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  manifest <- run_pipeline(
    counts_csv = a$counts,
    metadata_csv = a$metadata,
    outdir = a$outdir,
    species = a$species %||% "hsapiens",
    nt_label = a$nt_label %||% "NT",
    ko_label = a$ko_label %||% NULL,
    num_neighbors = as.integer(a$neighbors %||% 20),
    max_pcs = as.integer(a$max_pcs %||% 40),
    min_group_cells = as.integer(a$min_group_cells %||% 3),
    min_de_genes = as.integer(a$min_de_genes %||% 3),
    iter_num = as.integer(a$iter_num %||% 20),
    de_mode = a$de_mode %||% "auto",
    pseudobulk_method = a$pseudobulk_method %||% "limma_voom",
    sample_col = a$sample_col %||% "replicate",
    pseudobulk_min_replicates = as.integer(a$pseudobulk_min_replicates %||% 2),
    de_fdr = as.numeric(a$de_fdr %||% 0.05),
    logfc_threshold = as.numeric(a$logfc_threshold %||% 0),
    min_pct = as.numeric(a$min_pct %||% 0.1),
    test_use = a$test_use %||% "wilcox",
    enrich_backend = a$enrich_backend %||% "auto",
    gprof_sources = a$gprof_sources %||% "GO:BP,GO:MF,GO:CC,REAC,KEGG",
    enrichr_dbs = a$enrichr_dbs %||% "GO_Biological_Process_2021,GO_Molecular_Function_2021,GO_Cellular_Component_2021,Reactome_2022,KEGG_2021_Human",
    skip_enrichment = as_bool(a$skip_enrichment, FALSE),
    skip_plots = as_bool(a$skip_plots, FALSE)
  )
  invisible(manifest)
}, error = function(e) {
  status_path <- file.path(a$outdir, "run_status.json")
  if (!file.exists(status_path)) {
    writeLines(jsonlite::toJSON(list(status = "failed", message = e$message, updated_at = as.character(Sys.time())), auto_unbox = TRUE, pretty = TRUE), status_path)
  }
  stop(e)
})

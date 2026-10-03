#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

this_file <- tryCatch(
  normalizePath(sys.frames()[[1]]$ofile),
  error = function(e) NA_character_
)

core_candidates <- c(
  file.path(getwd(), "R", "crispr_mixscape_core.R"),
  file.path(dirname(this_file), "crispr_mixscape_core.R")
)

core_file <- core_candidates[
  file.exists(core_candidates)
][1]

if (is.na(core_file)) {
  stop("Could not find R/crispr_mixscape_core.R")
}

source(core_file)

parse_args <- function(args) {

  out <- list()
  i <- 1L

  while (i <= length(args)) {

    key <- args[[i]]

    if (!startsWith(key, "--")) {
      stop(
        paste(
          "Unexpected argument:",
          key
        )
      )
    }

    name <- sub(
      "^--",
      "",
      key
    )

    if (
      i == length(args) ||
      startsWith(
        args[[i + 1L]],
        "--"
      )
    ) {
      out[[name]] <- TRUE
      i <- i + 1L
    } else {
      out[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }

  out
}

args <- parse_args(
  commandArgs(
    trailingOnly = TRUE
  )
)

required <- c(
  "rds",
  "outdir",
  "ko_label"
)

missing <- setdiff(
  required,
  names(args)
)

if (length(missing)) {
  stop(
    "Missing required args: ",
    paste(missing, collapse = ", ")
  )
}

dir.create(
  args$outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

species <- species_map(
  args$species %||% "mmusculus"
)

nt_label <- args$nt_label %||% "CTRL"
ko_label <- args$ko_label
cell_type <- args$cell_type %||% NULL

de_mode <- args$de_mode %||% "auto"

pseudobulk_method <-
  args$pseudobulk_method %||%
  "limma_voom"

sample_col <-
  args$sample_col %||%
  "replicate"

pseudobulk_min_replicates <-
  as.integer(
    args$pseudobulk_min_replicates %||%
    2
  )

de_fdr <- as.numeric(
  args$de_fdr %||% 0.05
)

logfc_threshold <- as.numeric(
  args$logfc_threshold %||% 0
)

min_pct <- as.numeric(
  args$min_pct %||% 0.1
)

test_use <-
  args$test_use %||%
  "wilcox"

enrich_backend <-
  args$enrich_backend %||%
  "auto"

gprof_sources <-
  args$gprof_sources %||%
  "GO:BP,GO:MF,GO:CC,REAC,KEGG"

enrichr_dbs <-
  args$enrichr_dbs %||%
  paste(
    c(
      "GO_Biological_Process_2021",
      "GO_Molecular_Function_2021",
      "GO_Cellular_Component_2021",
      "Reactome_2022",
      "KEGG_2021_Human"
    ),
    collapse = ","
  )

parameters <- list(
  stage = "contrast",
  species = species,
  cell_type = cell_type,
  nt_label = nt_label,
  ko_label = ko_label,
  de_mode = de_mode,
  pseudobulk_method = pseudobulk_method,
  sample_col = sample_col,
  pseudobulk_min_replicates =
    pseudobulk_min_replicates,
  de_fdr = de_fdr,
  logfc_threshold = logfc_threshold,
  min_pct = min_pct,
  test_use = test_use,
  enrich_backend = enrich_backend,
  gprof_sources = gprof_sources,
  enrichr_dbs = enrichr_dbs
)

init_status(
  args$outdir,
  parameters
)

tryCatch({

  so <- record_step(
    "load_classified_object",
    {
      readRDS(args$rds)
    }
  )

  resolved_ko <- resolve_ko_label(
    so,
    ko_label
  )

  message(
    "Requested KO: ",
    ko_label,
    "; resolved Mixscape class: ",
    resolved_ko
  )

  de_tab <- run_de(
    so,
    ko_label = ko_label,
    nt_label = nt_label,
    de_mode = de_mode,
    pseudobulk_method =
      pseudobulk_method,
    sample_col = sample_col,
    pseudobulk_min_replicates =
      pseudobulk_min_replicates,
    min_pct = min_pct,
    logfc_threshold =
      logfc_threshold,
    test_use = test_use
  )

  manifest <- write_bundle(
    so,
    de_tab,
    outdir = args$outdir,
    species = species,
    parameters = parameters,
    skip_enrichment = FALSE,
    skip_plots = FALSE
  )

  write_crispr_literature_handoff(
    so,
    de_tab,
    outdir = args$outdir,
    species = species,
    parameters = parameters
  )

  provenance <- list(
    requested_perturbation =
      ko_label,
    resolved_mixscape_class =
      resolved_ko,
    de_comparison =
      if (
        "comparison" %in%
        colnames(de_tab) &&
        nrow(de_tab) > 0
      ) {
        as.character(
          de_tab$comparison[[1]]
        )
      } else {
        NA_character_
      }
  )

  jsonlite::write_json(
    provenance,
    file.path(
      args$outdir,
      "target_provenance.json"
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  .run_status$x$status <- "complete"
  .run_status$x$finished <-
    as.character(
      Sys.time()
    )

  write_status()

  invisible(manifest)

}, error = function(e) {

  if (!is.null(.run_status$x)) {

    .run_status$x$status <- "failed"

    .run_status$x$finished <-
      as.character(
        Sys.time()
      )

    .run_status$x$errors[[length(.run_status$x$errors) + 1L]] <- list(
      step = "contrast",
      message = conditionMessage(e)
    )

    write_status()
  }

  stop(e)
})

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
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

core_file <- core_candidates[file.exists(core_candidates)][1]

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
      stop(paste("Unexpected argument:", key))
    }

    name <- sub("^--", "", key)

    if (
      i == length(args) ||
      startsWith(args[[i + 1L]], "--")
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

read_one_col <- function(path, expected_name) {
  x <- read.delim(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) != 1L) {
    stop(path, " must contain exactly one column")
  }

  names(x)[1] <- expected_name
  as.character(x[[1]])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c(
  "counts",
  "cells",
  "genes",
  "metadata",
  "outdir"
)

missing <- setdiff(required, names(args))

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

num_neighbors <- as.integer(
  args$neighbors %||% 20
)

max_pcs <- as.integer(
  args$max_pcs %||% 40
)

min_group_cells <- as.integer(
  args$min_group_cells %||% 3
)

min_de_genes <- as.integer(
  args$min_de_genes %||% 3
)

iter_num <- as.integer(
  args$iter_num %||% 20
)

parameters <- list(
  stage = "classification",
  species = species,
  nt_label = nt_label,
  neighbors = num_neighbors,
  max_pcs = max_pcs,
  min_group_cells = min_group_cells,
  min_de_genes = min_de_genes,
  iter_num = iter_num
)

init_status(
  args$outdir,
  parameters
)

tryCatch({

  x <- record_step(
    "load_giladi_inputs",
    {

      cells <- read_one_col(
        args$cells,
        "cell_id"
      )

      genes <- read_one_col(
        args$genes,
        "gene"
      )

      counts <- Matrix::readMM(
        args$counts
      )

      counts <- as(
        counts,
        "dgCMatrix"
      )

      if (nrow(counts) != length(genes)) {
        stop(
          sprintf(
            "Matrix rows (%d) do not match genes (%d)",
            nrow(counts),
            length(genes)
          )
        )
      }

      if (ncol(counts) != length(cells)) {
        stop(
          sprintf(
            "Matrix columns (%d) do not match cells (%d)",
            ncol(counts),
            length(cells)
          )
        )
      }

      if (anyDuplicated(genes)) {
        stop("Gene table contains duplicated gene names.")
      }

      if (anyDuplicated(cells)) {
        stop("Cell table contains duplicated cell IDs.")
      }

      rownames(counts) <- genes
      colnames(counts) <- cells

      meta <- read.csv(
        args$metadata,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      needed <- c(
        "cell_id",
        "target",
        "guide",
        "sample"
      )

      miss <- setdiff(
        needed,
        colnames(meta)
      )

      if (length(miss)) {
        stop(
          "Metadata is missing required columns: ",
          paste(miss, collapse = ", ")
        )
      }

      if (anyDuplicated(meta$cell_id)) {
        stop(
          "Metadata contains duplicated cell_id values."
        )
      }

      idx <- match(
        cells,
        meta$cell_id
      )

      if (anyNA(idx)) {
        stop(
          sum(is.na(idx)),
          " count-matrix cells are absent from metadata."
        )
      }

      meta <- meta[
        idx,
        ,
        drop = FALSE
      ]

      rownames(meta) <- meta$cell_id

      meta$gene <- as.character(
        meta$target
      )

      meta$replicate <- as.character(
        meta$sample
      )

      meta$guide_ID <- as.character(
        meta$guide
      )

      meta$crispr <- as.character(
        meta$target
      )

      if (
        !any(
          toupper(meta$gene) ==
          toupper(nt_label)
        )
      ) {
        stop(
          "No control cells matching nt_label='",
          nt_label,
          "' were found."
        )
      }

      list(
        counts = counts,
        metadata = meta
      )
    }
  )

  write_input_qc(
    x$counts,
    x$metadata,
    args$outdir
  )

  so <- record_step(
    "normalize_pca_umap",
    {

      obj <- CreateSeuratObject(
        counts = x$counts,
        meta.data = x$metadata
      )

      obj@misc$input_n_cells <- ncol(x$counts)
      attr(obj, "input_n_cells") <- ncol(x$counts)

      DefaultAssay(obj) <- "RNA"

      obj <- .ensure_counts_layer(
        obj,
        "RNA"
      )

      obj <- NormalizeData(
        obj,
        verbose = FALSE
      ) |>
        FindVariableFeatures(
          verbose = FALSE
        ) |>
        ScaleData(
          verbose = FALSE
        )

      n_cells <- ncol(obj)
      n_genes <- nrow(obj)

      npcs <- min(
        max_pcs,
        n_cells - 1L,
        n_genes - 1L
      )

      if (npcs < 2L) {
        stop(
          "Not enough genes/cells for PCA."
        )
      }

      obj <- RunPCA(
        obj,
        npcs = npcs,
        verbose = FALSE
      )

      umap_neighbors <- min(
        num_neighbors,
        n_cells - 1L
      )

      if (umap_neighbors < 2L) {
        umap_neighbors <- 2L
      }

      obj <- RunUMAP(
        obj,
        dims = seq_len(npcs),
        n.neighbors = umap_neighbors,
        verbose = FALSE
      )

      attr(obj, "npcs") <- npcs
      attr(obj, "umap_neighbors") <- umap_neighbors

      obj
    }
  )

  so <- run_mixscape(
    so,
    num_neighbors = num_neighbors,
    min_de_genes = min_de_genes,
    iter_num = iter_num,
    nt_class = nt_label,
    min_group_cells = min_group_cells,
    max_pcs = max_pcs
  )

  saveRDS(
    so,
    file.path(
      args$outdir,
      "mixscape_classified.rds"
    )
  )

  write_class_tables(
    so,
    args$outdir
  )

  stats <- summary_stats(so)

  write.csv(
    stats,
    file.path(
      args$outdir,
      "classification_summary.csv"
    ),
    row.names = FALSE
  )

  classes <- ko_labels(so)

  jsonlite::write_json(
    list(
      status = "complete",
      module = "CRISPR Mixscape classification",
      stage = "classification",
      available_ko_classes = classes,
      files = list.files(
        args$outdir,
        recursive = TRUE
      )
    ),
    file.path(
      args$outdir,
      "manifest.json"
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  .run_status$x$status <- "complete"
  .run_status$x$finished <- as.character(
    Sys.time()
  )

  write_status()

}, error = function(e) {

  if (!is.null(.run_status$x)) {

    .run_status$x$status <- "failed"

    .run_status$x$finished <- as.character(
      Sys.time()
    )

    .run_status$x$errors[[length(.run_status$x$errors) + 1L]] <- list(
      step = "classification",
      message = conditionMessage(e)
    )

    write_status()
  }

  stop(e)
})

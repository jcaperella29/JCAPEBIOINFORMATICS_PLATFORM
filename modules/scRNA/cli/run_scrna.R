#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(Seurat)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(gprofiler2)
  library(enrichR)
  library(randomForest)
  library(caret)
  library(pROC)
  library(pwr)
})

VERSION <- "0.4.0"

option_list <- list(
  make_option(c("--counts"), type="character", help="Counts matrix CSV. Genes as rows, cells as columns."),
  make_option(c("--metadata"), type="character", help="Metadata CSV. Cells as rows. Must include condition/cell type columns for full pipeline."),
  make_option(c("--outdir"), type="character", default="scrna_results", help="Output directory [default %default]"),
  make_option(c("--species"), type="character", default="hsapiens", help="hsapiens/mmusculus/dmelanogaster/drerio or Human/Mouse/Fly/Zebrafish"),
  make_option(c("--id-type"), type="character", default="Symbol", help="Display only: Symbol/ENSEMBL/Entrez"),
  make_option(c("--condition-col"), type="character", default="stim", help="Metadata condition column [default %default]"),
  make_option(c("--celltype-col"), type="character", default="cell_type", help="Metadata cell type column [default %default]"),
  make_option(c("--annotation-col"), type="character", default="seurat_annotations", help="Preferred metadata annotation column for plots [default %default]"),
  make_option(c("--de-fdr"), type="double", default=0.05, help="Adjusted p-value cutoff for condition-only DE filtering [default %default]"),
  make_option(c("--logfc-threshold"), type="double", default=0.0, help="Minimum absolute log2FC for DE filtering [default %default]"),
  make_option(c("--min-pct"), type="double", default=0.1, help="Seurat FindMarkers min.pct [default %default]"),
  make_option(c("--test-use"), type="character", default="wilcox", help="Seurat FindMarkers test.use. Common choices: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2 [default %default]"),
  make_option(c("--de-mode"), type="character", default="cell", help="DE mode: cell or pseudobulk [default %default]"),
  make_option(c("--de-scope"), type="character", default="top_celltype", help="Condition DE scope: global, top_celltype, or target_celltype [default %default]"),
  make_option(c("--condition-a"), type="character", default=NULL, help="Condition level for ident/group 1. If unset, first sorted level is used."),
  make_option(c("--condition-b"), type="character", default=NULL, help="Condition level for ident/group 2. If unset, second sorted level is used."),
  make_option(c("--celltype-a"), type="character", default=NULL, help="Cell type level for ident/group 1. If unset, first sorted level is used."),
  make_option(c("--celltype-b"), type="character", default=NULL, help="Cell type level for ident/group 2. If unset, second sorted level is used."),
  make_option(c("--target-celltype"), type="character", default=NULL, help="Cell type used when --de-scope target_celltype is selected."),
  make_option(c("--sample-col"), type="character", default="orig.ident", help="Metadata sample/replicate column for pseudobulk and sample-aware reporting [default %default]"),
  make_option(c("--condition-only-mode"), type="character", default="setdiff", help="Condition-only mode: setdiff. Future: interaction/pseudobulk_covariate [default %default]"),
  make_option(c("--pseudobulk-min-replicates"), type="integer", default=2, help="Minimum pseudobulk samples per group [default %default]"),
  make_option(c("--pseudobulk-method"), type="character", default="limma_voom", help="Pseudobulk DE method: limma_voom or ttest [default %default]"),
  make_option(c("--enrich-backend"), type="character", default="gprof", help="gprof or enrichr [default %default]"),
  make_option(c("--gprof-sources"), type="character", default="GO:BP,GO:MF,GO:CC,KEGG,REAC", help="Comma-separated g:Profiler sources"),
  make_option(c("--enrichr-db"), type="character", default="GO_Biological_Process_2023", help="Enrichr library name"),
  make_option(c("--top-n-features"), type="integer", default=10, help="Top RF features for classifier/heatmap [default %default]"),
  make_option(c("--max-pcs"), type="integer", default=30, help="Maximum PCs/dims for PCA/UMAP [default %default]"),
  make_option(c("--skip-pca-umap"), action="store_true", default=FALSE, help="Skip PCA/UMAP outputs"),
  make_option(c("--skip-enrichment"), action="store_true", default=FALSE, help="Skip enrichment outputs"),
  make_option(c("--skip-classifier"), action="store_true", default=FALSE, help="Skip RF feature selection/classifier"),
  make_option(c("--classifier-validation"), type="character", default="sample_cv", help="Classifier validation: sample_cv, cell_split, or none [default %default]"),
  make_option(c("--classifier-model"), type="character", default="auto", help="Classifier model: auto, rf, or logistic [default %default]"),
  make_option(c("--classifier-selector"), type="character", default="auto", help="Feature selector inside training data: auto, rf, logistic, or none [default %default]"),
  make_option(c("--classifier-feature-source"), type="character", default="condition_only", help="Candidate feature source inside each training split: condition_only or condition_de [default %default]"),
  make_option(c("--classifier-min-samples-per-class"), type="integer", default=2, help="Minimum biological samples per class for sample-level CV [default %default]"),
  make_option(c("--skip-power"), action="store_true", default=FALSE, help="Skip power analysis"),
  make_option(c("--zip"), action="store_true", default=FALSE, help="Zip output directory at end")
)

opt <- parse_args(OptionParser(option_list=option_list))

write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, pretty=TRUE, auto_unbox=TRUE, null="null")
}

fail <- function(msg, code=1) {
  cat(sprintf("ERROR: %s\n", msg), file=stderr())
  quit(status=code)
}

if (is.null(opt$counts) || !file.exists(opt$counts)) fail("--counts is required and must exist")
if (is.null(opt$metadata) || !file.exists(opt$metadata)) fail("--metadata is required and must exist")
if (is.null(opt$outdir) || !nzchar(opt$outdir)) fail("--outdir must be a non-empty path")

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
tables_dir <- file.path(opt$outdir, "tables"); dir.create(tables_dir, recursive=TRUE, showWarnings=FALSE)
plots_dir <- file.path(opt$outdir, "plots"); dir.create(plots_dir, recursive=TRUE, showWarnings=FALSE)
objects_dir <- file.path(opt$outdir, "objects"); dir.create(objects_dir, recursive=TRUE, showWarnings=FALSE)

run_status <- list(
  module="scRNA",
  version=VERSION,
  started=as.character(Sys.time()),
  finished=NULL,
  status="running",
  outdir=normalizePath(opt$outdir, mustWork=FALSE),
  parameters=list(
    species=opt$species,
    id_type=opt$`id-type`,
    condition_col=opt$`condition-col`,
    celltype_col=opt$`celltype-col`,
    annotation_col=opt$`annotation-col`,
    de_fdr=opt$`de-fdr`,
    logfc_threshold=opt$`logfc-threshold`,
    min_pct=opt$`min-pct`,
    test_use=opt$`test-use`,
    de_mode=opt$`de-mode`,
    de_scope=opt$`de-scope`,
    condition_a=opt$`condition-a`,
    condition_b=opt$`condition-b`,
    celltype_a=opt$`celltype-a`,
    celltype_b=opt$`celltype-b`,
    target_celltype=opt$`target-celltype`,
    sample_col=opt$`sample-col`,
    condition_only_mode=opt$`condition-only-mode`,
    pseudobulk_min_replicates=opt$`pseudobulk-min-replicates`,
    pseudobulk_method=opt$`pseudobulk-method`,
    enrich_backend=opt$`enrich-backend`,
    gprof_sources=opt$`gprof-sources`,
    enrichr_db=opt$`enrichr-db`,
    top_n_features=opt$`top-n-features`,
    classifier_validation=opt$`classifier-validation`,
    classifier_model=opt$`classifier-model`,
    classifier_selector=opt$`classifier-selector`,
    classifier_feature_source=opt$`classifier-feature-source`,
    classifier_min_samples_per_class=opt$`classifier-min-samples-per-class`,
    max_pcs=opt$`max-pcs`
  ),
  steps=list(),
  warnings=list(),
  errors=list(),
  files=list(tables=list(), plots=list(), objects=list())
)

write_status <- function() {
  write_json_file(run_status, file.path(opt$outdir, "run_status.json"))
}

add_warning <- function(msg, step=NULL) {
  run_status$warnings[[length(run_status$warnings) + 1]] <<- list(
    time=as.character(Sys.time()),
    step=step,
    message=msg
  )
  write_status()
}

mark_skipped <- function(name, reason) {
  now <- as.character(Sys.time())
  run_status$steps[[name]] <<- list(
    status="skipped",
    started=now,
    finished=now,
    seconds=0,
    reason=reason
  )
  add_warning(paste0("Skipped ", name, ": ", reason), step=name)
}

record_step <- function(name, expr) {
  started <- Sys.time()
  run_status$steps[[name]] <<- list(status="running", started=as.character(started))
  write_status()
  tryCatch({
    value <- force(expr)
    finished <- Sys.time()
    run_status$steps[[name]] <<- list(
      status="complete",
      started=as.character(started),
      finished=as.character(finished),
      seconds=round(as.numeric(difftime(finished, started, units="secs")), 3)
    )
    write_status()
    value
  }, error=function(e) {
    finished <- Sys.time()
    run_status$steps[[name]] <<- list(
      status="failed",
      started=as.character(started),
      finished=as.character(finished),
      seconds=round(as.numeric(difftime(finished, started, units="secs")), 3),
      error=conditionMessage(e)
    )
    run_status$status <<- "failed"
    run_status$finished <<- as.character(finished)
    run_status$errors[[length(run_status$errors) + 1]] <<- list(step=name, message=conditionMessage(e))
    write_status()
    stop(e)
  })
}

write_status()

write_csv <- function(x, path, row_names=FALSE) {
  if (is.null(x)) return(invisible(NULL))
  write.csv(x, path, row.names=row_names)
  invisible(path)
}

species_map <- function(x) {
  x0 <- x
  x <- tolower(trimws(x))
  switch(x,
    "human"="hsapiens", "hsapiens"="hsapiens",
    "mouse"="mmusculus", "mmusculus"="mmusculus",
    "fly"="dmelanogaster", "drosophila"="dmelanogaster", "dmelanogaster"="dmelanogaster",
    "zebrafish"="drerio", "drerio"="drerio",
    stop(paste0("Unsupported species: ", x0))
  )
}

map_gprof_sources <- function(srcs) {
  srcs <- trimws(unlist(strsplit(srcs, ",")))
  srcs <- gsub("^MSigDB H$", "MSigDBH", srcs)
  srcs <- gsub("^MSigDB C2$", "MSigDBC2", srcs)
  srcs[nzchar(srcs)]
}

strip_bg <- function(x) unique(x[!grepl("^BG\\d+$", x)])

find_gene_list_col <- function(df) {
  candidates <- c("Genes", "genes", "GENES", "intersection", "intersections")
  for (nm in candidates) if (nm %in% colnames(df)) return(as.character(df[[nm]]))
  idx <- grep("gene", colnames(df), ignore.case=TRUE)
  if (length(idx) > 0) return(as.character(df[[idx[1]]]))
  rep(NA_character_, nrow(df))
}

fmt_enrich_for_table <- function(res) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  if (all(c("term_name", "term_size", "intersection_size", "p_value") %in% colnames(res))) {
    adj_col <- intersect(colnames(res), c("p_adjusted", "adjusted_p_value", "p_adj"))
    adj <- if (length(adj_col)) res[[adj_col[1]]] else p.adjust(res$p_value, method="fdr")
    df <- data.frame(
      Term=paste0(res$term_name, " (", res$source, ")"),
      Overlap=paste0(res$intersection_size, "/", res$term_size),
      P.value=as.numeric(res$p_value),
      Adjusted.P.value=as.numeric(adj),
      Genes=find_gene_list_col(res),
      stringsAsFactors=FALSE,
      check.names=FALSE
    )
    df <- df[order(df$Adjusted.P.value, df$P.value), , drop=FALSE]
    rownames(df) <- NULL
    return(df)
  }
  if (all(c("Term", "Overlap", "P.value", "Adjusted.P.value") %in% colnames(res))) {
    term <- if ("source" %in% colnames(res)) paste0(res$Term, " (", res$source, ")") else paste0(res$Term, " (ENRICHR)")
    df <- data.frame(
      Term=term,
      Overlap=res$Overlap,
      P.value=as.numeric(res$P.value),
      Adjusted.P.value=as.numeric(res$Adjusted.P.value),
      Genes=find_gene_list_col(res),
      stringsAsFactors=FALSE,
      check.names=FALSE
    )
    df <- df[order(df$Adjusted.P.value, df$P.value), , drop=FALSE]
    rownames(df) <- NULL
    return(df)
  }
  NULL
}

make_edges <- function(enrich_tbl) {
  if (is.null(enrich_tbl) || !"Genes" %in% colnames(enrich_tbl)) return(NULL)
  tmp <- enrich_tbl[!is.na(enrich_tbl$Genes) & nchar(enrich_tbl$Genes) > 0, c("Term", "Genes"), drop=FALSE]
  if (nrow(tmp) == 0) return(NULL)
  edges_list <- lapply(seq_len(nrow(tmp)), function(i) {
    genes <- unlist(strsplit(as.character(tmp$Genes[i]), "[,;]"))
    data.frame(Term=tmp$Term[i], Gene=trimws(genes), stringsAsFactors=FALSE)
  })
  edges <- do.call(rbind, edges_list)
  edges <- edges[nchar(edges$Gene) > 0, , drop=FALSE]
  rownames(edges) <- NULL
  edges
}

safe_gprof <- function(genes, org, sources, retries=3) {
  genes <- unique(genes)
  waits <- c(0.6, 1.2, 2.4)
  for (i in seq_len(retries)) {
    out <- tryCatch(
      gprofiler2::gost(query=genes, organism=org, sources=sources, correction_method="fdr", evcodes=TRUE),
      error=function(e) {
        add_warning(paste("g:Profiler attempt failed:", conditionMessage(e)), step="enrichment")
        NULL
      }
    )
    if (!is.null(out$result) && nrow(out$result) > 0) return(out$result)
    if (i < retries) Sys.sleep(waits[i])
  }
  NULL
}

do_enrichment <- function(genes, species, backend, sources, db) {
  genes <- strip_bg(genes)
  if (length(genes) < 3) return(NULL)
  backend <- tolower(trimws(backend))
  if (backend == "gprof") {
    return(fmt_enrich_for_table(safe_gprof(genes, species, sources)))
  }
  if (backend == "enrichr") {
    er <- tryCatch(enrichR::enrichr(unique(genes), db), error=function(e) {
      add_warning(paste("Enrichr failed:", conditionMessage(e)), step="enrichment")
      NULL
    })
    if (is.null(er) || is.null(er[[db]]) || nrow(er[[db]]) == 0) return(NULL)
    return(fmt_enrich_for_table(er[[db]]))
  }
  add_warning(paste("Unknown enrichment backend:", backend), step="enrichment")
  NULL
}

save_plotly <- function(p, path) {
  ok <- tryCatch({
    htmlwidgets::saveWidget(plotly::ggplotly(p), path, selfcontained=TRUE)
    TRUE
  }, error=function(e) {
    add_warning(paste("Could not save plot", basename(path), ":", conditionMessage(e)), step="plotting")
    FALSE
  })
  invisible(ok)
}



annotate_de_table <- function(df, comparison, de_mode, de_method,
                              pseudobulk_samples_1=NA, pseudobulk_samples_2=NA) {
  if (is.null(df)) {
    df <- data.frame()
  }

  n <- nrow(df)

  if (!"gene" %in% colnames(df)) {
    df$gene <- rownames(df)
  }

  df$comparison <- rep(as.character(comparison), n)
  df$de_mode <- rep(as.character(de_mode), n)
  df$de_method <- rep(as.character(de_method), n)

  if (!is.na(pseudobulk_samples_1) || !is.na(pseudobulk_samples_2)) {
    df$pseudobulk_samples_1 <- rep(pseudobulk_samples_1, n)
    df$pseudobulk_samples_2 <- rep(pseudobulk_samples_2, n)
  }

  df
}


get_assay_matrix <- function(seu, assay="RNA", layer="data") {
  tryCatch(
    {
      Seurat::GetAssayData(seu, assay=assay, layer=layer)
    },
    error=function(e1) {
      Seurat::GetAssayData(seu, assay=assay, slot=layer)
    }
  )
}

validate_inputs <- function(counts, meta, condition_col, celltype_col) {
  if (anyDuplicated(rownames(counts))) stop("Counts matrix has duplicated gene IDs. Make row names unique before running.")
  if (anyDuplicated(colnames(counts))) stop("Counts matrix has duplicated cell IDs in column names.")
  if (anyDuplicated(rownames(meta))) stop("Metadata has duplicated cell IDs in row names.")
  if (nrow(counts) < 2) stop("Counts matrix must contain at least 2 genes.")
  if (ncol(counts) < 4) stop("Counts matrix must contain at least 4 cells for this pipeline.")

  common_cells <- intersect(colnames(counts), rownames(meta))
  if (length(common_cells) < 4) stop("Counts columns and metadata rownames must share at least 4 cell IDs.")

  numeric_ok <- vapply(counts, is.numeric, logical(1))
  if (!all(numeric_ok)) stop("Counts matrix contains non-numeric columns. Check CSV formatting.")

  missing_cols <- setdiff(c(condition_col, celltype_col), colnames(meta))
  if (length(missing_cols) > 0) {
    add_warning(paste("Missing metadata columns:", paste(missing_cols, collapse=", ")), step="validate_inputs")
  }

  if (condition_col %in% colnames(meta)) {
    cond_n <- length(unique(na.omit(trimws(as.character(meta[[condition_col]])))))
    if (cond_n < 2) add_warning(paste("Condition column has fewer than 2 groups:", condition_col), step="validate_inputs")
  }

  if (celltype_col %in% colnames(meta)) {
    ct_n <- length(unique(na.omit(trimws(as.character(meta[[celltype_col]])))))
    if (ct_n < 2) add_warning(paste("Cell type column has fewer than 2 groups:", celltype_col), step="validate_inputs")
  }

  invisible(common_cells)
}

table_info <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read.csv(path, check.names=FALSE), error=function(e) NULL)
  list(
    name=tools::file_path_sans_ext(basename(path)),
    path=file.path("tables", basename(path)),
    label=label,
    rows=if (is.null(df)) NA else nrow(df),
    columns=if (is.null(df)) NA else ncol(df),
    bytes=file.info(path)$size
  )
}

plot_info <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  list(
    name=tools::file_path_sans_ext(basename(path)),
    path=file.path("plots", basename(path)),
    label=label,
    type="html",
    bytes=file.info(path)$size
  )
}

object_info <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  list(
    name=tools::file_path_sans_ext(basename(path)),
    path=file.path("objects", basename(path)),
    label=label,
    type=tools::file_ext(path),
    bytes=file.info(path)$size
  )
}


validate_de_options <- function() {
  de_mode <- tolower(trimws(opt$`de-mode`))
  if (!(de_mode %in% c("cell", "pseudobulk"))) stop("--de-mode must be one of: cell, pseudobulk")
  pb_method <- tolower(trimws(opt$`pseudobulk-method`))
  if (!(pb_method %in% c("limma_voom", "ttest"))) stop("--pseudobulk-method must be one of: limma_voom, ttest")
  de_scope <- tolower(trimws(opt$`de-scope`))
  if (!(de_scope %in% c("global", "top_celltype", "target_celltype"))) stop("--de-scope must be one of: global, top_celltype, target_celltype")
  if (de_scope == "target_celltype" && (is.null(opt$`target-celltype`) || !nzchar(opt$`target-celltype`))) {
    stop("--target-celltype is required when --de-scope target_celltype is selected")
  }
  if (tolower(trimws(opt$`condition-only-mode`)) != "setdiff") {
    stop("Only --condition-only-mode setdiff is currently implemented")
  }
  allowed_tests <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
  if (!(opt$`test-use` %in% allowed_tests)) {
    add_warning(paste0("Uncommon Seurat test.use='", opt$`test-use`, "'. Seurat may error if unsupported."), step="validate_inputs")
  }
  invisible(TRUE)
}

pick_two_levels <- function(values, a=NULL, b=NULL, label="group", step="differential_expression") {
  lv <- sort(unique(na.omit(trimws(as.character(values)))))
  if (!is.null(a) && nzchar(a) && !is.null(b) && nzchar(b)) {
    if (!(a %in% lv)) stop(paste0(label, " level not found: ", a))
    if (!(b %in% lv)) stop(paste0(label, " level not found: ", b))
    if (identical(a, b)) stop(paste0(label, " contrast levels must differ"))
    return(c(a, b))
  }
  if (length(lv) < 2) stop(paste0("Need at least two levels for ", label))
  add_warning(paste0("No explicit ", label, " contrast supplied; using sorted levels: ", lv[1], " vs ", lv[2], "."), step=step)
  c(lv[1], lv[2])
}

select_condition_cells <- function(seu, de_scope, celltype_col, target_celltype=NULL) {
  de_scope <- tolower(trimws(de_scope))
  if (de_scope == "global") return(colnames(seu))
  if (!(celltype_col %in% colnames(seu@meta.data))) stop("Cell type column is required for non-global DE scope")
  ct <- trimws(as.character(seu[[celltype_col]][,1]))
  if (de_scope == "top_celltype") {
    top_ct <- names(sort(table(ct), decreasing=TRUE))[1]
    add_warning(paste0("Condition DE scoped to most abundant cell type: ", top_ct), step="differential_expression")
    return(colnames(seu)[ct == top_ct])
  }
  if (de_scope == "target_celltype") {
    if (!(target_celltype %in% unique(ct))) stop(paste0("--target-celltype not found in metadata: ", target_celltype))
    return(colnames(seu)[ct == target_celltype])
  }
  colnames(seu)
}

run_cell_findmarkers <- function(seu, group_col, group_a=NULL, group_b=NULL, cells=NULL, label="DE") {
  if (!(group_col %in% colnames(seu@meta.data))) stop(paste0("Missing metadata column: ", group_col))

  obj <- seu
  if (!is.null(cells)) {
    cells <- intersect(cells, colnames(seu))
    if (length(cells) < 4) stop(paste0(label, ": fewer than 4 cells after subsetting"))
    obj <- subset(seu, cells=cells)
  }

  obj[[group_col]][,1] <- factor(trimws(as.character(obj[[group_col]][,1])))
  contrast <- pick_two_levels(obj[[group_col]][,1], group_a, group_b, label=label)

  tab <- table(obj[[group_col]][,1])
  if (any(tab[contrast] < 2)) {
    add_warning(
      paste0(label, " has a group with fewer than 2 cells; p-values may be unstable."),
      step="differential_expression"
    )
  }

  Idents(obj) <- group_col

  de <- tryCatch(
    {
      FindMarkers(
        obj,
        ident.1=contrast[1],
        ident.2=contrast[2],
        test.use=opt$`test-use`,
        min.pct=opt$`min-pct`,
        logfc.threshold=opt$`logfc-threshold`,
        verbose=FALSE
      )
    },
    error=function(e) {
      stop(e)
    }
  )

  if (is.null(de) || nrow(de) == 0) {
    add_warning(
      paste0(label, " returned 0 DE rows for contrast: ", contrast[1], " vs ", contrast[2]),
      step="differential_expression"
    )

    de <- data.frame(
      p_val=numeric(0),
      avg_log2FC=numeric(0),
      pct.1=numeric(0),
      pct.2=numeric(0),
      p_val_adj=numeric(0),
      stringsAsFactors=FALSE
    )
  }

  # Critical: annotation must work even when de has 0 rows.
  de <- annotate_de_table(
    de,
    comparison=paste0(label, ": ", contrast[1], " vs ", contrast[2]),
    de_mode="cell",
    de_method=paste0("Seurat FindMarkers test.use=", opt$`test-use`)
  )

  de
}

run_pseudobulk_ttest <- function(pb_counts, pb_meta, contrast, label) {
  lib_size <- colSums(pb_counts)
  keep_lib <- is.finite(lib_size) & lib_size > 0
  pb_counts <- pb_counts[, keep_lib, drop=FALSE]
  pb_meta <- pb_meta[colnames(pb_counts), , drop=FALSE]
  lib_size <- lib_size[keep_lib]
  cpm <- t(t(pb_counts) / lib_size * 1e6)
  logcpm <- log2(cpm + 1)

  g1 <- which(pb_meta$.group == contrast[1])
  g2 <- which(pb_meta$.group == contrast[2])
  pvals <- apply(logcpm, 1, function(x) {
    x1 <- x[g1]; x2 <- x[g2]
    if (length(unique(x1)) < 2 && length(unique(x2)) < 2) return(1)
    tryCatch(stats::t.test(x1, x2)$p.value, error=function(e) NA_real_)
  })
  avg1 <- rowMeans(logcpm[, g1, drop=FALSE])
  avg2 <- rowMeans(logcpm[, g2, drop=FALSE])
  pct1 <- rowMeans(pb_counts[, g1, drop=FALSE] > 0)
  pct2 <- rowMeans(pb_counts[, g2, drop=FALSE] > 0)
  df <- data.frame(
    p_val=pvals,
    avg_log2FC=avg1 - avg2,
    pct.1=pct1,
    pct.2=pct2,
    p_val_adj=p.adjust(pvals, method="BH"),
    gene=rownames(logcpm),
    comparison=paste0(label, ": ", contrast[1], " vs ", contrast[2]),
    de_mode="pseudobulk",
    de_method="pseudobulk logCPM Welch t-test",
    pseudobulk_method="ttest",
    pseudobulk_samples_1=length(g1),
    pseudobulk_samples_2=length(g2),
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
  df <- df[order(df$p_val_adj, df$p_val), , drop=FALSE]
  rownames(df) <- df$gene
  df
}

run_pseudobulk_limma_voom <- function(pb_counts, pb_meta, contrast, label) {
  if (!requireNamespace("edgeR", quietly=TRUE)) {
    stop("--pseudobulk-method limma_voom requires the Bioconductor package edgeR. Install with BiocManager::install('edgeR').")
  }
  if (!requireNamespace("limma", quietly=TRUE)) {
    stop("--pseudobulk-method limma_voom requires the Bioconductor package limma. Install with BiocManager::install('limma').")
  }

  lib_size <- colSums(pb_counts)
  keep_lib <- is.finite(lib_size) & lib_size > 0
  pb_counts <- pb_counts[, keep_lib, drop=FALSE]
  pb_meta <- pb_meta[colnames(pb_counts), , drop=FALSE]

  group <- factor(pb_meta$.group, levels=contrast)
  g1 <- which(group == contrast[1])
  g2 <- which(group == contrast[2])

  dge <- edgeR::DGEList(counts=pb_counts, group=group)
  keep_genes <- edgeR::filterByExpr(dge, group=group)
  if (sum(keep_genes) < 2) {
    stop(paste0(label, ": fewer than 2 genes remained after edgeR::filterByExpr."))
  }
  dge <- dge[keep_genes, , keep.lib.sizes=FALSE]
  dge <- edgeR::calcNormFactors(dge)

  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- make.names(levels(group))
  contrast_expr <- paste0(make.names(contrast[1]), "-", make.names(contrast[2]))

  residual_df <- ncol(dge$counts) - qr(design)$rank
  if (!is.finite(residual_df) || residual_df <= 0) {
    stop(paste0(
      label, " limma/voom has no residual degrees of freedom: samples=",
      ncol(dge$counts), ", design_rank=", qr(design)$rank,
      ". Add biological replicates per group or use --de-mode cell for exploratory cell-level analysis."
    ))
  }

  v <- limma::voom(dge, design=design, plot=FALSE)
  cm <- limma::makeContrasts(contrasts=contrast_expr, levels=design)
  fit <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit, cm)
  fit <- limma::eBayes(fit, trend=TRUE)
  tt <- limma::topTable(fit, coef=1, number=Inf, sort.by="P")

  if (is.null(tt) || nrow(tt) == 0) {
    add_warning(paste0(label, " limma/voom returned 0 DE rows."), step="differential_expression")
    tt <- data.frame(logFC=numeric(0), AveExpr=numeric(0), t=numeric(0), P.Value=numeric(0), adj.P.Val=numeric(0), B=numeric(0))
  }

  genes <- rownames(tt)
  pct1 <- if (length(genes)) rowMeans(pb_counts[genes, g1, drop=FALSE] > 0) else numeric(0)
  pct2 <- if (length(genes)) rowMeans(pb_counts[genes, g2, drop=FALSE] > 0) else numeric(0)

  df <- data.frame(
    p_val=as.numeric(tt$P.Value),
    avg_log2FC=as.numeric(tt$logFC),
    pct.1=as.numeric(pct1),
    pct.2=as.numeric(pct2),
    p_val_adj=as.numeric(tt$adj.P.Val),
    gene=genes,
    AveExpr=if ("AveExpr" %in% colnames(tt)) as.numeric(tt$AveExpr) else NA_real_,
    t=if ("t" %in% colnames(tt)) as.numeric(tt$t) else NA_real_,
    B=if ("B" %in% colnames(tt)) as.numeric(tt$B) else NA_real_,
    comparison=paste0(label, ": ", contrast[1], " vs ", contrast[2]),
    de_mode="pseudobulk",
    de_method="edgeR::DGEList + edgeR::filterByExpr + edgeR::calcNormFactors + limma::voom + limma::eBayes(trend=TRUE)",
    pseudobulk_method="limma_voom",
    pseudobulk_samples_1=length(g1),
    pseudobulk_samples_2=length(g2),
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
  df <- df[order(df$p_val_adj, df$p_val), , drop=FALSE]
  rownames(df) <- df$gene
  df
}

run_pseudobulk_de <- function(seu, group_col, group_a=NULL, group_b=NULL, cells=NULL, label="Pseudobulk DE") {
  sample_col <- opt$`sample-col`
  method <- tolower(trimws(opt$`pseudobulk-method`))
  if (!(method %in% c("limma_voom", "ttest"))) stop("--pseudobulk-method must be one of: limma_voom, ttest")
  if (!(sample_col %in% colnames(seu@meta.data))) stop(paste0("Pseudobulk requires sample column: ", sample_col))
  if (!(group_col %in% colnames(seu@meta.data))) stop(paste0("Missing metadata column for pseudobulk grouping: ", group_col))

  use_cells <- if (is.null(cells)) colnames(seu) else intersect(cells, colnames(seu))
  if (length(use_cells) < 4) stop(paste0(label, ": fewer than 4 cells after subsetting"))
  md <- seu@meta.data[use_cells, , drop=FALSE]
  md$.sample <- trimws(as.character(md[[sample_col]]))
  md$.group <- trimws(as.character(md[[group_col]]))
  keep <- nzchar(md$.sample) & nzchar(md$.group) & !is.na(md$.sample) & !is.na(md$.group)
  md <- md[keep, , drop=FALSE]
  use_cells <- rownames(md)
  if (length(use_cells) < 4) stop(paste0(label, ": fewer than 4 usable cells after removing missing sample/group values"))

  contrast <- pick_two_levels(md$.group, group_a, group_b, label=label)
  md <- md[md$.group %in% contrast, , drop=FALSE]
  use_cells <- rownames(md)
  if (length(use_cells) < 4) stop(paste0(label, ": fewer than 4 cells in selected contrast"))

  md$.pb_id <- paste(md$.sample, md$.group, sep="__")
  raw_counts <- as.matrix(get_assay_matrix(seu, assay="RNA", layer="counts")[, use_cells, drop=FALSE])
  pb_counts_t <- rowsum(t(raw_counts), group=md$.pb_id, reorder=FALSE)
  pb_counts <- t(pb_counts_t)
  pb_meta <- unique(md[, c(".pb_id", ".sample", ".group"), drop=FALSE])
  rownames(pb_meta) <- pb_meta$.pb_id
  pb_meta <- pb_meta[colnames(pb_counts), , drop=FALSE]

  rep_counts <- table(pb_meta$.group)
  min_reps <- opt$`pseudobulk-min-replicates`

  # limma/voom needs biological replication to estimate residual variance.
  # Do not allow --pseudobulk-min-replicates 1 to sneak into limma/voom.
  effective_min_reps <- min_reps
  if (method == "limma_voom") {
    effective_min_reps <- max(2, min_reps)
  }

  if (any(rep_counts[contrast] < effective_min_reps)) {
    stop(paste0(
      label, " using --pseudobulk-method ", method,
      " requires at least ", effective_min_reps,
      " pseudobulk samples per group. Observed: ",
      paste(names(rep_counts), rep_counts, sep="=", collapse=", "),
      ". For limma/voom, this is required because the model needs residual degrees of freedom to estimate variance."
    ))
  }

  if (method == "limma_voom") {
    return(run_pseudobulk_limma_voom(pb_counts, pb_meta, contrast, label))
  }
  run_pseudobulk_ttest(pb_counts, pb_meta, contrast, label)
}


statistical_notes <- function() {
  notes <- c(
    paste0("DE mode: ", opt$`de-mode`, "."),
    paste0("Condition-only mode: ", opt$`condition-only-mode`, "; this is a set difference, not a formal condition-by-cell-type interaction test."),
    "Adjusted p-values are corrected within each output table, not across the entire workflow.",
    "PCA/UMAP are exploratory visualization outputs and should not be interpreted as formal statistical tests.",
    "Enrichment depends on the submitted DE gene list and default service background unless a future custom background option is added.",
    "Classifier feature selection is exploratory unless validation is sample-aware and feature selection is performed inside each training split.",
    "Cell-level classifier splits can be optimistic if cells from the same biological sample are split across train/test; sample-level validation is preferred."
  )
  if (tolower(trimws(opt$`de-mode`)) == "cell") {
    notes <- c(notes, "Cell-level DE treats cells as observations; for replicate-level inference, use --de-mode pseudobulk with a valid --sample-col.")
  } else {
    notes <- c(notes, paste0("Pseudobulk DE aggregates raw counts by sample/replicate and group. Current pseudobulk method: ", opt$`pseudobulk-method`, "."))
    if (tolower(trimws(opt$`pseudobulk-method`)) == "limma_voom") {
      notes <- c(notes, "limma/voom pseudobulk uses edgeR filtering/TMM normalization, voom precision weights, and limma empirical Bayes moderation for replicate-level inference.")
    } else {
      notes <- c(notes, "The pseudobulk t-test method is a lightweight fallback on logCPM values; limma_voom is preferred for serious pseudobulk DE.")
    }
  }
  notes
}

validate_classifier_options <- function() {
  validation <- tolower(trimws(opt$`classifier-validation`))
  model <- tolower(trimws(opt$`classifier-model`))
  selector <- tolower(trimws(opt$`classifier-selector`))
  source <- tolower(trimws(opt$`classifier-feature-source`))

  if (!(validation %in% c("sample_cv", "cell_split", "none"))) {
    stop("--classifier-validation must be one of: sample_cv, cell_split, none")
  }
  if (!(model %in% c("auto", "rf", "logistic"))) {
    stop("--classifier-model must be one of: auto, rf, logistic")
  }
  if (!(selector %in% c("auto", "rf", "logistic", "none"))) {
    stop("--classifier-selector must be one of: auto, rf, logistic, none")
  }
  if (!(source %in% c("condition_only", "condition_de"))) {
    stop("--classifier-feature-source must be one of: condition_only, condition_de")
  }
  invisible(TRUE)
}

safe_binary_levels <- function(y) {
  y <- droplevels(as.factor(y))
  if (length(levels(y)) != 2) return(NULL)
  levels(y)
}

make_safe_feature_frame <- function(mat) {
  df <- as.data.frame(mat, check.names=FALSE)
  original <- colnames(df)
  safe <- make.names(original, unique=TRUE)
  colnames(df) <- safe
  attr(df, "feature_map") <- data.frame(Safe=safe, Gene=original, stringsAsFactors=FALSE)
  df
}

rank_features_training_only <- function(x_train, y_train, method="auto", top_n=10) {
  y_train <- droplevels(as.factor(y_train))
  fmap <- attr(x_train, "feature_map")
  if (is.null(fmap)) {
    fmap <- data.frame(Safe=colnames(x_train), Gene=colnames(x_train), stringsAsFactors=FALSE)
  }

  keep <- vapply(x_train, function(z) is.numeric(z) && stats::sd(z, na.rm=TRUE) > 0, logical(1))
  x_train <- x_train[, keep, drop=FALSE]
  fmap <- fmap[fmap$Safe %in% colnames(x_train), , drop=FALSE]

  if (ncol(x_train) == 0) return(fmap[0, , drop=FALSE])

  method <- tolower(trimws(method))
  if (method == "auto") {
    method <- if (length(levels(y_train)) == 2 && nrow(x_train) < 80) "logistic" else "rf"
  }

  if (method == "none") {
    out <- fmap
    out$Score <- NA_real_
    out$Selector <- "none"
    return(head(out, top_n))
  }

  if (method == "logistic" && length(levels(y_train)) == 2) {
    y01 <- as.numeric(y_train == levels(y_train)[2])
    scores <- vapply(colnames(x_train), function(nm) {
      dat <- data.frame(y=y01, x=x_train[[nm]])
      fit <- tryCatch(stats::glm(y ~ x, data=dat, family=stats::binomial()), error=function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      co <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
      if (is.null(co) || nrow(co) < 2) return(NA_real_)
      z <- suppressWarnings(abs(co[2, "z value"]))
      if (!is.finite(z)) return(NA_real_)
      z
    }, numeric(1))
    out <- merge(fmap, data.frame(Safe=names(scores), Score=as.numeric(scores), stringsAsFactors=FALSE), by="Safe", all.x=FALSE)
    out <- out[order(out$Score, decreasing=TRUE, na.last=NA), , drop=FALSE]
    out$Selector <- "training-only univariate logistic"
    rownames(out) <- NULL
    return(head(out, top_n))
  }

  rf <- randomForest::randomForest(x=x_train, y=y_train, importance=TRUE, ntree=500)
  imp <- as.data.frame(randomForest::importance(rf))
  imp$Safe <- rownames(imp)
  score_col <- intersect(colnames(imp), c("MeanDecreaseGini", "MeanDecreaseAccuracy"))
  if (length(score_col) == 0) {
    numeric_cols <- colnames(imp)[vapply(imp, is.numeric, logical(1))]
    score_col <- numeric_cols[1]
  }
  imp$Score <- as.numeric(imp[[score_col[1]]])
  out <- merge(fmap, imp[, c("Safe", "Score"), drop=FALSE], by="Safe", all.x=FALSE)
  out <- out[order(out$Score, decreasing=TRUE, na.last=NA), , drop=FALSE]
  out$Selector <- "training-only random forest"
  rownames(out) <- NULL
  head(out, top_n)
}

train_predict_classifier <- function(x_train, y_train, x_test, model="auto") {
  y_train <- droplevels(as.factor(y_train))
  model <- tolower(trimws(model))
  if (model == "auto") {
    model <- if (length(levels(y_train)) == 2 && nrow(x_train) >= (ncol(x_train) + 5)) "logistic" else "rf"
  }

  if (model == "logistic" && length(levels(y_train)) == 2) {
    dat <- data.frame(y=as.numeric(y_train == levels(y_train)[2]), x_train, check.names=FALSE)
    fit <- tryCatch(stats::glm(y ~ ., data=dat, family=stats::binomial()), error=function(e) NULL)
    if (!is.null(fit)) {
      prob <- tryCatch(as.numeric(stats::predict(fit, newdata=x_test, type="response")), error=function(e) NULL)
      if (!is.null(prob) && length(prob) == nrow(x_test) && all(is.finite(prob))) {
        pred <- factor(ifelse(prob >= 0.5, levels(y_train)[2], levels(y_train)[1]), levels=levels(y_train))
        return(list(pred=pred, prob=prob, model="logistic"))
      }
    }
  }

  rf <- randomForest::randomForest(x=x_train, y=y_train, ntree=500)
  pred <- predict(rf, newdata=x_test, type="response")
  probs <- tryCatch({
    pp <- predict(rf, newdata=x_test, type="prob")
    if (ncol(pp) >= 2) as.numeric(pp[, 2]) else rep(NA_real_, nrow(x_test))
  }, error=function(e) rep(NA_real_, nrow(x_test)))
  list(pred=pred, prob=probs, model="rf")
}

select_training_de_candidates <- function(seu, train_cells, condition_col, celltype_col, feature_source="condition_only") {
  train_cells <- intersect(train_cells, colnames(seu))
  condition_train <- run_cell_findmarkers(
    seu,
    group_col=condition_col,
    group_a=opt$`condition-a`,
    group_b=opt$`condition-b`,
    cells=train_cells,
    label="classifier training-fold condition DE"
  )

  cond_genes <- rownames(subset(
    condition_train,
    p_val_adj < opt$`de-fdr` & abs(avg_log2FC) >= opt$`logfc-threshold`
  ))

  feature_source <- tolower(trimws(feature_source))
  if (feature_source == "condition_de" || !(celltype_col %in% colnames(seu@meta.data))) {
    return(unique(cond_genes))
  }

  celltype_train <- tryCatch(
    run_cell_findmarkers(
      seu,
      group_col=celltype_col,
      group_a=opt$`celltype-a`,
      group_b=opt$`celltype-b`,
      cells=train_cells,
      label="classifier training-fold cell type DE"
    ),
    error=function(e) {
      add_warning(paste("Training-fold cell type DE failed; using condition DE candidates:", conditionMessage(e)), step="classifier")
      NULL
    }
  )

  if (is.null(celltype_train)) return(unique(cond_genes))

  ct_genes <- rownames(subset(
    celltype_train,
    p_val_adj < opt$`de-fdr` & abs(avg_log2FC) >= opt$`logfc-threshold`
  ))

  unique(setdiff(cond_genes, ct_genes))
}

sample_cv_folds <- function(meta, sample_col, condition_col, min_samples_per_class=2) {
  if (!(sample_col %in% colnames(meta))) stop(paste0("Missing sample column for sample-level classifier validation: ", sample_col))
  md <- meta[, c(sample_col, condition_col), drop=FALSE]
  colnames(md) <- c(".sample", ".condition")
  md$.sample <- trimws(as.character(md$.sample))
  md$.condition <- trimws(as.character(md$.condition))
  md <- md[nzchar(md$.sample) & nzchar(md$.condition) & !is.na(md$.sample) & !is.na(md$.condition), , drop=FALSE]

  by_sample <- split(md$.condition, md$.sample)
  sample_condition <- vapply(by_sample, function(z) {
    tab <- sort(table(z), decreasing=TRUE)
    names(tab)[1]
  }, character(1))

  if (any(vapply(by_sample, function(z) length(unique(z)) > 1, logical(1)))) {
    add_warning("Some samples contain multiple condition labels; sample-level CV uses each sample's majority condition.", step="classifier")
  }

  sample_df <- data.frame(sample=names(sample_condition), condition=factor(sample_condition), stringsAsFactors=FALSE)
  class_counts <- table(sample_df$condition)
  if (length(class_counts) < 2 || any(class_counts < min_samples_per_class)) {
    return(list(ok=FALSE, sample_df=sample_df, reason=paste0(
      "Sample-level validation requires at least ", min_samples_per_class,
      " biological samples per condition. Observed: ",
      paste(names(class_counts), as.integer(class_counts), sep="=", collapse=", ")
    )))
  }

  k <- min(5, as.integer(min(class_counts)))
  set.seed(420)
  folds <- caret::createFolds(sample_df$condition, k=k, returnTrain=FALSE)
  list(ok=TRUE, sample_df=sample_df, folds=folds, k=k, reason=NULL)
}

write_classifier_skip <- function(reason) {
  out <- data.frame(
    Status="skipped_trusted_validation",
    Reason=reason,
    RecommendedAction="Use at least 2 independent biological samples per condition for sample-level CV, or run --classifier-validation cell_split for exploratory cell-level behavior only.",
    stringsAsFactors=FALSE
  )
  write_csv(out, file.path(tables_dir, "classifier_validation_status.csv"), row_names=FALSE)
  add_warning(reason, step="classifier")
  invisible(out)
}

run_full_data_feature_prioritization <- function(seu, features, condition_col, selector="auto", top_n=10) {
  features <- intersect(unique(features), rownames(seu))
  if (length(features) < 1) return(NULL)
  expr <- t(as.matrix(get_assay_matrix(seu, assay="RNA", layer="data")[features, , drop=FALSE]))
  x <- make_safe_feature_frame(expr)
  y <- droplevels(as.factor(seu[[condition_col]][,1]))
  ranked <- rank_features_training_only(x, y, method=selector, top_n=max(top_n, length(features)))
  if (nrow(ranked) > 0) {
    ranked$Analysis <- "full-data feature prioritization only; not validation metrics"
    ranked$LeakageSafePerformance <- FALSE
    write_csv(ranked, file.path(tables_dir, "classifier_feature_priority_table.csv"), row_names=FALSE)
  }
  ranked
}



main <- function() {
  species_code <- record_step("resolve_species", {
    species_map(opt$species)
  })

  inputs <- record_step("load_inputs", {
    counts <- read.csv(opt$counts, row.names=1, check.names=FALSE)
    meta <- read.csv(opt$metadata, row.names=1, check.names=FALSE)
    list(counts=counts, meta=meta)
  })

  common_cells <- record_step("validate_inputs", {
    validate_inputs(inputs$counts, inputs$meta, opt$`condition-col`, opt$`celltype-col`)
  })

  counts <- inputs$counts[, common_cells, drop=FALSE]
  meta <- inputs$meta[common_cells, , drop=FALSE]

  if (opt$`celltype-col` %in% colnames(meta) && !(opt$`annotation-col` %in% colnames(meta))) {
    meta[[opt$`annotation-col`]] <- meta[[opt$`celltype-col`]]
  }

  seu <- record_step("create_seurat", {
    obj <- CreateSeuratObject(counts=counts, meta.data=meta)
    obj <- NormalizeData(obj, verbose=FALSE)
    if (opt$`celltype-col` %in% colnames(obj@meta.data)) obj[[opt$`celltype-col`]][,1] <- factor(trimws(as.character(obj[[opt$`celltype-col`]][,1])))
    if (opt$`condition-col` %in% colnames(obj@meta.data)) obj[[opt$`condition-col`]][,1] <- factor(trimws(as.character(obj[[opt$`condition-col`]][,1])))
    obj
  })

  summary <- list(
    app="scrna_api_cli",
    module="scRNA",
    version=VERSION,
    species=species_code,
    id_type=opt$`id-type`,
    n_genes=nrow(seu),
    n_cells=ncol(seu),
    metadata_columns=colnames(seu@meta.data),
    parameters=run_status$parameters,
    statistical_notes=statistical_notes()
  )
  write_json_file(summary, file.path(opt$outdir, "run_summary.json"))

  group_col <- if (opt$`annotation-col` %in% colnames(seu@meta.data)) {
    opt$`annotation-col`
  } else if (opt$`celltype-col` %in% colnames(seu@meta.data)) {
    opt$`celltype-col`
  } else {
    NULL
  }

  if (!opt$`skip-pca-umap`) {
    seu <- record_step("pca_umap", {
      nfeat <- min(2000, nrow(seu))
      seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=nfeat, verbose=FALSE)
      seu <- ScaleData(seu, verbose=FALSE)
      npcs <- max(2, min(opt$`max-pcs`, length(VariableFeatures(seu)), ncol(seu)-1))
      if (npcs < 2) stop("Not enough cells/features to compute PCA.")
      seu <- RunPCA(seu, features=VariableFeatures(seu), npcs=npcs, verbose=FALSE)

      pca_df <- as.data.frame(Embeddings(seu, "pca"))
      pca_df$cell_id <- rownames(pca_df)
      pca_df <- cbind(pca_df, seu@meta.data[rownames(pca_df), , drop=FALSE])
      write_csv(pca_df, file.path(tables_dir, "pca_coordinates.csv"), row_names=FALSE)

      if (!is.null(group_col)) {
        save_plotly(
          DimPlot(seu, reduction="pca", group.by=group_col, label=TRUE) + theme_minimal(),
          file.path(plots_dir, "pca_plot.html")
        )
      }

      if (ncol(seu) >= 3) {
        use_dims <- 1:min(opt$`max-pcs`, ncol(Embeddings(seu, "pca")))
        n_neighbors <- max(2, min(30, ncol(seu) - 1))
        seu <- RunUMAP(seu, reduction="pca", dims=use_dims, n.neighbors=n_neighbors, verbose=FALSE)

        umap_df <- as.data.frame(Embeddings(seu, "umap"))
        umap_df$cell_id <- rownames(umap_df)
        umap_df <- cbind(umap_df, seu@meta.data[rownames(umap_df), , drop=FALSE])
        write_csv(umap_df, file.path(tables_dir, "umap_coordinates.csv"), row_names=FALSE)

        if (!is.null(group_col)) {
          save_plotly(
            DimPlot(seu, reduction="umap", group.by=group_col, label=TRUE) + theme_minimal(),
            file.path(plots_dir, "umap_plot.html")
          )
        }
      } else {
        add_warning("Skipping UMAP: need at least 3 cells.", step="pca_umap")
      }
      seu
    })
  } else {
    mark_skipped("pca_umap", "--skip-pca-umap was set")
  }

  condition_de <- NULL
  celltype_de <- NULL
  condition_only <- NULL

  if (opt$`condition-col` %in% colnames(seu@meta.data)) {
    de_results <- record_step("differential_expression", {
      validate_de_options()
      condition_col <- opt$`condition-col`
      celltype_col <- opt$`celltype-col`
      de_mode <- tolower(trimws(opt$`de-mode`))
      de_scope <- tolower(trimws(opt$`de-scope`))

      condition_de <- NULL
      celltype_de <- NULL
      condition_only <- NULL

      condition_cells <- if (de_scope == "global") {
        colnames(seu)
      } else if (celltype_col %in% colnames(seu@meta.data)) {
        select_condition_cells(seu, de_scope, celltype_col, opt$`target-celltype`)
      } else {
        add_warning("Cell type column missing; falling back to global condition DE scope.", step="differential_expression")
        colnames(seu)
      }

      if (de_mode == "cell") {
        condition_de <- run_cell_findmarkers(
          seu,
          group_col=condition_col,
          group_a=opt$`condition-a`,
          group_b=opt$`condition-b`,
          cells=condition_cells,
          label=paste0("condition DE (", de_scope, ")")
        )
      } else {
        condition_de <- run_pseudobulk_de(
          seu,
          group_col=condition_col,
          group_a=opt$`condition-a`,
          group_b=opt$`condition-b`,
          cells=condition_cells,
          label=paste0("condition pseudobulk DE (", de_scope, ")")
        )
      }

      if (!is.null(condition_de)) {
        write_csv(condition_de, file.path(tables_dir, "condition_de_table.csv"), row_names=TRUE)
      }

      if (celltype_col %in% colnames(seu@meta.data)) {
        if (de_mode == "cell") {
          celltype_de <- run_cell_findmarkers(
            seu,
            group_col=celltype_col,
            group_a=opt$`celltype-a`,
            group_b=opt$`celltype-b`,
            cells=NULL,
            label="cell type DE"
          )
        } else {
          celltype_de <- run_pseudobulk_de(
            seu,
            group_col=celltype_col,
            group_a=opt$`celltype-a`,
            group_b=opt$`celltype-b`,
            cells=NULL,
            label="cell type pseudobulk DE"
          )
        }
        if (!is.null(celltype_de)) {
          write_csv(celltype_de, file.path(tables_dir, "celltype_de_table.csv"), row_names=TRUE)
        }
      } else {
        add_warning("Skipping cell type DE: cell type column is missing.", step="differential_expression")
      }

      if (!is.null(condition_de) && !is.null(celltype_de)) {
        cond_genes <- rownames(subset(condition_de, p_val_adj < opt$`de-fdr` & abs(avg_log2FC) >= opt$`logfc-threshold`))
        celltype_genes <- rownames(subset(celltype_de, p_val_adj < opt$`de-fdr` & abs(avg_log2FC) >= opt$`logfc-threshold`))
        only_cond <- setdiff(cond_genes, celltype_genes)
        condition_only <- condition_de[rownames(condition_de) %in% only_cond, , drop=FALSE]
        condition_only$condition_only_mode <- rep(opt$`condition-only-mode`, nrow(condition_only))
        write_csv(condition_only, file.path(tables_dir, "condition_only_de_table.csv"), row_names=TRUE)
        if (nrow(condition_only) == 0) {
          add_warning(sprintf("condition_only_de_table.csv is empty at FDR < %.4g and abs(log2FC) >= %.4g.", opt$`de-fdr`, opt$`logfc-threshold`), step="differential_expression")
        } else {
          volcano <- condition_only
          volcano$log10p <- -log10(volcano$p_val_adj + 1e-300)
          volcano$gene <- rownames(volcano)
          write_csv(volcano, file.path(tables_dir, "condition_only_volcano_table.csv"), row_names=FALSE)
          vp <- ggplot(volcano, aes(x=avg_log2FC, y=log10p, text=gene)) +
            geom_point(alpha=0.6) +
            theme_minimal() +
            labs(title="Volcano Plot: Condition-Only DE Genes", x="Log2 Fold Change", y="-log10 adjusted p-value")
          save_plotly(vp, file.path(plots_dir, "condition_only_volcano.html"))
        }
      } else {
        add_warning("Condition-only filtering skipped because condition DE or cell type DE was unavailable.", step="differential_expression")
      }

      list(condition_de=condition_de, celltype_de=celltype_de, condition_only=condition_only, seu=seu)
    })
    condition_de <- de_results$condition_de
    celltype_de <- de_results$celltype_de
    condition_only <- de_results$condition_only
    seu <- de_results$seu
  } else {
    mark_skipped("differential_expression", "Required condition metadata column is missing")
  }

  if (!opt$`skip-enrichment` && !is.null(condition_only) && nrow(condition_only) > 0) {
    record_step("enrichment", {
      sources <- map_gprof_sources(opt$`gprof-sources`)
      enrich_sets <- list(
        all=rownames(condition_only),
        up=rownames(condition_only[condition_only$avg_log2FC > 0, , drop=FALSE]),
        down=rownames(condition_only[condition_only$avg_log2FC < 0, , drop=FALSE])
      )
      for (nm in names(enrich_sets)) {
        tbl <- do_enrichment(enrich_sets[[nm]], species_code, opt$`enrich-backend`, sources, opt$`enrichr-db`)
        if (!is.null(tbl)) {
          write_csv(tbl, file.path(tables_dir, paste0("enrichment_", nm, ".csv")), row_names=FALSE)
          edges <- make_edges(tbl)
          write_csv(edges, file.path(tables_dir, paste0("enrichment_", nm, "_edges.csv")), row_names=FALSE)
          keep <- tbl[is.finite(tbl$Adjusted.P.value) & tbl$Adjusted.P.value > 0, , drop=FALSE]
          if (nrow(keep) > 0) {
            keep$mlog10 <- -log10(keep$Adjusted.P.value)
            keep <- head(keep[order(keep$Adjusted.P.value), ], 10)
            bp <- ggplot(keep, aes(x=mlog10, y=reorder(Term, mlog10), text=Genes)) +
              geom_col() +
              theme_minimal() +
              labs(title=paste("Top enriched terms:", nm), x="-log10(FDR)", y="")
            save_plotly(bp, file.path(plots_dir, paste0("enrichment_", nm, "_barplot.html")))
          }
        } else {
          add_warning(paste("No enrichment terms returned for", nm), step="enrichment")
        }
      }
      TRUE
    })
  } else if (opt$`skip-enrichment`) {
    mark_skipped("enrichment", "--skip-enrichment was set")
  } else {
    mark_skipped("enrichment", "No condition-only DE genes available")
  }

  if (!opt$`skip-classifier` && !is.null(condition_only) && nrow(condition_only) >= 3 && opt$`condition-col` %in% colnames(seu@meta.data)) {
    record_step("classifier", local({
      validate_classifier_options()

      validation <- tolower(trimws(opt$`classifier-validation`))
      model_choice <- tolower(trimws(opt$`classifier-model`))
      selector_choice <- tolower(trimws(opt$`classifier-selector`))
      feature_source <- tolower(trimws(opt$`classifier-feature-source`))

      labels <- droplevels(as.factor(seu[[opt$`condition-col`]][,1]))
      if (length(levels(labels)) < 2 || min(table(labels)) < 2) {
        write_classifier_skip("Need at least 2 condition groups with at least 2 cells each.")
        return(FALSE)
      }

      # Full-data feature prioritization is allowed, but it is explicitly not used for validation metrics.
      # This preserves useful RF/logistic feature ranking from DE hits without pretending it is test-set performance.
      run_full_data_feature_prioritization(
        seu=seu,
        features=rownames(condition_only),
        condition_col=opt$`condition-col`,
        selector=selector_choice,
        top_n=max(opt$`top-n-features`, nrow(condition_only))
      )

      if (validation == "none") {
        write_classifier_skip("--classifier-validation none was selected; wrote feature-prioritization table only.")
        return(TRUE)
      }

      prediction_rows <- list()
      importance_rows <- list()

      if (validation == "sample_cv") {
        folds_obj <- sample_cv_folds(
          seu@meta.data,
          sample_col=opt$`sample-col`,
          condition_col=opt$`condition-col`,
          min_samples_per_class=opt$`classifier-min-samples-per-class`
        )

        if (!isTRUE(folds_obj$ok)) {
          write_classifier_skip(folds_obj$reason)
          return(TRUE)
        }

        sample_df <- folds_obj$sample_df

        for (fold_i in seq_along(folds_obj$folds)) {
          test_samples <- sample_df$sample[folds_obj$folds[[fold_i]]]
          train_samples <- setdiff(sample_df$sample, test_samples)

          train_cells <- rownames(seu@meta.data)[trimws(as.character(seu@meta.data[[opt$`sample-col`]])) %in% train_samples]
          test_cells <- rownames(seu@meta.data)[trimws(as.character(seu@meta.data[[opt$`sample-col`]])) %in% test_samples]

          candidates <- select_training_de_candidates(
            seu=seu,
            train_cells=train_cells,
            condition_col=opt$`condition-col`,
            celltype_col=opt$`celltype-col`,
            feature_source=feature_source
          )
          candidates <- intersect(candidates, rownames(seu))
          if (length(candidates) < 3) {
            add_warning(paste0("Fold ", fold_i, " skipped: fewer than 3 training-only DE candidate features."), step="classifier")
            next
          }

          expr_train <- t(as.matrix(get_assay_matrix(seu, assay="RNA", layer="data")[candidates, train_cells, drop=FALSE]))
          expr_test <- t(as.matrix(get_assay_matrix(seu, assay="RNA", layer="data")[candidates, test_cells, drop=FALSE]))
          x_train_all <- make_safe_feature_frame(expr_train)
          fmap <- attr(x_train_all, "feature_map")
          colnames(expr_test) <- fmap$Safe[match(colnames(expr_test), fmap$Gene)]
          x_test_all <- as.data.frame(expr_test, check.names=FALSE)

          y_train <- droplevels(as.factor(seu@meta.data[train_cells, opt$`condition-col`]))
          actual <- droplevels(as.factor(seu@meta.data[test_cells, opt$`condition-col`]))

          ranked <- rank_features_training_only(x_train_all, y_train, method=selector_choice, top_n=opt$`top-n-features`)
          if (nrow(ranked) < 1) {
            add_warning(paste0("Fold ", fold_i, " skipped: feature selector returned no features."), step="classifier")
            next
          }

          keep_safe <- ranked$Safe
          x_train <- x_train_all[, keep_safe, drop=FALSE]
          x_test <- x_test_all[, keep_safe, drop=FALSE]

          pred_obj <- train_predict_classifier(x_train, y_train, x_test, model=model_choice)

          prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
            Fold=fold_i,
            Cell=rownames(x_test),
            Sample=seu@meta.data[rownames(x_test), opt$`sample-col`],
            Actual=as.character(actual),
            Predicted=as.character(pred_obj$pred),
            Probability=pred_obj$prob,
            Model=pred_obj$model,
            Validation="sample_cv",
            stringsAsFactors=FALSE
          )

          ranked$Fold <- fold_i
          ranked$Validation <- "sample_cv"
          ranked$ModelUsed <- pred_obj$model
          importance_rows[[length(importance_rows) + 1]] <- ranked
        }
      }

      if (validation == "cell_split") {
        add_warning("Using cell-level split. This is exploratory only and can be optimistic when cells share biological samples.", step="classifier")

        set.seed(42)
        idx <- caret::createDataPartition(labels, p=0.7, list=FALSE)
        train_cells <- colnames(seu)[idx]
        test_cells <- colnames(seu)[-idx]

        candidates <- select_training_de_candidates(
          seu=seu,
          train_cells=train_cells,
          condition_col=opt$`condition-col`,
          celltype_col=opt$`celltype-col`,
          feature_source=feature_source
        )
        candidates <- intersect(candidates, rownames(seu))
        if (length(candidates) < 3) {
          write_classifier_skip("Cell-level classifier skipped: fewer than 3 training-only DE candidate features.")
          return(TRUE)
        }

        expr_train <- t(as.matrix(get_assay_matrix(seu, assay="RNA", layer="data")[candidates, train_cells, drop=FALSE]))
        expr_test <- t(as.matrix(get_assay_matrix(seu, assay="RNA", layer="data")[candidates, test_cells, drop=FALSE]))
        x_train_all <- make_safe_feature_frame(expr_train)
        fmap <- attr(x_train_all, "feature_map")
        colnames(expr_test) <- fmap$Safe[match(colnames(expr_test), fmap$Gene)]
        x_test_all <- as.data.frame(expr_test, check.names=FALSE)

        y_train <- droplevels(as.factor(seu@meta.data[train_cells, opt$`condition-col`]))
        actual <- droplevels(as.factor(seu@meta.data[test_cells, opt$`condition-col`]))

        ranked <- rank_features_training_only(x_train_all, y_train, method=selector_choice, top_n=opt$`top-n-features`)
        keep_safe <- ranked$Safe
        x_train <- x_train_all[, keep_safe, drop=FALSE]
        x_test <- x_test_all[, keep_safe, drop=FALSE]

        pred_obj <- train_predict_classifier(x_train, y_train, x_test, model=model_choice)

        prediction_rows[[1]] <- data.frame(
          Fold=1,
          Cell=rownames(x_test),
          Sample=if (opt$`sample-col` %in% colnames(seu@meta.data)) seu@meta.data[rownames(x_test), opt$`sample-col`] else NA,
          Actual=as.character(actual),
          Predicted=as.character(pred_obj$pred),
          Probability=pred_obj$prob,
          Model=pred_obj$model,
          Validation="cell_split_exploratory",
          stringsAsFactors=FALSE
        )

        ranked$Fold <- 1
        ranked$Validation <- "cell_split_exploratory"
        ranked$ModelUsed <- pred_obj$model
        importance_rows[[1]] <- ranked
      }

      if (length(prediction_rows) == 0) {
        write_classifier_skip("No classifier validation folds produced predictions.")
        return(TRUE)
      }

      preds <- do.call(rbind, prediction_rows)
      write_csv(preds, file.path(tables_dir, "rf_predictions_table.csv"), row_names=FALSE)
      write_csv(preds, file.path(tables_dir, "classifier_predictions_table.csv"), row_names=FALSE)

      imps <- do.call(rbind, importance_rows)
      write_csv(imps, file.path(tables_dir, "rf_importance_table.csv"), row_names=FALSE)
      write_csv(imps, file.path(tables_dir, "classifier_selected_features_table.csv"), row_names=FALSE)

      pred_factor <- factor(preds$Predicted, levels=levels(labels))
      actual_factor <- factor(preds$Actual, levels=levels(labels))
      conf <- caret::confusionMatrix(pred_factor, actual_factor)

      metrics_df <- as.data.frame(t(conf$byClass))
      metrics_df$Metric <- rownames(metrics_df)
      rownames(metrics_df) <- NULL
      metrics_df$Accuracy <- as.numeric(conf$overall["Accuracy"])
      metrics_df$Kappa <- as.numeric(conf$overall["Kappa"])
      metrics_df$Validation <- unique(preds$Validation)[1]
      metrics_df$FeatureSelection <- paste0("DE candidates selected inside training split; selector=", selector_choice)
      metrics_df$LeakageSafe <- validation == "sample_cv"

      if (length(levels(actual_factor)) == 2 && !all(is.na(preds$Probability))) {
        metrics_df$AUC <- as.numeric(pROC::auc(pROC::roc(actual_factor, preds$Probability, quiet=TRUE)))
      }

      write_csv(metrics_df, file.path(tables_dir, "rf_metrics_table.csv"), row_names=FALSE)
      write_csv(metrics_df, file.path(tables_dir, "classifier_metrics_table.csv"), row_names=FALSE)

      stability <- aggregate(
        list(SelectionCount=imps$Gene),
        by=list(Gene=imps$Gene, Selector=imps$Selector, Validation=imps$Validation),
        FUN=length
      )
      stability <- stability[order(stability$SelectionCount, decreasing=TRUE), , drop=FALSE]
      write_csv(stability, file.path(tables_dir, "classifier_feature_stability_table.csv"), row_names=FALSE)

      TRUE
    }))
  } else if (opt$`skip-classifier`) {
    mark_skipped("classifier", "--skip-classifier was set")
  } else {
    mark_skipped("classifier", "Need at least 3 condition-only DE genes and a condition column")
  }

  if (!opt$`skip-power` && !is.null(condition_only) && nrow(condition_only) > 0 && opt$`condition-col` %in% colnames(seu@meta.data)) {
    record_step("power_analysis", {
      de_tab <- condition_only
      top_genes <- head(rownames(de_tab[order(de_tab$p_val_adj), , drop=FALSE]), 10)
      labels <- seu[[opt$`condition-col`]][,1]
      features <- get_assay_matrix(seu, assay="RNA", layer="data")
      power_rows <- list()
      curve_rows <- list()
      for (gene in top_genes) {
        expr <- as.numeric(features[gene, ])
        lv <- unique(labels)
        if (length(lv) < 2) next
        g1 <- expr[labels == lv[1]]
        g2 <- expr[labels == lv[2]]
        if (length(g1) < 2 || length(g2) < 2 || (var(g1) + var(g2)) == 0) next
        eff <- abs(mean(g1) - mean(g2)) / sqrt((var(g1) + var(g2)) / 2)
        pwr_res <- pwr::pwr.t.test(n=min(length(g1), length(g2)), d=eff, sig.level=0.05, type="two.sample", alternative="two.sided")
        n_seq <- seq(5, max(100, max(length(g1), length(g2)) * 2), by=1)
        powers <- sapply(n_seq, function(n) pwr::pwr.t.test(n=n, d=eff, sig.level=0.05, type="two.sample", alternative="two.sided")$power)
        min_n <- if (any(powers >= 0.8)) n_seq[which(powers >= 0.8)[1]] else NA
        power_rows[[gene]] <- data.frame(Gene=gene, EffectSize=round(eff,3), ObservedN1=length(g1), ObservedN2=length(g2), ObservedPower=round(pwr_res$power,3), N_for_0.8=min_n)
        curve_rows[[gene]] <- data.frame(Gene=gene, SampleSizePerGroup=n_seq, Power=powers)
      }
      if (length(power_rows)) {
        write_csv(do.call(rbind, power_rows), file.path(tables_dir, "power_table.csv"), row_names=FALSE)
      } else {
        add_warning("No valid genes for power analysis after variance/group-size checks.", step="power_analysis")
      }
      if (length(curve_rows)) write_csv(do.call(rbind, curve_rows), file.path(tables_dir, "power_curves_long.csv"), row_names=FALSE)
      TRUE
    })
  } else if (opt$`skip-power`) {
    mark_skipped("power_analysis", "--skip-power was set")
  } else {
    mark_skipped("power_analysis", "No condition-only DE genes available")
  }

  record_step("save_outputs", {
    saveRDS(seu, file.path(objects_dir, "seurat_object.rds"))

    table_specs <- list(
      c("condition_de_table.csv", "Condition DE Table"),
      c("celltype_de_table.csv", "Cell Type DE Table"),
      c("condition_only_de_table.csv", "Condition-Only DE Table"),
      c("condition_only_volcano_table.csv", "Condition-Only Volcano Table"),
      c("pca_coordinates.csv", "PCA Coordinates"),
      c("umap_coordinates.csv", "UMAP Coordinates"),
      c("enrichment_all.csv", "Enrichment: All"),
      c("enrichment_up.csv", "Enrichment: Up"),
      c("enrichment_down.csv", "Enrichment: Down"),
      c("enrichment_all_edges.csv", "Enrichment All Term-Gene Edges"),
      c("enrichment_up_edges.csv", "Enrichment Up Term-Gene Edges"),
      c("enrichment_down_edges.csv", "Enrichment Down Term-Gene Edges"),
      c("rf_importance_table.csv", "Random Forest Importance"),
      c("rf_predictions_table.csv", "Random Forest Predictions"),
      c("rf_metrics_table.csv", "Random Forest Metrics"),
      c("classifier_validation_status.csv", "Classifier Validation Status"),
      c("classifier_feature_priority_table.csv", "Classifier Full-Data Feature Prioritization"),
      c("classifier_predictions_table.csv", "Classifier Predictions"),
      c("classifier_selected_features_table.csv", "Classifier Selected Features by Fold"),
      c("classifier_metrics_table.csv", "Classifier Metrics"),
      c("classifier_feature_stability_table.csv", "Classifier Feature Stability"),
      c("power_table.csv", "Power Table"),
      c("power_curves_long.csv", "Power Curves Long Table")
    )
    plot_specs <- list(
      c("pca_plot.html", "PCA Plot"),
      c("umap_plot.html", "UMAP Plot"),
      c("condition_only_volcano.html", "Condition-Only Volcano Plot"),
      c("enrichment_all_barplot.html", "Enrichment All Barplot"),
      c("enrichment_up_barplot.html", "Enrichment Up Barplot"),
      c("enrichment_down_barplot.html", "Enrichment Down Barplot")
    )

    tables <- Filter(Negate(is.null), lapply(table_specs, function(x) table_info(file.path(tables_dir, x[1]), x[2])))
    plots <- Filter(Negate(is.null), lapply(plot_specs, function(x) plot_info(file.path(plots_dir, x[1]), x[2])))
    objects <- Filter(Negate(is.null), list(object_info(file.path(objects_dir, "seurat_object.rds"), "Seurat Object")))

    manifest <- list(
      module="scRNA",
      app="scrna_api_cli",
      version=VERSION,
      created=as.character(Sys.time()),
      summary=summary,
      parameters=run_status$parameters,
      statistical_notes=statistical_notes(),
      tables=tables,
      plots=plots,
      objects=objects,
      warnings=run_status$warnings,
      steps=run_status$steps
    )
    write_json_file(manifest, file.path(opt$outdir, "manifest.json"))
    TRUE
  })

  run_status$finished <<- as.character(Sys.time())
  run_status$status <<- "complete"
  run_status$files <<- list(
    tables=sort(list.files(tables_dir)),
    plots=sort(list.files(plots_dir, pattern="\\.html$")),
    objects=sort(list.files(objects_dir)),
    manifest="manifest.json"
  )
  write_status()

  # Refresh manifest after final status/files are known, so the dashboard never sees save_outputs as still running.
  manifest_path <- file.path(opt$outdir, "manifest.json")
  if (file.exists(manifest_path)) {
    manifest <- jsonlite::fromJSON(manifest_path, simplifyVector=FALSE)
    manifest$status <- run_status$status
    manifest$finished <- run_status$finished
    manifest$files <- run_status$files
    manifest$warnings <- run_status$warnings
    manifest$steps <- run_status$steps
    manifest$statistical_notes <- statistical_notes()
    write_json_file(manifest, manifest_path)
  }

  if (opt$zip) {
    old <- getwd()
    setwd(dirname(opt$outdir))
    on.exit(setwd(old), add=TRUE)
    zipfile <- paste0(basename(opt$outdir), ".zip")
    utils::zip(zipfile, files=basename(opt$outdir))
    cat(normalizePath(zipfile), "\n")
  }

  cat("scRNA-seq API/CLI run complete: ", normalizePath(opt$outdir), "\n", sep="")
}

tryCatch(
  main(),
  error=function(e) {
    run_status$status <<- "failed"
    run_status$finished <<- as.character(Sys.time())
    run_status$errors[[length(run_status$errors) + 1]] <<- list(step="main", message=conditionMessage(e))
    write_status()
    cat(sprintf("ERROR: %s\n", conditionMessage(e)), file=stderr())
    quit(status=1)
  }
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
  library(gprofiler2)
  library(enrichR)
})

VERSION <- "0.2.0-stat-upgrade"

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) == 1 && is.na(a)) return(b)
  a
}

.run_status <- new.env(parent = emptyenv())
.run_status$x <- NULL
.run_status$outdir <- NULL

write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

write_status <- function() {
  if (!is.null(.run_status$x) && !is.null(.run_status$outdir)) {
    write_json_file(.run_status$x, file.path(.run_status$outdir, "run_status.json"))
  }
}

init_status <- function(outdir, parameters) {
  .run_status$outdir <- outdir
  .run_status$x <- list(
    module = "CRISPR Mixscape",
    version = VERSION,
    started = as.character(Sys.time()),
    finished = NULL,
    status = "running",
    outdir = normalizePath(outdir, mustWork = FALSE),
    parameters = parameters,
    steps = list(),
    warnings = list(),
    errors = list(),
    files = list(tables = list(), plots = list(), objects = list())
  )
  write_status()
}

add_warning <- function(msg, step = NULL) {
  .run_status$x$warnings[[length(.run_status$x$warnings) + 1]] <- list(
    time = as.character(Sys.time()),
    step = step,
    message = msg
  )
  write_status()
}

record_step <- function(name, expr) {
  started <- Sys.time()
  .run_status$x$steps[[name]] <- list(status = "running", started = as.character(started))
  write_status()
  tryCatch({
    value <- force(expr)
    finished <- Sys.time()
    .run_status$x$steps[[name]] <- list(
      status = "complete",
      started = as.character(started),
      finished = as.character(finished),
      seconds = round(as.numeric(difftime(finished, started, units = "secs")), 3)
    )
    write_status()
    value
  }, error = function(e) {
    finished <- Sys.time()
    .run_status$x$status <- "failed"
    .run_status$x$finished <- as.character(finished)
    .run_status$x$steps[[name]] <- list(
      status = "failed",
      started = as.character(started),
      finished = as.character(finished),
      seconds = round(as.numeric(difftime(finished, started, units = "secs")), 3),
      error = conditionMessage(e)
    )
    .run_status$x$errors[[length(.run_status$x$errors) + 1]] <- list(step = name, message = conditionMessage(e))
    write_status()
    stop(e)
  })
}

register_file <- function(path, kind = "tables", label = NULL) {
  if (!file.exists(path)) return(invisible(NULL))
  rel <- tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) path)
  root <- tryCatch(normalizePath(.run_status$outdir, mustWork = FALSE), error = function(e) .run_status$outdir)
  rel <- sub(paste0("^", gsub("([\\\\.\"])", "\\\\\\1", root), "/?"), "", rel)
  info <- list(
    name = tools::file_path_sans_ext(basename(path)),
    path = rel,
    label = label %||% basename(path),
    bytes = file.info(path)$size
  )
  .run_status$x$files[[kind]][[length(.run_status$x$files[[kind]]) + 1]] <- info
  write_status()
  invisible(info)
}

write_csv <- function(x, path, row.names = FALSE, label = NULL) {
  write.csv(x, path, row.names = row.names)
  register_file(path, "tables", label)
  invisible(path)
}

species_map <- function(x) {
  x0 <- x
  x <- tolower(trimws(x))
  switch(x,
    "human" = "hsapiens", "hsapiens" = "hsapiens",
    "mouse" = "mmusculus", "mmusculus" = "mmusculus",
    "fly" = "dmelanogaster", "drosophila" = "dmelanogaster", "dmelanogaster" = "dmelanogaster",
    "zebrafish" = "drerio", "drerio" = "drerio",
    stop(paste0("Unsupported species: ", x0))
  )
}

parse_csv_arg <- function(x) {
  trimws(unlist(strsplit(x %||% "", ",")))[nzchar(trimws(unlist(strsplit(x %||% "", ","))))]
}

.is_assay5 <- function(so, assay) inherits(so[[assay]], "Assay5")

get_assay_matrix <- function(seu, assay = "RNA", layer = "data") {
  tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = layer),
    error = function(e1) Seurat::GetAssayData(seu, assay = assay, slot = layer)
  )
}

.ensure_counts_layer <- function(so, assay = DefaultAssay(so)) {
  if (!assay %in% Assays(so)) return(so)
  has_counts <- tryCatch({
    lyr <- GetAssayData(so, assay = assay, layer = "counts")
    !is.null(lyr) && length(lyr) > 0
  }, error = function(...) FALSE)
  if (!has_counts) {
    dat <- tryCatch(GetAssayData(so, assay = assay, layer = "data"), error = function(...) NULL)
    if (is.null(dat) || nrow(dat) == 0) dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
    if (!is.null(dat) && nrow(dat) > 0) {
      cnt <- Matrix::Matrix(round(pmax(exp(dat) - 1, 0)), sparse = TRUE)
      dimnames(cnt) <- list(rownames(dat), colnames(dat))
    } else {
      rn <- rownames(so[[assay]]); cn <- colnames(so)
      cnt <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE, dimnames = list(rn, cn))
    }
    so <- tryCatch(SetAssayData(so, assay = assay, layer = "counts", new.data = cnt), error = function(e) { so[[assay]]@counts <- cnt; so })
  }
  has_data <- tryCatch({
    lyr <- GetAssayData(so, assay = assay, layer = "data")
    !is.null(lyr) && length(lyr) > 0
  }, error = function(...) FALSE)
  if (!has_data) {
    dat <- tryCatch(GetAssayData(so, assay = assay, slot = "data"), error = function(...) NULL)
    if (!is.null(dat) && length(dat) > 0) {
      so <- tryCatch(SetAssayData(so, assay = assay, layer = "data", new.data = dat), error = function(e) so)
    }
  }
  so
}

.harden_prtb <- function(so) {
  if (!"PRTB" %in% Assays(so)) return(so)
  dat <- tryCatch(GetAssayData(so, assay = "PRTB", layer = "data"), error = function(...) NULL)
  if (is.null(dat) || nrow(dat) == 0) dat <- tryCatch(GetAssayData(so, assay = "PRTB", slot = "data"), error = function(...) NULL)
  rn <- if (!is.null(dat) && nrow(dat) > 0) rownames(dat) else rownames(so[["PRTB"]])
  cn <- colnames(so)
  zeros <- Matrix::Matrix(0, nrow = length(rn), ncol = length(cn), sparse = TRUE, dimnames = list(rn, cn))
  so <- tryCatch(SetAssayData(so, assay = "PRTB", layer = "counts", new.data = zeros), error = function(e) { so[["PRTB"]]@counts <- zeros; so })
  if (!is.null(dat) && nrow(dat) > 0) so <- tryCatch(SetAssayData(so, assay = "PRTB", layer = "data", new.data = dat), error = function(e) so)
  sd <- tryCatch(so[["PRTB"]]@scale.data, error = function(...) NULL)
  if (!is.null(sd) && length(sd) && (is.null(dat) || !identical(rownames(sd), rownames(dat)))) {
    so[["PRTB"]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
  }
  DefaultAssay(so) <- "PRTB"
  so <- ScaleData(so, assay = "PRTB", verbose = FALSE)
  DefaultAssay(so) <- "RNA"
  so
}

validate_inputs <- function(counts, metadata, nt_label = "NT") {
  if (anyDuplicated(rownames(counts))) stop("Counts matrix has duplicated gene IDs. Make row names unique before running.")
  if (anyDuplicated(colnames(counts))) stop("Counts matrix has duplicated cell IDs in column names.")
  if (anyDuplicated(rownames(metadata))) stop("Metadata has duplicated cell IDs in row names.")
  if (nrow(counts) < 2) stop("Counts matrix must contain at least 2 genes.")
  if (ncol(counts) < 4) stop("Counts matrix must contain at least 4 cells.")
  numeric_ok <- vapply(counts, is.numeric, logical(1))
  if (!all(numeric_ok)) stop("Counts matrix contains non-numeric columns. Check CSV formatting.")
  shared <- intersect(colnames(counts), rownames(metadata))
  if (length(shared) < 4) stop("Counts columns and metadata rownames must share at least 4 cell IDs.")
  missing_cols <- setdiff(c("gene", "replicate"), colnames(metadata))
  if (length(missing_cols)) stop(paste("Metadata must contain columns:", paste(missing_cols, collapse = ", ")))
  gene <- toupper(trimws(as.character(metadata[shared, "gene"])))
  rep <- trimws(as.character(metadata[shared, "replicate"]))
  if (!any(gene == toupper(nt_label))) stop(paste0("No '", nt_label, "' cells present in metadata column gene."))
  if (length(unique(gene[nzchar(gene)])) < 2) stop("Need at least two perturbation labels in metadata column gene.")
  if (length(unique(rep[nzchar(rep)])) < 2) add_warning("Only one replicate level detected; cell-level DE will be exploratory and pseudobulk DE will not be possible.", step = "validate_inputs")
  invisible(shared)
}

read_inputs <- function(counts_csv, metadata_csv, nt_label = "NT") {
  counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
  metadata <- read.csv(metadata_csv, row.names = 1, check.names = FALSE)
  shared <- validate_inputs(counts, metadata, nt_label = nt_label)
  list(counts = counts[, shared, drop = FALSE], metadata = metadata[shared, , drop = FALSE])
}

write_input_qc <- function(counts, metadata, outdir) {
  qc <- data.frame(
    Metric = c("Input genes", "Input cells", "Metadata rows", "Shared cells", "Perturbation labels", "Replicates", "NT/control cells"),
    Value = c(nrow(counts), ncol(counts), nrow(metadata), length(intersect(colnames(counts), rownames(metadata))), length(unique(metadata$gene)), length(unique(metadata$replicate)), sum(toupper(as.character(metadata$gene)) == "NT")),
    stringsAsFactors = FALSE
  )
  write_csv(qc, file.path(outdir, "input_qc.csv"), row.names = FALSE, label = "Input QC")
  invisible(qc)
}

statistical_notes <- function(parameters) {
  notes <- c(
    paste0("DE mode: ", parameters$de_mode, ". In auto mode the pipeline tries replicate-aware pseudobulk DE first and falls back to cell-level DE only if pseudobulk requirements are not met."),
    "Mixscape classification is performed at the cell level after perturbation signature calculation; downstream KO-vs-NT expression testing should be interpreted according to the selected DE mode.",
    "Pseudobulk DE aggregates raw RNA counts by biological sample/replicate and Mixscape class, then tests replicate-level groups. This is preferred for inferential claims when enough replicates exist.",
    "Cell-level FindMarkers DE treats cells as observations and is best treated as exploratory when cells come from shared biological replicates.",
    "Adjusted p-values are corrected within each output table, not across every table generated by the workflow.",
    "UMAP/PCA and Mixscape posterior plots are exploratory visualizations, not formal statistical tests.",
    "Enrichment depends on the submitted DE gene list and default service background unless a future custom background option is added."
  )
  notes
}

run_basic <- function(counts_csv, metadata_csv, num_neighbors = 20, max_pcs = 40, nt_label = "NT") {
  x <- record_step("load_and_validate_inputs", {
    read_inputs(counts_csv, metadata_csv, nt_label = nt_label)
  })
  write_input_qc(x$counts, x$metadata, .run_status$outdir)
  so <- record_step("normalize_pca_umap", {
    obj <- CreateSeuratObject(counts = x$counts, meta.data = x$metadata)
    DefaultAssay(obj) <- "RNA"
    obj <- .ensure_counts_layer(obj, "RNA")
    obj <- NormalizeData(obj, verbose = FALSE) |> FindVariableFeatures(verbose = FALSE) |> ScaleData(verbose = FALSE)
    n_cells <- ncol(obj); n_genes <- nrow(obj)
    npcs <- min(as.integer(max_pcs), n_cells - 1, n_genes - 1)
    if (npcs < 2) stop("Not enough genes/cells for PCA.")
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
    umap_neighbors <- min(as.integer(num_neighbors), n_cells - 1)
    if (umap_neighbors < 2) umap_neighbors <- 2
    obj <- RunUMAP(obj, dims = 1:npcs, n.neighbors = umap_neighbors, verbose = FALSE)
    attr(obj, "npcs") <- npcs
    attr(obj, "umap_neighbors") <- umap_neighbors
    obj
  })
  so
}

run_mixscape <- function(so, num_neighbors = 20, min_de_genes = 3, iter_num = 20, nt_class = "NT", min_group_cells = 3, max_pcs = 40) {
  record_step("mixscape", {
    DefaultAssay(so) <- "RNA"
    so$gene <- toupper(trimws(as.character(so$gene)))
    so$replicate <- trimws(as.character(so$replicate))
    nt_class <- toupper(nt_class)
    keep <- rownames(so@meta.data)[nzchar(so$gene) & nzchar(so$replicate)]
    if (length(keep) < 50) stop("Too few valid cells after metadata cleanup.")
    so <- subset(so, cells = keep)
    repeat {
      gc <- as.matrix(table(so$gene, so$replicate))
      tiny <- gc < min_group_cells & gc > 0
      if (!any(tiny)) break
      drop_cells <- rownames(so@meta.data)[mapply(function(g, r) gc[g, r, drop = FALSE] < min_group_cells, so$gene, so$replicate)]
      if (!length(drop_cells)) break
      so <- subset(so, cells = setdiff(Cells(so), drop_cells))
    }
    gc <- as.matrix(table(so$gene, so$replicate))
    reps_with_nt <- colnames(gc)[gc[nt_class, ] > 0]
    if (!length(reps_with_nt)) stop("No replicate contains NT after filtering.")
    so <- subset(so, cells = rownames(so@meta.data)[so$replicate %in% reps_with_nt])
    so <- .ensure_counts_layer(so, "RNA")
    so <- NormalizeData(so, verbose = FALSE) |> FindVariableFeatures(verbose = FALSE)
    so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
    n_cells <- ncol(so); n_genes <- nrow(so)
    npcs <- max(2, min(as.integer(max_pcs), n_cells - 1, n_genes - 1))
    so <- RunPCA(so, npcs = npcs, verbose = FALSE)
    ndims <- max(2, min(npcs, ncol(Embeddings(so, "pca"))))
    gc <- as.matrix(table(so$gene, so$replicate))
    min_group_size_per_rep <- apply(gc, 2, function(col) suppressWarnings(min(col[col > 0])))
    global_min_group_size <- suppressWarnings(min(min_group_size_per_rep[is.finite(min_group_size_per_rep)]))
    nt_n <- sum(so$gene == nt_class)
    safe_neighbors <- max(2, min(as.integer(num_neighbors), global_min_group_size - 1, nt_n - 1))
    Idents(so) <- factor(so$gene)
    so <- CalcPerturbSig(object = so, assay = "RNA", slot = "data", gd.class = "gene", nt.cell.class = nt_class, reduction = "pca", ndims = ndims, num.neighbors = safe_neighbors, new.assay.name = "PRTB", split.by = "replicate", verbose = TRUE)
    if (!"PRTB" %in% Assays(so)) stop("CalcPerturbSig did not create PRTB.")
    so <- .harden_prtb(so)
    so <- RunMixscape(object = so, assay = "PRTB", slot = "scale.data", labels = "gene", nt.class.name = nt_class, min.de.genes = as.integer(min_de_genes), iter.num = as.integer(iter_num), de.assay = "RNA", verbose = FALSE, prtb.type = "KO")
    class_candidates <- intersect(c("mixscape_class", "mixscape_class_global"), colnames(so@meta.data))
    if (!length(class_candidates)) stop("Mixscape finished but no class column was added.")
    cc_div <- sapply(class_candidates, function(x) length(na.omit(unique(so@meta.data[[x]]))))
    class_col <- class_candidates[which.max(cc_div)]
    so@meta.data$mixscape_class_active <- so@meta.data[[class_col]]
    post_candidates <- intersect(c("mixscape_class_p_ko", "mixscape_class_global_p_ko"), colnames(so@meta.data))
    if (length(post_candidates)) so@meta.data$mixscape_p_active <- so@meta.data[[post_candidates[1]]]
    attr(so, "mixscape_class_col") <- class_col
    attr(so, "safe_neighbors") <- safe_neighbors
    attr(so, "npcs_mixscape") <- npcs
    so
  })
}

summary_stats <- function(so) {
  meta <- so@meta.data
  cls <- meta$mixscape_class_active %||% meta$mixscape_class %||% rep(NA_character_, nrow(meta))
  data.frame(
    Metric = c("Total Cells", "Total Genes", "KO Cells", "NT Cells", "NP Cells", "Replicates", "sgRNAs", "Variable Genes", "Mixscape class column", "Safe neighbors"),
    Value = c(nrow(meta), nrow(so), sum(grepl(" KO$|^KO$", cls), na.rm = TRUE), sum(cls == "NT", na.rm = TRUE), sum(grepl(" NP$|^NP$", cls), na.rm = TRUE), length(unique(meta$replicate %||% NA)), length(unique(meta$guide_ID %||% NA)), length(VariableFeatures(so)), attr(so, "mixscape_class_col") %||% "mixscape_class_active", attr(so, "safe_neighbors") %||% NA),
    stringsAsFactors = FALSE
  )
}

write_class_tables <- function(so, outdir) {
  meta <- so@meta.data
  cls_tbl <- as.data.frame(table(gene = meta$gene, replicate = meta$replicate, mixscape_class = meta$mixscape_class_active), stringsAsFactors = FALSE) |>
    dplyr::group_by(gene, replicate) |>
    dplyr::mutate(percent = Freq / sum(Freq)) |>
    dplyr::ungroup()
  write_csv(cls_tbl, file.path(outdir, "mixscape_class_counts_by_replicate.csv"), row.names = FALSE, label = "Mixscape class counts by replicate")
  invisible(cls_tbl)
}

ko_labels <- function(so) {
  cls <- as.character(so@meta.data$mixscape_class_active)
  sort(unique(cls[!is.na(cls) & (cls == "KO" | grepl(" KO$", cls))]))
}

annotate_de_table <- function(df, comparison, de_mode, de_method, pseudobulk_samples_1 = NA, pseudobulk_samples_2 = NA) {
  if (is.null(df)) df <- data.frame()
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  n <- nrow(df)
  df$comparison <- rep(comparison, n)
  df$de_mode <- rep(de_mode, n)
  df$de_method <- rep(de_method, n)
  if (!is.na(pseudobulk_samples_1) || !is.na(pseudobulk_samples_2)) {
    df$pseudobulk_samples_1 <- rep(pseudobulk_samples_1, n)
    df$pseudobulk_samples_2 <- rep(pseudobulk_samples_2, n)
  }
  df
}

run_cell_de <- function(so, ko_label = NULL, nt_label = "NT", min_pct = 0.1, logfc_threshold = 0, test_use = "wilcox") {
  choices <- ko_labels(so)
  if (!length(choices)) return(data.frame(message = "No KO class detected."))
  ko <- ko_label %||% choices[1]
  if (!(ko %in% choices)) ko <- choices[1]
  tab <- FindMarkers(so, ident.1 = ko, ident.2 = nt_label, group.by = "mixscape_class_active", assay = "RNA", min.pct = min_pct, logfc.threshold = logfc_threshold, test.use = test_use, verbose = FALSE)
  annotate_de_table(tab, comparison = paste0(ko, " vs ", nt_label), de_mode = "cell", de_method = paste0("Seurat FindMarkers test.use=", test_use, "; cell-level exploratory"))
}

run_pseudobulk_limma_voom <- function(pb_counts, pb_meta, contrast, label) {
  if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
    stop("Pseudobulk limma_voom requires edgeR and limma. Install with BiocManager::install(c('edgeR','limma')).")
  }
  group <- factor(pb_meta$.group, levels = contrast)
  dge <- edgeR::DGEList(counts = pb_counts, group = group)
  keep_genes <- edgeR::filterByExpr(dge, group = group)
  if (sum(keep_genes) < 2) stop("Fewer than 2 genes remained after edgeR::filterByExpr.")
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- make.names(levels(group))
  residual_df <- ncol(dge$counts) - qr(design)$rank
  if (!is.finite(residual_df) || residual_df <= 0) stop("limma/voom has no residual degrees of freedom. Add biological replicates or use --de-mode cell for exploratory analysis.")
  contrast_expr <- paste0(make.names(contrast[1]), "-", make.names(contrast[2]))
  v <- limma::voom(dge, design = design, plot = FALSE)
  cm <- limma::makeContrasts(contrasts = contrast_expr, levels = design)
  fit <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit, cm)
  fit <- limma::eBayes(fit, trend = TRUE)
  tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "P")
  genes <- rownames(tt)
  g1 <- which(group == contrast[1]); g2 <- which(group == contrast[2])
  df <- data.frame(
    p_val = as.numeric(tt$P.Value),
    avg_log2FC = as.numeric(tt$logFC),
    pct.1 = if (length(genes)) rowMeans(pb_counts[genes, g1, drop = FALSE] > 0) else numeric(0),
    pct.2 = if (length(genes)) rowMeans(pb_counts[genes, g2, drop = FALSE] > 0) else numeric(0),
    p_val_adj = as.numeric(tt$adj.P.Val),
    gene = genes,
    AveExpr = if ("AveExpr" %in% colnames(tt)) as.numeric(tt$AveExpr) else NA_real_,
    t = if ("t" %in% colnames(tt)) as.numeric(tt$t) else NA_real_,
    B = if ("B" %in% colnames(tt)) as.numeric(tt$B) else NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  df <- df[order(df$p_val_adj, df$p_val), , drop = FALSE]
  rownames(df) <- df$gene
  annotate_de_table(df, comparison = label, de_mode = "pseudobulk", de_method = "edgeR::filterByExpr + TMM + limma::voom + limma::eBayes(trend=TRUE)", pseudobulk_samples_1 = length(g1), pseudobulk_samples_2 = length(g2))
}

run_pseudobulk_de <- function(so, ko_label = NULL, nt_label = "NT", sample_col = "replicate", min_reps = 2, method = "limma_voom") {
  choices <- ko_labels(so)
  if (!length(choices)) stop("No KO class detected.")
  ko <- ko_label %||% choices[1]
  if (!(ko %in% choices)) ko <- choices[1]
  if (!(sample_col %in% colnames(so@meta.data))) stop(paste0("Pseudobulk requires sample column: ", sample_col))
  md <- so@meta.data
  md$.sample <- trimws(as.character(md[[sample_col]]))
  md$.group <- trimws(as.character(md$mixscape_class_active))
  md <- md[md$.group %in% c(ko, nt_label) & nzchar(md$.sample), , drop = FALSE]
  if (nrow(md) < 4) stop("Too few cells in selected KO/NT classes for pseudobulk.")
  md$.pb_id <- paste(md$.sample, md$.group, sep = "__")
  raw_counts <- as.matrix(get_assay_matrix(so, assay = "RNA", layer = "counts")[, rownames(md), drop = FALSE])
  pb_counts <- t(rowsum(t(raw_counts), group = md$.pb_id, reorder = FALSE))
  pb_meta <- unique(md[, c(".pb_id", ".sample", ".group"), drop = FALSE])
  rownames(pb_meta) <- pb_meta$.pb_id
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  rep_counts <- table(pb_meta$.group)
  if (any(rep_counts[c(ko, nt_label)] < min_reps)) {
    stop(paste0("Pseudobulk requires at least ", min_reps, " samples per class. Observed: ", paste(names(rep_counts), as.integer(rep_counts), sep = "=", collapse = ", ")))
  }
  run_pseudobulk_limma_voom(pb_counts, pb_meta, contrast = c(ko, nt_label), label = paste0(ko, " vs ", nt_label))
}

run_de <- function(so, ko_label = NULL, nt_label = "NT", de_mode = "auto", pseudobulk_method = "limma_voom", sample_col = "replicate", pseudobulk_min_replicates = 2, min_pct = 0.1, logfc_threshold = 0, test_use = "wilcox") {
  record_step("differential_expression", {
    de_mode <- tolower(trimws(de_mode))
    if (!(de_mode %in% c("auto", "pseudobulk", "cell"))) stop("--de-mode must be one of: auto, pseudobulk, cell")
    if (de_mode %in% c("auto", "pseudobulk")) {
      pb <- tryCatch(run_pseudobulk_de(so, ko_label, nt_label, sample_col, pseudobulk_min_replicates, pseudobulk_method), error = function(e) {
        if (de_mode == "pseudobulk") stop(e)
        add_warning(paste("Pseudobulk DE unavailable; falling back to cell-level exploratory DE:", conditionMessage(e)), step = "differential_expression")
        NULL
      })
      if (!is.null(pb)) return(pb)
    }
    run_cell_de(so, ko_label, nt_label, min_pct = min_pct, logfc_threshold = logfc_threshold, test_use = test_use)
  })
}

pick_de_genes <- function(tab, mode = c("all", "up", "down"), fdr = 0.05, logfc_threshold = 0) {
  mode <- match.arg(mode)
  if (is.null(tab) || !nrow(tab) || "message" %in% names(tab)) return(character(0))
  fc_col <- if ("avg_log2FC" %in% names(tab)) "avg_log2FC" else if ("avg_logFC" %in% names(tab)) "avg_logFC" else NA
  padj_col <- if ("p_val_adj" %in% names(tab)) "p_val_adj" else if ("adj.P.Val" %in% names(tab)) "adj.P.Val" else if ("p_adj" %in% names(tab)) "p_adj" else NA
  if (is.na(fc_col) || is.na(padj_col)) return(character(0))
  keep <- tab[is.finite(tab[[padj_col]]) & tab[[padj_col]] < fdr & abs(tab[[fc_col]]) >= logfc_threshold, , drop = FALSE]
  if (!nrow(keep)) return(character(0))
  if (mode == "up") keep <- keep[keep[[fc_col]] > 0, , drop = FALSE]
  if (mode == "down") keep <- keep[keep[[fc_col]] < 0, , drop = FALSE]
  unique(rownames(keep) %||% keep$gene)
}

fmt_enrich_for_table <- function(res) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  if (all(c("term_name", "term_size", "intersection_size", "p_value") %in% colnames(res))) {
    adj <- if ("p_adjusted" %in% colnames(res)) res$p_adjusted else p.adjust(res$p_value, "fdr")
    df <- data.frame(Term = paste0(res$term_name, " (", res$source, ")"), Overlap = paste0(res$intersection_size, "/", res$term_size), P.value = as.numeric(res$p_value), Adjusted.P.value = as.numeric(adj), Genes = if ("intersection" %in% colnames(res)) res$intersection else NA_character_, stringsAsFactors = FALSE, check.names = FALSE)
    return(df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE])
  }
  if (all(c("Term", "Overlap", "P.value", "Adjusted.P.value") %in% colnames(res))) {
    df <- data.frame(Term = res$Term, Overlap = res$Overlap, P.value = as.numeric(res$P.value), Adjusted.P.value = as.numeric(res$Adjusted.P.value), Genes = if ("Genes" %in% colnames(res)) res$Genes else NA_character_, stringsAsFactors = FALSE, check.names = FALSE)
    return(df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE])
  }
  NULL
}

run_enrichment_set <- function(genes, species, backend, gprof_sources, enrichr_dbs) {
  genes <- unique(genes[!grepl("^BG\\d+$", genes)])
  if (length(genes) < 3) return(data.frame(message = "Fewer than 3 DE genes passed thresholds."))
  backend <- tolower(backend)
  if (backend == "auto") backend <- if (species == "hsapiens") "enrichr" else "gprof"
  if (backend == "gprof") {
    gp <- tryCatch(gprofiler2::gost(query = genes, organism = species, sources = gprof_sources, correction_method = "fdr", evcodes = TRUE), error = function(e) { add_warning(paste("g:Profiler failed:", conditionMessage(e)), step = "enrichment"); NULL })
    tbl <- fmt_enrich_for_table(gp$result %||% NULL)
    if (is.null(tbl)) return(data.frame(message = "No enrichment terms returned."))
    return(tbl)
  }
  if (backend == "enrichr") {
    out <- tryCatch(enrichR::enrichr(genes, enrichr_dbs), error = function(e) { add_warning(paste("Enrichr failed:", conditionMessage(e)), step = "enrichment"); NULL })
    if (is.null(out)) return(data.frame(message = "No enrichment terms returned."))
    tbl <- do.call(rbind, lapply(names(out), function(db) { df <- out[[db]]; if (is.null(df) || nrow(df) == 0) return(NULL); df$source <- db; df }))
    if (is.null(tbl) || nrow(tbl) == 0) return(data.frame(message = "No enrichment terms returned."))
    std <- fmt_enrich_for_table(tbl)
    std$Term <- paste0(std$Term, " (", tbl$source[match(std$Term, tbl$Term)], ")")
    return(std)
  }
  stop("--enrich-backend must be one of: auto, gprof, enrichr")
}

run_enrichment <- function(de_tab, outdir, species, backend, gprof_sources, enrichr_dbs, de_fdr = 0.05, logfc_threshold = 0) {
  record_step("enrichment", {
    for (mode in c("all", "up", "down")) {
      genes <- pick_de_genes(de_tab, mode = mode, fdr = de_fdr, logfc_threshold = logfc_threshold)
      tbl <- run_enrichment_set(genes, species, backend, gprof_sources, enrichr_dbs)
      write_csv(tbl, file.path(outdir, paste0("enrichment_", mode, ".csv")), row.names = FALSE, label = paste("Enrichment", mode))
    }
    TRUE
  })
}

make_plots <- function(so, de_tab, outdir) {
  record_step("plots", {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    p1 <- DimPlot(so, group.by = "replicate", pt.size = 0.2, reduction = "umap", repel = TRUE)
    ggsave(file.path(outdir, "umap_replicate.png"), p1, width = 8, height = 6, dpi = 150)
    register_file(file.path(outdir, "umap_replicate.png"), "plots", "UMAP by replicate")
    if ("Phase" %in% colnames(so@meta.data)) {
      p2 <- DimPlot(so, group.by = "Phase", pt.size = 0.2, reduction = "umap", repel = TRUE)
      ggsave(file.path(outdir, "umap_phase.png"), p2, width = 8, height = 6, dpi = 150)
      register_file(file.path(outdir, "umap_phase.png"), "plots", "UMAP by phase")
    }
    if ("crispr" %in% colnames(so@meta.data)) {
      p3 <- DimPlot(so, group.by = "crispr", pt.size = 0.2, reduction = "umap", repel = TRUE)
      ggsave(file.path(outdir, "umap_crispr.png"), p3, width = 8, height = 6, dpi = 150)
      register_file(file.path(outdir, "umap_crispr.png"), "plots", "UMAP by CRISPR label")
    }
    if ("mixscape_p_active" %in% colnames(so@meta.data)) {
      p4 <- VlnPlot(so, features = "mixscape_p_active", group.by = "mixscape_class_active", pt.size = 0) + ggtitle("Mixscape posterior") + NoLegend()
      ggsave(file.path(outdir, "mixscape_posterior_violin.png"), p4, width = 10, height = 6, dpi = 150)
      register_file(file.path(outdir, "mixscape_posterior_violin.png"), "plots", "Mixscape posterior violin")
    }
    TRUE
  })
}

write_bundle <- function(so, de_tab, outdir, species, parameters, skip_enrichment = FALSE, skip_plots = FALSE) {
  record_step("write_outputs", {
    write_csv(so@meta.data, file.path(outdir, "cell_metadata_mixscape.csv"), row.names = TRUE, label = "Cell metadata with Mixscape classes")
    write_csv(summary_stats(so), file.path(outdir, "summary_stats.csv"), row.names = FALSE, label = "Summary stats")
    write_csv(de_tab, file.path(outdir, "ko_genes.csv"), row.names = TRUE, label = "KO vs NT DE genes")
    write_class_tables(so, outdir)
    notes <- data.frame(Note = statistical_notes(parameters), stringsAsFactors = FALSE)
    write_csv(notes, file.path(outdir, "statistical_notes.csv"), row.names = FALSE, label = "Statistical notes")
    summary <- list(app = "crispr_mixscape_api_cli", module = "CRISPR Mixscape", version = VERSION, species = species, n_genes = nrow(so), n_cells = ncol(so), metadata_columns = colnames(so@meta.data), parameters = parameters, statistical_notes = statistical_notes(parameters))
    write_json_file(summary, file.path(outdir, "run_summary.json"))
    register_file(file.path(outdir, "run_summary.json"), "tables", "Run summary JSON")
    saveRDS(so, file.path(outdir, "mixscape_seurat.rds"))
    register_file(file.path(outdir, "mixscape_seurat.rds"), "objects", "Seurat object")
    TRUE
  })
  if (!skip_enrichment) run_enrichment(de_tab, outdir, species, parameters$enrich_backend, parse_csv_arg(parameters$gprof_sources), parse_csv_arg(parameters$enrichr_dbs), parameters$de_fdr, parameters$logfc_threshold) else add_warning("Skipped enrichment because --skip-enrichment was set.", step = "enrichment")
  if (!skip_plots) make_plots(so, de_tab, file.path(outdir, "plots")) else add_warning("Skipped plots because --skip-plots was set.", step = "plots")
  manifest <- list(status = "complete", module = "CRISPR Mixscape", version = VERSION, files = list.files(outdir, recursive = TRUE), class_column = attr(so, "mixscape_class_col") %||% "mixscape_class_active", safe_neighbors = attr(so, "safe_neighbors") %||% NA, de_mode_used = if ("de_mode" %in% colnames(de_tab)) unique(de_tab$de_mode)[1] else NA)
  write_json_file(manifest, file.path(outdir, "manifest.json"))
  register_file(file.path(outdir, "manifest.json"), "tables", "Output manifest")
  manifest
}

run_pipeline <- function(counts_csv, metadata_csv, outdir, species = "hsapiens", nt_label = "NT", ko_label = NULL, num_neighbors = 20, max_pcs = 40, min_group_cells = 3, min_de_genes = 3, iter_num = 20, de_mode = "auto", pseudobulk_method = "limma_voom", sample_col = "replicate", pseudobulk_min_replicates = 2, de_fdr = 0.05, logfc_threshold = 0, min_pct = 0.1, test_use = "wilcox", enrich_backend = "auto", gprof_sources = "GO:BP,GO:MF,GO:CC,REAC,KEGG", enrichr_dbs = "GO_Biological_Process_2021,GO_Molecular_Function_2021,GO_Cellular_Component_2021,Reactome_2022,KEGG_2021_Human", skip_enrichment = FALSE, skip_plots = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  species <- species_map(species)
  parameters <- list(species = species, nt_label = nt_label, ko_label = ko_label, neighbors = num_neighbors, max_pcs = max_pcs, min_group_cells = min_group_cells, min_de_genes = min_de_genes, iter_num = iter_num, de_mode = de_mode, pseudobulk_method = pseudobulk_method, sample_col = sample_col, pseudobulk_min_replicates = pseudobulk_min_replicates, de_fdr = de_fdr, logfc_threshold = logfc_threshold, min_pct = min_pct, test_use = test_use, enrich_backend = enrich_backend, gprof_sources = gprof_sources, enrichr_dbs = enrichr_dbs)
  init_status(outdir, parameters)
  so <- run_basic(counts_csv, metadata_csv, num_neighbors = num_neighbors, max_pcs = max_pcs, nt_label = nt_label)
  so <- run_mixscape(so, num_neighbors = num_neighbors, min_de_genes = min_de_genes, iter_num = iter_num, nt_class = nt_label, min_group_cells = min_group_cells, max_pcs = max_pcs)
  de_tab <- run_de(so, ko_label = ko_label, nt_label = nt_label, de_mode = de_mode, pseudobulk_method = pseudobulk_method, sample_col = sample_col, pseudobulk_min_replicates = pseudobulk_min_replicates, min_pct = min_pct, logfc_threshold = logfc_threshold, test_use = test_use)
  manifest <- write_bundle(so, de_tab, outdir = outdir, species = species, parameters = parameters, skip_enrichment = skip_enrichment, skip_plots = skip_plots)
  .run_status$x$status <- "complete"
  .run_status$x$finished <- as.character(Sys.time())
  write_status()
  manifest
}

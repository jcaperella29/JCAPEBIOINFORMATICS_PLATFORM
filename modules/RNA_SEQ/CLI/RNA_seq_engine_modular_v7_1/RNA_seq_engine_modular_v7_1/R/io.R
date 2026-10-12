# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

read_matrix_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    read.csv(path, row.names = 1, check.names = FALSE)
  } else {
    read.table(path, header = TRUE, row.names = 1, check.names = FALSE, sep = "", quote = "\"", comment.char = "")
  }
}

strip_ens_version <- function(x) sub("\\..*", "", x)

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


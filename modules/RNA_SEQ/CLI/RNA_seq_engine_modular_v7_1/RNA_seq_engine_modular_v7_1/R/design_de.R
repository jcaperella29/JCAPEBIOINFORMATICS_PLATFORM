# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

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


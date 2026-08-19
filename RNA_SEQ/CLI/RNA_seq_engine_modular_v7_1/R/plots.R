# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

save_plotly_html <- function(plot_obj, path) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Package htmlwidgets is required.")
  htmlwidgets::saveWidget(plot_obj, file = path, selfcontained = TRUE)
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


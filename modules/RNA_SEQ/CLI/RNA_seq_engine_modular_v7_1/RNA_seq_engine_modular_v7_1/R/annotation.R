# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

bm_dataset_for <- function(sp) {
  switch(
    sp,
    "hsapiens"      = "hsapiens_gene_ensembl",
    "mmusculus"     = "mmusculus_gene_ensembl",
    "drerio"        = "drerio_gene_ensembl",
    "dmelanogaster" = "dmelanogaster_gene_ensembl",
    "hsapiens_gene_ensembl"
  )
}

annotate_genes <- function(ensembl_ids, species_code) {
  ens_ids <- unique(strip_ens_version(ensembl_ids))

  # Prefer an offline Bioconductor annotation database for human genes.
  # This avoids runtime dependence on the Ensembl BioMart service.
  if (species_code == "hsapiens") {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      stop("Package 'AnnotationDbi' is required for human gene annotation.")
    }
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is required for human gene annotation.")
    }

    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = ens_ids,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )

    descriptions <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = ens_ids,
      keytype = "ENSEMBL",
      column = "GENENAME",
      multiVals = "first"
    )

    return(data.frame(
      ensembl_gene_id = ens_ids,
      symbol = unname(symbols[ens_ids]),
      description = unname(descriptions[ens_ids]),
      stringsAsFactors = FALSE
    ))
  }

  # Existing BioMart behavior remains available for other supported species.
  dataset <- bm_dataset_for(species_code)
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  ann <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = ens_ids,
    mart = mart
  )
  names(ann) <- c("ensembl_gene_id", "symbol", "description")
  ann
}

annotate_de_table <- function(de_table, species_code, skip_annotation = FALSE) {
  de_table$Ensembl_IDs <- strip_ens_version(de_table$Ensembl_IDs)
  de_table$.row_order <- seq_len(nrow(de_table))
  if (skip_annotation) {
    de_table$symbol <- de_table$Ensembl_IDs
    de_table$description <- NA_character_
    de_table <- de_table[order(de_table$.row_order), , drop = FALSE]
    de_table$.row_order <- NULL
    return(de_table)
  }

  log_msg("Annotating genes with biomaRt...")
  ann <- tryCatch(
    annotate_genes(de_table$Ensembl_IDs, species_code),
    error = function(e) {
      log_msg("Gene annotation failed; continuing with Ensembl IDs as symbols: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(ann) || nrow(ann) == 0) {
    de_table$symbol <- de_table$Ensembl_IDs
    de_table$description <- NA_character_
    de_table <- de_table[order(de_table$.row_order), , drop = FALSE]
    de_table$.row_order <- NULL
    return(de_table)
  }

  # biomaRt can return multiple rows for one Ensembl gene. Collapse to exactly
  # one annotation row per gene before joining so annotation cannot duplicate DE rows.
  clean_nonempty <- function(x) {
    x <- unique(trimws(as.character(x)))
    sort(x[!is.na(x) & nzchar(x)])
  }
  first_nonempty <- function(x, fallback = NA_character_) {
    x <- clean_nonempty(x)
    if (length(x) == 0) return(fallback)
    x[[1]]
  }
  collapse_nonempty <- function(x, fallback = NA_character_) {
    x <- clean_nonempty(x)
    if (length(x) == 0) return(fallback)
    paste(x, collapse = " | ")
  }
  ann_ids <- clean_nonempty(ann$ensembl_gene_id)
  ann_collapsed <- do.call(rbind, lapply(ann_ids, function(id) {
    z <- ann[ann$ensembl_gene_id == id, , drop = FALSE]
    data.frame(
      ensembl_gene_id = id,
      symbol = first_nonempty(z$symbol, fallback = id),
      description = collapse_nonempty(z$description, fallback = NA_character_),
      stringsAsFactors = FALSE
    )
  }))

  n_before <- nrow(de_table)
  out <- merge(de_table, ann_collapsed, by.x = "Ensembl_IDs", by.y = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$.row_order), , drop = FALSE]
  out$.row_order <- NULL
  out$symbol[is.na(out$symbol) | out$symbol == ""] <- out$Ensembl_IDs[is.na(out$symbol) | out$symbol == ""]
  if (nrow(out) != n_before) {
    stop("Annotation join changed DE row count (", n_before, " -> ", nrow(out), "). Annotation must be one-row-per-tested-gene.")
  }
  out
}


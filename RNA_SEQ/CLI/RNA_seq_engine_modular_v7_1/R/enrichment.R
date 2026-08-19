# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

gp_org_for <- function(sp) {
  switch(sp,
         "hsapiens" = "hsapiens",
         "mmusculus" = "mmusculus",
         "drerio" = "drerio",
         "dmelanogaster" = "dmelanogaster",
         "hsapiens")
}

get_species_enrich_choices <- function(species_code) {
  if (species_code == "hsapiens") {
    c("GO_Biological_Process_2021", "KEGG_2021_Human", "WikiPathway_2021_Human", "Reactome_2022")
  } else {
    c("GO: Biological Process" = "GO:BP", "KEGG" = "KEGG", "Reactome" = "REAC", "WikiPathways" = "WP")
  }
}

default_enrich_db <- function(species_code) unname(get_species_enrich_choices(species_code)[[1]])

is_ensembl_like <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^ENS[A-Z]*G[0-9]+", x)) >= 0.5
}

gprofiler_source_for_db <- function(db) {
  db_txt <- toupper(as.character(db))
  if (db_txt %in% c("GO:BP", "GO_BP", "GOBP") || grepl("GO", db_txt) || grepl("BIOLOGICAL", db_txt)) return("GO:BP")
  if (db_txt %in% c("KEGG") || grepl("KEGG", db_txt)) return("KEGG")
  if (db_txt %in% c("REAC", "REACTOME") || grepl("REACTOME", db_txt)) return("REAC")
  if (db_txt %in% c("WP", "WIKIPATHWAYS") || grepl("WIKI", db_txt)) return("WP")
  "GO:BP"
}

format_gprofiler_enrichment <- function(g, source_label = "g:Profiler") {
  if (is.null(g) || is.null(g$result) || nrow(g$result) == 0) return(NULL)
  df <- g$result
  comb <- -log10(df$p_value) * (df$intersection_size / df$term_size)
  genes_col <- vapply(seq_len(nrow(df)), function(i) {
    x <- df$intersection[[i]]
    if (is.null(x) || length(x) == 0 || all(is.na(x))) return("")
    paste(unique(as.character(x)), collapse = ";")
  }, character(1))
  gene_count_col <- if ("intersection_size" %in% names(df)) df$intersection_size else vapply(strsplit(genes_col, ";", fixed = TRUE), function(x) sum(x != ""), integer(1))
  out <- data.frame(
    Term = df$term_name,
    Genes = genes_col,
    Gene.Count = gene_count_col,
    Adjusted.P.value = df$p_value,
    Combined.Score = comb,
    Source = source_label,
    stringsAsFactors = FALSE
  )
  out[order(out$Adjusted.P.value), ]
}

run_gprofiler_enrichment <- function(gene_list, db, species_code) {
  source <- gprofiler_source_for_db(db)
  org <- gp_org_for(species_code)
  log_msg("Running g:Profiler enrichment fallback/source ", source, " with ", length(gene_list), " identifiers...")
  g <- tryCatch(
    gprofiler2::gost(query = gene_list, organism = org, sources = source, correction_method = "g_SCS", evcodes = TRUE),
    error = function(e) {
      log_msg("g:Profiler enrichment failed for ", source, ": ", conditionMessage(e))
      NULL
    }
  )
  format_gprofiler_enrichment(g, source_label = paste0("g:Profiler ", source))
}

perform_enrichment <- function(gene_list, db, species_code) {
  gene_list <- unique(na.omit(as.character(gene_list)))
  gene_list <- gene_list[gene_list != ""]
  if (length(gene_list) < 2) return(NULL)

  # Production behavior:
  # - Enrichr generally prefers human gene symbols.
  # - If BioMart annotation fails and we only have Ensembl IDs, use g:Profiler,
  #   which can resolve Ensembl IDs directly. This keeps enrichment functional
  #   without making annotation a hard dependency.
  if (species_code == "hsapiens" && !is_ensembl_like(gene_list)) {
    enrichr_res <- tryCatch(enrichR::enrichr(gene_list, databases = db)[[1]], error = function(e) {
      log_msg("Human Enrichr failed for ", db, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
      if (!"Source" %in% names(enrichr_res)) enrichr_res$Source <- paste0("Enrichr ", db)
      return(enrichr_res)
    }
    log_msg("Enrichr returned no usable results; trying g:Profiler fallback.")
    return(run_gprofiler_enrichment(gene_list, db, species_code))
  }

  run_gprofiler_enrichment(gene_list, db, species_code)
}

perform_enrichment_v6 <- function(gene_list, db, species_code, backend = "auto", gprof_sources = c("GO:BP"), enrichr_db = "GO_Biological_Process_2023", background_genes = NULL) {
  gene_list <- unique(na.omit(as.character(gene_list)))
  gene_list <- gene_list[gene_list != ""]
  background_genes <- unique(na.omit(as.character(background_genes)))
  background_genes <- background_genes[background_genes != ""]
  if (length(gene_list) < 2) return(NULL)
  if (length(background_genes) < 2) stop("Enrichment requires the tested-gene background; fewer than two background genes were supplied.")
  backend <- tolower(backend)

  enrichr_ids_ok <- species_code == "hsapiens" &&
    !any(grepl("^ENS[A-Z]*G[0-9]+", gene_list)) &&
    !any(grepl("^ENS[A-Z]*G[0-9]+", background_genes))
  if (backend %in% c("auto", "enrichr") && enrichr_ids_ok) {
    er_db <- if (!is.null(enrichr_db) && nzchar(enrichr_db)) enrichr_db else db
    enrichr_res <- tryCatch(enrichR::enrichr(gene_list, databases = er_db, background = background_genes)[[1]], error = function(e) {
      log_msg("Enrichr failed for ", er_db, " with tested-gene background: ", conditionMessage(e))
      NULL
    })
    if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
      if (!"Source" %in% names(enrichr_res)) enrichr_res$Source <- paste0("Enrichr ", er_db)
      return(enrichr_res)
    }
    if (backend == "enrichr") return(NULL)
    log_msg("Enrichr returned no usable results; trying g:Profiler fallback with the same tested-gene background.")
  }

  if (backend == "enrichr" && !enrichr_ids_ok) {
    log_msg("Enrichr backend requested but tested/query IDs include Ensembl identifiers; using g:Profiler so the full tested-gene background remains valid.")
  }

  org <- gp_org_for(species_code)
  sources <- unique(trimws(gprof_sources))
  sources <- sources[nzchar(sources)]
  if (length(sources) < 1) sources <- gprof_source_for_db(db)
  log_msg("Running g:Profiler enrichment source(s) ", paste(sources, collapse = ","), " with ", length(gene_list), " identifiers and ", length(background_genes), " tested-background identifiers...")
  g <- tryCatch(
    gprofiler2::gost(query = gene_list, organism = org, sources = sources, correction_method = "g_SCS", evcodes = TRUE, custom_bg = background_genes, domain_scope = "custom"),
    error = function(e) {
      log_msg("g:Profiler enrichment failed: ", conditionMessage(e))
      NULL
    }
  )
  format_gprofiler_enrichment(g, source_label = paste0("g:Profiler ", paste(sources, collapse = ",")))
}

normalize_col_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

find_first_col <- function(df, candidates) {
  norm_names <- normalize_col_name(names(df))
  norm_candidates <- normalize_col_name(candidates)
  hit <- match(norm_candidates, norm_names, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) == 0) return(NULL)
  names(df)[hit[[1]]]
}

make_network_edges <- function(enrich_df) {
  empty <- data.frame(gene = character(), term = character(), adjusted_pvalue = numeric(), stringsAsFactors = FALSE)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(empty)
  if ("Message" %in% names(enrich_df)) return(empty)

  term_col <- find_first_col(enrich_df, c("term", "Term", "term_name", "name", "description"))
  genes_col <- find_first_col(enrich_df, c("genes", "Genes", "intersection", "overlapping genes", "overlap_genes", "core_enrichment"))
  padj_col <- find_first_col(enrich_df, c("adjusted_pvalue", "Adjusted.P.value", "Adjusted P-value", "adjusted p value", "padj", "p.adjust", "p_adjust", "qvalue", "fdr", "p_value", "pvalue"))

  if (is.null(term_col) || is.null(genes_col) || is.null(padj_col)) {
    stop("Cannot create network edges. Need term, genes/intersection, and adjusted p-value columns. Found: ", paste(names(enrich_df), collapse = ", "))
  }

  rows <- list()
  k <- 1L
  for (i in seq_len(nrow(enrich_df))) {
    term <- as.character(enrich_df[[term_col]][[i]])
    padj <- suppressWarnings(as.numeric(enrich_df[[padj_col]][[i]]))
    raw_genes <- enrich_df[[genes_col]][[i]]
    if (is.null(raw_genes) || length(raw_genes) == 0 || is.na(term) || term == "") next
    genes_text <- paste(as.character(raw_genes), collapse = ";")
    genes <- unlist(strsplit(genes_text, "[;,/|]+"))
    genes <- trimws(genes)
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) next
    for (gene in unique(genes)) {
      rows[[k]] <- data.frame(gene = gene, term = term, adjusted_pvalue = padj, stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }

  if (length(rows) == 0) return(empty)
  unique(do.call(rbind, rows))
}


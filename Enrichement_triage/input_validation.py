from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd

# Canonical fields used by triage.py. The triage layer still accepts aliases;
# this module detects the incoming schema earlier and gives user-friendly errors.
ALIASES = {
    "term": ["Term", "term", "pathway", "Pathway", "name", "Name", "Description", "Gene_set", "gene_set", "target", "label"],
    "genes": ["Genes", "genes", "gene", "Gene", "overlap_genes", "gene_list", "Genes;", "leadingEdge", "leading_edge", "core_enrichment", "source"],
    "padj": ["Adjusted.P.value", "Adjusted P-value", "Adjusted.P-val", "FDR", "qvalue", "q_value", "padj", "adj_p", "Adjusted P-value", "adjusted_pvalue", "adjusted_p_value", "best_adj_p", "weight_raw"],
    "pval": ["P.value", "pvalue", "P-value", "Pval", "p_val", "p.value", "p-value"],
    "overlap": ["Overlap", "overlap"],
    "odds_ratio": ["Odds.Ratio", "Odds Ratio", "odds_ratio", "OR"],
    "combined_score": ["Combined.Score", "Combined Score", "combined_score", "weight_plot", "priority_score", "consensus_score", "max_evidence", "followup_score"],
}

SCHEMA_HINTS = {
    "enrichr_like": {"required": ["term", "genes", "padj"], "nice": ["overlap", "pval", "odds_ratio", "combined_score"]},
    "clusterprofiler_like": {"required": ["term", "genes", "padj"], "nice": ["pval", "overlap"]},
    "gsea_like": {"required": ["term", "genes", "padj"], "nice": ["pval"]},
}


def _find_col(columns: List[str], candidates: List[str]) -> Optional[str]:
    exact = set(columns)
    for c in candidates:
        if c in exact:
            return c
    lower = {c.lower(): c for c in columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def detect_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    cols = list(df.columns)
    return {key: _find_col(cols, aliases) for key, aliases in ALIASES.items()}


def detect_schema(mapped: Dict[str, Optional[str]]) -> str:
    # Enrichr commonly has overlap + combined score + odds ratio.
    if mapped.get("overlap") and mapped.get("combined_score"):
        return "enrichr_like"
    # clusterProfiler commonly has core_enrichment or geneID mapped as genes.
    if mapped.get("genes") and mapped.get("pval") and mapped.get("padj"):
        return "clusterprofiler_or_gsea_like"
    return "unknown_or_minimal"


def validate_enrichment_df(df: pd.DataFrame, *, max_rows: int = 25000) -> Dict[str, Any]:
    if df is None or df.empty:
        return {"ok": False, "error": "The CSV is empty.", "columns": [], "mapped_columns": {}, "missing_required": ["term", "genes", "padj"]}

    mapped = detect_columns(df)
    required = ["term", "genes", "padj"]
    missing = [field for field in required if not mapped.get(field)]

    warnings: List[str] = []
    if len(df) > max_rows:
        warnings.append(f"Large file: {len(df)} rows. The app will run, but consider filtering to top enrichment results for faster demos.")
    if not mapped.get("overlap"):
        warnings.append("No Overlap column detected. The triage layer will infer overlap size from the gene list.")
    if not mapped.get("combined_score"):
        warnings.append("No combined score detected. Ranking will rely more heavily on adjusted p-value, overlap size, and BioFit rules.")

    ok = len(missing) == 0
    return {
        "ok": ok,
        "schema": detect_schema(mapped),
        "n_rows": int(len(df)),
        "columns": list(df.columns),
        "mapped_columns": mapped,
        "missing_required": missing,
        "warnings": warnings,
        "accepted_formats": [
            "Enrichr-like CSV",
            "clusterProfiler-like CSV",
            "GSEA-like CSV with term/gene/FDR columns",
            "Network Studio LLM triage bundle .zip",
            "Long-format gene/term/adjusted_pvalue CSV",
        ],
        "error": "" if ok else f"Missing required enrichment fields: {missing}. Required: term, genes, adjusted p-value/FDR.",
    }

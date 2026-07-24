from __future__ import annotations

from io import BytesIO
from typing import Any, Dict, Tuple
import os
import zipfile

import numpy as np
import pandas as pd


GENE_COLS = ["Genes", "genes", "gene", "Gene", "item", "Item", "source_gene"]
TERM_COLS = ["Term", "term", "pathway", "Pathway", "name", "Name", "Description", "Gene_set", "gene_set", "target_term"]
PADJ_COLS = ["Adjusted.P.value", "Adjusted P-value", "Adjusted.P-val", "FDR", "qvalue", "q_value", "padj", "adj_p", "adjusted_pvalue", "adjusted_p_value", "best_adj_p", "weight_raw"]
PVAL_COLS = ["P.value", "pvalue", "P-value", "Pval", "p_val", "p.value", "p-value"]
SCORE_COLS = ["Combined.Score", "Combined Score", "combined_score", "score", "priority_score", "consensus_score", "weight_plot", "max_evidence", "followup_score"]


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    exact = set(df.columns)
    for c in candidates:
        if c in exact:
            return c
    lower = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def _strip_network_prefix(value: Any) -> str:
    s = "" if value is None else str(value).strip()
    if "::" in s:
        return s.split("::", 1)[1].strip()
    return s


def _safe_numeric(series: pd.Series, default: float | None = None) -> pd.Series:
    out = pd.to_numeric(series, errors="coerce")
    if default is not None:
        out = out.fillna(default)
    return out


def _looks_like_long_gene_term_table(df: pd.DataFrame) -> bool:
    gene_col = _find_col(df, GENE_COLS)
    term_col = _find_col(df, TERM_COLS)
    genes_col = _find_col(df, ["Genes", "genes", "gene_list", "overlap_genes", "core_enrichment"])
    # If it has a singular gene column and a term column, aggregate it.
    # If it has a canonical Genes column, leave it as enrichment-like.
    return bool(gene_col and term_col and gene_col != genes_col)


def _aggregate_long_gene_term_table(df: pd.DataFrame, *, source_table: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    gene_col = _find_col(df, GENE_COLS)
    term_col = _find_col(df, TERM_COLS)
    padj_col = _find_col(df, PADJ_COLS)
    pval_col = _find_col(df, PVAL_COLS)
    score_col = _find_col(df, SCORE_COLS)

    if not gene_col or not term_col:
        raise ValueError(
            "Could not normalize long-format gene/term table. Expected columns like gene + term, "
            f"but found: {list(df.columns)}"
        )

    work = df.copy()
    work["__gene"] = work[gene_col].map(_strip_network_prefix)
    work["__term"] = work[term_col].map(_strip_network_prefix)
    work = work[(work["__gene"] != "") & (work["__term"] != "")]

    if work.empty:
        raise ValueError("The long-format gene/term table had no usable gene-term rows after cleaning.")

    if padj_col:
        work["__padj"] = _safe_numeric(work[padj_col], default=1.0)
    else:
        work["__padj"] = 1.0

    if pval_col:
        work["__pval"] = _safe_numeric(work[pval_col], default=np.nan)
    else:
        work["__pval"] = work["__padj"]

    if score_col:
        work["__score"] = _safe_numeric(work[score_col], default=np.nan)
    else:
        # Higher means stronger in the current triage scoring, so use -log10(adjusted p)
        work["__score"] = -np.log10(work["__padj"].clip(lower=1e-300))

    rows = []
    for term, grp in work.groupby("__term", sort=False):
        genes = []
        seen = set()
        for g in grp["__gene"].astype(str):
            if g and g not in seen:
                genes.append(g)
                seen.add(g)
        padj = float(grp["__padj"].min()) if len(grp) else 1.0
        pval = float(grp["__pval"].min()) if len(grp) else padj
        score = float(grp["__score"].max()) if len(grp) else 0.0
        rows.append({
            "Term": term,
            "Genes": ";".join(genes),
            "Adjusted.P.value": padj,
            "P.value": pval,
            "Overlap": f"{len(genes)}/{len(genes)}",
            "Combined.Score": score,
        })

    out = pd.DataFrame(rows)
    info = {
        "normalized": True,
        "source_table": source_table,
        "normalization": "long_gene_term_to_enrichment_like",
        "n_long_rows": int(len(work)),
        "n_terms_after_aggregation": int(len(out)),
        "columns_used": {
            "gene": gene_col,
            "term": term_col,
            "adjusted_pvalue": padj_col,
            "pvalue": pval_col,
            "score": score_col,
        },
    }
    return out, info


def _read_csv_from_filestorage(file_storage) -> pd.DataFrame:
    file_storage.stream.seek(0)
    return pd.read_csv(file_storage.stream)


def _read_bundle_zip(file_storage) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    file_storage.stream.seek(0)
    payload = file_storage.read()

    with zipfile.ZipFile(BytesIO(payload)) as zf:
        names = zf.namelist()
        csv_names = [n for n in names if n.lower().endswith(".csv")]

        # Prefer the bundle's explicit input preview if available because it already
        # expresses the gene/term/adjusted-pvalue table intended for LLM triage.
        preferred = [
            n for n in csv_names
            if os.path.basename(n).lower() in {"input_preview.csv", "main_edges.csv"}
        ]
        if not preferred:
            raise ValueError(
                "Bundle zip did not contain input_preview.csv or main_edges.csv. "
                f"CSV files found: {csv_names[:10]}"
            )

        # input_preview.csv first, main_edges.csv second.
        preferred.sort(key=lambda n: 0 if os.path.basename(n).lower() == "input_preview.csv" else 1)
        chosen = preferred[0]

        with zf.open(chosen) as fh:
            df = pd.read_csv(fh)

    if _looks_like_long_gene_term_table(df):
        out, info = _aggregate_long_gene_term_table(df, source_table=chosen)
    elif {"source", "target"}.issubset(set(df.columns)):
        # Network edge table fallback: source=item::GENE, target=group::TERM.
        edge_df = df.rename(columns={"source": "gene", "target": "term"}).copy()
        out, info = _aggregate_long_gene_term_table(edge_df, source_table=chosen)
    else:
        out = df
        info = {"normalized": False, "source_table": chosen, "normalization": "none"}

    info.update({
        "input_type": "llm_triage_bundle",
        "filename": getattr(file_storage, "filename", ""),
        "bundle_csv_files": csv_names,
    })
    return out, info


def load_analysis_input(file_storage) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Load either:
      1) a normal enrichment CSV, preserving the existing workflow, or
      2) a Network Studio LLM triage bundle zip, normalized into enrichment-like rows.

    Returns (df, input_info). df is always suitable for validate_enrichment_df() and
    run_enrichment_pipeline().
    """
    filename = (getattr(file_storage, "filename", "") or "").strip()
    lower = filename.lower()

    if lower.endswith(".zip"):
        return _read_bundle_zip(file_storage)

    if lower.endswith(".csv") or lower.endswith(".tsv") or lower.endswith(".txt") or not lower:
        df = _read_csv_from_filestorage(file_storage)
        info: Dict[str, Any] = {
            "input_type": "csv",
            "filename": filename,
            "normalized": False,
            "normalization": "none",
        }
        if _looks_like_long_gene_term_table(df):
            df, norm_info = _aggregate_long_gene_term_table(df, source_table=filename or "uploaded_csv")
            info.update(norm_info)
            info["input_type"] = "long_gene_term_csv"
        return df, info

    raise ValueError("Unsupported input file type. Please upload a .csv enrichment table or a .zip Network Studio bundle.")

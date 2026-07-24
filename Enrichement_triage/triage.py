# triage.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import math
import re


import numpy as np
import pandas as pd

from biofit import biofit_score, BioFitConfig

# -------------------------
# Column mapping / parsing
# -------------------------

_COL_ALIASES = {
    "term": ["Term", "term", "pathway", "Pathway", "name", "Name", "Description", "Gene_set", "gene_set", "target", "label"],
    "overlap": ["Overlap", "overlap"],
    "pval": ["P.value", "pvalue", "P-value", "Pval", "p_val", "p.value", "p-value"],
    "padj": ["Adjusted.P.value", "Adjusted P-value", "Adjusted.P-val", "FDR", "qvalue", "q_value", "padj", "adj_p", "adjusted_pvalue", "adjusted_p_value", "best_adj_p", "weight_raw"],
    "odds_ratio": ["Odds.Ratio", "Odds Ratio", "odds_ratio", "OR"],
    "combined_score": ["Combined.Score", "Combined Score", "combined_score", "weight_plot", "priority_score", "consensus_score", "max_evidence", "followup_score"],
    "genes": ["Genes", "genes", "gene", "Gene", "overlap_genes", "gene_list", "Genes;", "leadingEdge", "leading_edge", "core_enrichment", "source"],
}

_GENE_SPLIT_RE = re.compile(r"[;,]\s*")


def _pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    # tolerate case differences
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        lc = c.lower()
        if lc in lower_map:
            return lower_map[lc]
    return None


def _safe_float(x) -> Optional[float]:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    try:
        v = float(x)
        if np.isnan(v) or np.isinf(v):
            return None
        return v
    except Exception:
        return None


def parse_overlap(overlap: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Enrichr overlap typically looks like "3/171".
    Returns (k, n) = (overlap genes, gene set size).
    """
    if overlap is None or (isinstance(overlap, float) and np.isnan(overlap)):
        return None, None
    s = str(overlap).strip()
    if "/" not in s:
        return None, None
    left, right = s.split("/", 1)
    try:
        return int(left), int(right)
    except Exception:
        return None, None


def parse_genes(genes: str) -> List[str]:
    """
    Enrichr 'Genes' is usually "GENE1;GENE2;GENE3".
    """
    if genes is None or (isinstance(genes, float) and np.isnan(genes)):
        return []
    s = str(genes).strip()
    if not s:
        return []
    parts = [g.strip() for g in _GENE_SPLIT_RE.split(s) if g.strip()]
    # de-dupe while preserving order
    seen = set()
    out = []
    for g in parts:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


# -------------------------
# Triage scoring
# -------------------------

def _log1p_clip(x: Optional[float], lo: float = 0.0, hi: float = 10.0) -> float:
    if x is None:
        return 0.0
    return float(np.clip(math.log1p(max(x, 0.0)), lo, hi))


def _neglog10_p(p: Optional[float]) -> float:
    if p is None:
        return 0.0
    # floor to prevent inf
    p = max(p, 1e-300)
    return -math.log10(p)


def overlap_weight(k: int) -> float:
    """
    Penalize tiny overlaps hard, then saturate.
    """
    if k <= 0:
        return 0.0
    # 1 gene => 0.25, 2 => 0.45, 3 => 0.6, 5 => 0.8, 8+ => ~1.0
    return float(np.clip(1.0 - math.exp(-k / 3.0), 0.0, 1.0))


def triage_score(
    padj: Optional[float],
    odds_ratio: Optional[float],
    combined_score: Optional[float],
    overlap_k: int,
) -> float:
    """
    A stats-aware, biology-friendly pre-score.
    Not 'biological plausibility' yet — this is the filter that prevents junk from dominating
    before GPT-5 does phenotype grounding.
    """
    stat = _neglog10_p(padj)                       # backbone
    or_term = _log1p_clip(odds_ratio, hi=6.0)      # strong but capped
    cs_term = _log1p_clip(combined_score, hi=8.0)  # strong but capped
    ow = overlap_weight(overlap_k)

    # Weighting: padj dominates, OR and combined provide lift, overlap gates everything.
    raw = (1.2 * stat + 0.9 * or_term + 0.6 * cs_term) * (0.35 + 0.65 * ow)

    # Normalize to 0-100-ish for UI convenience
    return float(np.clip(raw * 8.0, 0.0, 100.0))


# -------------------------
# Redundancy clustering
# -------------------------

def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


@dataclass
class TriageConfig:
    # cluster terms if they share enough gene overlap
    jaccard_threshold: float = 0.35
    # maximum rows to cluster (keep fast)
    max_rows_for_clustering: int = 4000
    # keep only top N by triage for clustering (optional speed)
    cluster_top_n: int = 500


def cluster_by_gene_overlap(
    rows: List[Dict],
    cfg: TriageConfig,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Simple greedy clustering:
    - sort by combined_pre_gpt_score (fallback triage_score) desc
    - each term becomes a new cluster seed unless it matches an existing cluster by Jaccard
    - cluster membership uses gene sets
    Returns:
      clusters: [{cluster_id, seed_term, members, cluster_score, union_genes}]
      updated_rows: rows annotated with cluster_id
    """
    if len(rows) == 0:
        return [], rows

    score_key = lambda r: r.get("combined_pre_gpt_score", r["triage_score"])

    if len(rows) > cfg.max_rows_for_clustering:
        # avoid surprise quadratic runtime
        rows = sorted(rows, key=score_key, reverse=True)[: cfg.cluster_top_n]

    rows_sorted = sorted(rows, key=score_key, reverse=True)

    clusters = []
    for r in rows_sorted:
        gset = set(r["genes_list"])
        best_idx = None
        best_sim = 0.0

        for i, c in enumerate(clusters):
            sim = jaccard(gset, c["union_genes"])
            if sim > best_sim:
                best_sim = sim
                best_idx = i

        score_here = score_key(r)

        if best_idx is not None and best_sim >= cfg.jaccard_threshold:
            clusters[best_idx]["members"].append(r["row_id"])
            clusters[best_idx]["union_genes"] |= gset
            clusters[best_idx]["cluster_score"] = max(clusters[best_idx]["cluster_score"], score_here)

            r["cluster_id"] = clusters[best_idx]["cluster_id"]
            r["cluster_sim_to_seed"] = best_sim
        else:
            cid = f"C{len(clusters)+1:03d}"
            cluster = {
                "cluster_id": cid,
                "seed_term": r["term"],
                "members": [r["row_id"]],
                "cluster_score": score_here,
                "union_genes": set(gset),
            }
            clusters.append(cluster)
            r["cluster_id"] = cid
            r["cluster_sim_to_seed"] = 1.0

    # clean cluster output for JSON
    clusters_out = []
    for c in clusters:
        clusters_out.append({
            "cluster_id": c["cluster_id"],
            "seed_term": c["seed_term"],
            "member_count": len(c["members"]),
            "members": c["members"],
            "cluster_score": float(c["cluster_score"]),
            "union_genes": sorted(list(c["union_genes"])),
        })

    return clusters_out, rows_sorted


# -------------------------
# Main entrypoint
# -------------------------

def load_enrichr_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df

from typing import Dict, List, Optional
import pandas as pd

# add these imports at top of triage.py
# from biofit import biofit_score, BioFitConfig


def triage_enrichment_table(
    df: pd.DataFrame,
    cfg: Optional[TriageConfig] = None,
    *,
    phenotype: str = "",
    context: Optional[Dict] = None,
    bio_cfg: Optional["BioFitConfig"] = None,
) -> Dict:
    cfg = cfg or TriageConfig()
    context = context or {}

    # find columns
    col_term = _pick_col(df, _COL_ALIASES["term"])
    col_overlap = _pick_col(df, _COL_ALIASES["overlap"])
    col_padj = _pick_col(df, _COL_ALIASES["padj"])
    col_or = _pick_col(df, _COL_ALIASES["odds_ratio"])
    col_cs = _pick_col(df, _COL_ALIASES["combined_score"])
    col_genes = _pick_col(df, _COL_ALIASES["genes"])
    col_p = _pick_col(df, _COL_ALIASES["pval"])

    missing = [k for k, v in [
        ("term", col_term),
        ("genes", col_genes),
        ("padj", col_padj),
    ] if v is None]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {list(df.columns)}")

    has_pheno = bool((phenotype or "").strip())
    # If phenotype exists, BioFit should matter more.
    w_bio = 0.60 if has_pheno else 0.35
    w_tri = 1.0 - w_bio

    rows: List[Dict] = []
    for i, rec in df.iterrows():
        term = str(rec[col_term]).strip()
        genes_list = parse_genes(rec[col_genes])

        k, n = parse_overlap(rec[col_overlap]) if col_overlap else (len(genes_list), None)

        padj = _safe_float(rec[col_padj]) if col_padj else None
        pval = _safe_float(rec[col_p]) if col_p else None
        odds_ratio = _safe_float(rec[col_or]) if col_or else None
        combined_score = _safe_float(rec[col_cs]) if col_cs else None

        overlap_k = k if k is not None else len(genes_list)

        # 1) Stats-aware triage score (0-100)
        tri_score = triage_score(
            padj=padj,
            odds_ratio=odds_ratio,
            combined_score=combined_score,
            overlap_k=overlap_k,
        )

        # 2) Deterministic biological plausibility (0-100)
        bf = biofit_score(
            term=term,
            genes=genes_list,
            overlap_k=overlap_k,
            phenotype=phenotype or "",
            context=context,
            cfg=bio_cfg,
        )

        # 3) Combined "pre-GPT" ranking score
        combined_pre_gpt = float(max(0.0, min(100.0, w_tri * tri_score + w_bio * bf["biofit_score"])))

        # Lightweight flags (fast, interpretable)
        flags = []
        if overlap_k <= 2:
            flags.append("tiny_overlap")
        if padj is not None and padj > 0.2:
            flags.append("weak_stats")
        if any(x in term.lower() for x in ["metabolic process", "regulation of", "cellular process"]) and overlap_k <= 2:
            flags.append("generic_term_small_overlap")

        # include biofit flags too
        flags.extend(bf.get("flags", []))

        rows.append({
            "row_id": f"R{i+1:04d}",
            "term": term,
            "p_value": pval,
            "adjusted_p_value": padj,
            "odds_ratio": odds_ratio,
            "combined_score": combined_score,
            "overlap_k": overlap_k,
            "overlap_n": n,
            "genes_list": genes_list,

            "triage_score": float(tri_score),

            # BioFit outputs (deterministic biology checks)
            "biofit_score": float(bf["biofit_score"]),
            "biofit_program": bf.get("best_program"),
            "biofit_components": bf.get("components", {}),
            "biofit_artifact_reasons": bf.get("artifact_reasons", []),

            # final sort key before GPT-5
            "combined_pre_gpt_score": combined_pre_gpt,

            "flags": flags,
        })

    # Sort by combined score (this is what you’ll feed into clustering + GPT)
    rows_sorted = sorted(rows, key=lambda r: r["combined_pre_gpt_score"], reverse=True)

    # Redundancy clustering:
    # Option A (recommended): cluster using combined score
    clusters, rows_clustered = cluster_by_gene_overlap(
        rows_sorted, cfg
    )
    # If you want clustering to also rank by combined score, change cluster_by_gene_overlap
    # to reference r["combined_pre_gpt_score"] instead of r["triage_score"].

    return {
        "rows": rows_clustered,
        "clusters": clusters,
        "meta": {
            "n_rows": len(rows),
            "n_clusters": len(clusters),
            "has_phenotype": has_pheno,
            "weights": {"w_triage": w_tri, "w_biofit": w_bio},
            "config": {
                "jaccard_threshold": cfg.jaccard_threshold,
                "max_rows_for_clustering": cfg.max_rows_for_clustering,
                "cluster_top_n": cfg.cluster_top_n,
            },
        },
    }


    # basic sorting
    rows_sorted = sorted(rows, key=lambda r: r.get("combined_pre_gpt_score", r["triage_score"]), reverse=True)

    # redundancy clustering
    clusters, rows_clustered = cluster_by_gene_overlap(rows_sorted, cfg)

    # lightweight flags (useful to show before GPT-5)
    for r in rows_clustered:
        flags = []
        if r["overlap_k"] <= 2:
            flags.append("tiny_overlap")
        if r["adjusted_p_value"] is not None and r["adjusted_p_value"] > 0.2:
            flags.append("weak_stats")
        # very generic-looking GO terms (cheap heuristic; GPT-5 will do the real job later)
        if any(x in r["term"].lower() for x in ["metabolic process", "regulation of", "cellular process"]) and r["overlap_k"] <= 2:
            flags.append("generic_term_small_overlap")
        r["flags"] = flags

    return {
        "rows": rows_clustered,
        "clusters": clusters,
        "meta": {
            "n_rows": len(rows),
            "n_clusters": len(clusters),
            "config": {
                "jaccard_threshold": cfg.jaccard_threshold,
                "max_rows_for_clustering": cfg.max_rows_for_clustering,
                "cluster_top_n": cfg.cluster_top_n,
            },
        },
    }

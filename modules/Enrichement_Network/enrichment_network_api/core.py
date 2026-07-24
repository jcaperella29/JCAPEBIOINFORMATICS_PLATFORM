from __future__ import annotations
import math
from typing import Dict, Tuple
from typing import Dict
import networkx as nx
import numpy as np
import base64
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import io
import json
import zipfile
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

from io import StringIO



# ----------------------------
# Helpers
# ----------------------------
REQUIRED_RAW_COLS_NOTE = (
    "Raw CSV should be long-format: one row per membership edge (gene ↔ pathway). "
    "You’ll map the item + group columns after upload."
)

import plotly.express as px
import math


def edge_weight_from_adj_p(p: float) -> float:
    """Convert adjusted p-value to a 'bigger = more significant' weight."""
    if p is None:
        return 0.0
    try:
        p = float(p)
    except Exception:
        return 0.0
    # clamp to avoid log(0)
    p = max(p, 1e-300)
    return max(0.0, -math.log10(p))


def minmax_scale(x: float, xmin: float, xmax: float, out_min: float, out_max: float) -> float:
    """Map x from [xmin,xmax] to [out_min,out_max]."""
    if xmax == xmin:
        return (out_min + out_max) / 2.0
    return out_min + (x - xmin) * (out_max - out_min) / (xmax - xmin)


def term_color_map(terms):
    palette = px.colors.qualitative.Set3 + px.colors.qualitative.Dark24
    return {t: palette[i % len(palette)] for i, t in enumerate(sorted(terms))}


def parse_upload(contents: str) -> pd.DataFrame:
    """Parse Dash upload contents into a pandas DataFrame (CSV)."""
    content_type, content_string = contents.split(",", 1)
    decoded = base64.b64decode(content_string)
    return pd.read_csv(io.StringIO(decoded.decode("utf-8")))

DEMO_CSV_PATH = Path(__file__).parent / "TEST.csv"


def load_demo_dataframe() -> pd.DataFrame:
    """
    Load the bundled demo file and convert Enrichr-style term -> semicolon genes
    into long-format gene -> term rows for the network builder.
    """
    df = pd.read_csv(DEMO_CSV_PATH)

    required = {"term", "genes", "adjusted_pvalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Demo CSV is missing required columns: {sorted(missing)}")

    long_df = (
        df.assign(gene=df["genes"].astype(str).str.split(";"))
          .explode("gene")
          .assign(gene=lambda x: x["gene"].astype(str).str.strip())
    )

    long_df = long_df[long_df["gene"].ne("")]
    long_df = long_df[["gene", "term", "adjusted_pvalue"]].drop_duplicates()

    return long_df


COLUMN_MAPPING_PRESETS = [
    {"label": "Custom long-format CSV", "value": "custom_long"},
    {"label": "Enrichr-style: term + semicolon genes", "value": "enrichr"},
    {"label": "g:Profiler / gprofiler2-style", "value": "gprofiler"},
    {"label": "clusterProfiler-style", "value": "clusterprofiler"},
    {"label": "GSEA / MSigDB-style", "value": "gsea_msigdb"},
]


def _normalize_col_name(name: str) -> str:
    """Normalize column names for forgiving preset matching."""
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Find a dataframe column from likely names, ignoring case/punctuation."""
    normalized = {_normalize_col_name(c): c for c in df.columns}
    for candidate in candidates:
        hit = normalized.get(_normalize_col_name(candidate))
        if hit is not None:
            return hit
    return None


def _guess_long_format_columns(df: pd.DataFrame) -> tuple[str | None, str | None, str | None]:
    """Guess item/group/weight columns for already-long enrichment edge tables."""
    item_col = _find_col(
        df,
        [
            "gene", "genes", "item", "items", "symbol", "gene_symbol", "geneid",
            "gene_id", "target", "feature", "node", "protein",
        ],
    )
    group_col = _find_col(
        df,
        [
            "term", "pathway", "pathways", "group", "description", "name",
            "term_name", "termid", "term_id", "geneset", "gene_set",
        ],
    )
    weight_col = _find_col(
        df,
        [
            "adjusted_pvalue", "adjusted_p_value", "adjusted p-value", "adjusted p value",
            "adjusted.p.value", "padj", "p.adjust", "p_adjust", "qvalue", "q_value",
            "fdr", "fdr q-val", "fdr_q_val", "p_value", "pvalue", "p.val", "pval",
        ],
    )
    return item_col, group_col, weight_col


def _split_gene_memberships(value) -> list[str]:
    """Split common enrichment gene-list cells into individual gene symbols/items."""
    import re

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []

    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return []

    # clusterProfiler often uses '/', Enrichr often ';', many CSV exports use ','.
    # Also handle pipes and whitespace as secondary separators.
    parts = re.split(r"[;,/|]+", text)
    genes = [p.strip() for p in parts if p and p.strip()]
    return genes


def _expand_membership_table(
    df: pd.DataFrame,
    term_col: str,
    genes_col: str,
    weight_col: str | None,
    preset_label: str,
) -> tuple[pd.DataFrame, str, str, str | None, str]:
    """Convert term -> gene-list tables into long gene -> term membership rows."""
    records = []

    for _, row in df.iterrows():
        term = row.get(term_col)
        genes = _split_gene_memberships(row.get(genes_col))
        if term is None or pd.isna(term) or not str(term).strip():
            continue
        for gene in genes:
            rec = {"gene": gene, "term": str(term).strip()}
            if weight_col:
                rec["adjusted_pvalue"] = row.get(weight_col)
            records.append(rec)

    if not records:
        raise ValueError(
            f"{preset_label} preset found columns, but no gene memberships could be expanded. "
            "Check the gene-list separator or choose Custom long-format CSV."
        )

    long_df = pd.DataFrame(records).drop_duplicates()
    item_col = "gene"
    group_col = "term"
    out_weight_col = "adjusted_pvalue" if "adjusted_pvalue" in long_df.columns else None

    msg = (
        f"Applied {preset_label} preset: expanded {len(long_df):,} gene-term edges "
        f"from {long_df['gene'].nunique():,} genes/items and {long_df['term'].nunique():,} pathways/terms."
    )
    return long_df, item_col, group_col, out_weight_col, msg


def apply_column_mapping_preset(df: pd.DataFrame, preset: str) -> tuple[pd.DataFrame, str | None, str | None, str | None, str]:
    """
    Apply a column mapping preset.

    For term -> gene-list formats, this converts the table into the long-format
    gene ↔ pathway edge table used by the network builder.
    """
    preset = preset or "custom_long"

    if preset == "custom_long":
        item_col, group_col, weight_col = _guess_long_format_columns(df)
        msg = "Applied Custom long-format preset."
        if item_col and group_col:
            msg += f" Guessed item={item_col}, group={group_col}"
            if weight_col:
                msg += f", weight={weight_col}"
            msg += "."
        else:
            msg += " Select item and group columns manually."
        return df, item_col, group_col, weight_col, msg

    if preset == "enrichr":
        term_col = _find_col(df, ["term", "Term", "name", "pathway", "description"])
        genes_col = _find_col(df, ["genes", "Genes", "overlapping genes", "overlap_genes"])
        weight_col = _find_col(df, ["adjusted_pvalue", "Adjusted.P.value", "Adjusted P-value", "adjusted p value", "padj", "fdr", "pvalue", "p_value"])
        if not term_col or not genes_col:
            raise ValueError("Enrichr preset needs a term column and a genes column.")
        return _expand_membership_table(df, term_col, genes_col, weight_col, "Enrichr-style")

    if preset == "gprofiler":
        term_col = _find_col(df, ["term_name", "name", "term", "native", "term_id", "description"])
        genes_col = _find_col(df, ["intersection", "intersections", "genes", "query", "members"])
        weight_col = _find_col(df, ["p_value", "pvalue", "adjusted_pvalue", "padj", "fdr", "qvalue"])
        if not term_col or not genes_col:
            raise ValueError("g:Profiler preset needs a term/name column and an intersection/genes column.")
        return _expand_membership_table(df, term_col, genes_col, weight_col, "g:Profiler/gprofiler2-style")

    if preset == "clusterprofiler":
        term_col = _find_col(df, ["Description", "description", "term", "pathway", "ID", "id"])
        genes_col = _find_col(df, ["geneID", "gene_id", "genes", "Genes", "core_enrichment"])
        weight_col = _find_col(df, ["p.adjust", "p_adjust", "padj", "qvalue", "pvalue", "p_value"])
        if not term_col or not genes_col:
            raise ValueError("clusterProfiler preset needs Description/ID and geneID/core_enrichment columns.")
        return _expand_membership_table(df, term_col, genes_col, weight_col, "clusterProfiler-style")

    if preset == "gsea_msigdb":
        term_col = _find_col(df, ["NAME", "name", "Description", "description", "term", "pathway", "geneset", "gene_set"])
        genes_col = _find_col(df, ["core_enrichment", "leading_edge", "leading edge", "genes", "Genes", "members"])
        weight_col = _find_col(df, ["FDR q-val", "FDR.q.val", "fdr", "qvalue", "q_value", "p.adjust", "padj", "pvalue", "p_value"])
        if not term_col or not genes_col:
            raise ValueError("GSEA/MSigDB preset needs a NAME/term column and core_enrichment/genes column.")
        return _expand_membership_table(df, term_col, genes_col, weight_col, "GSEA/MSigDB-style")

    raise ValueError(f"Unknown column mapping preset: {preset}")


import math
import pandas as pd
import networkx as nx


def build_bipartite_graph(
    df: pd.DataFrame,
    item_col: str,
    group_col: str,
    weight_col: str | None = None,
) -> nx.Graph:
    """
    Build an undirected bipartite graph from long-form rows.

    If weight_col is provided (e.g. adjusted p-value),
    edge weight is transformed to:
        weight = -log10(adjusted_pvalue)

    This means:
      - smaller p-values → stronger edges
      - filtering with weight >= threshold behaves correctly
    """
    g = nx.Graph()

    cols = [item_col, group_col] + ([weight_col] if weight_col else [])
    tmp = df[cols].copy()

    # Clean strings
    tmp[item_col] = tmp[item_col].astype(str).str.strip()
    tmp[group_col] = tmp[group_col].astype(str).str.strip()

    # Drop junk
    tmp = tmp.dropna(subset=[item_col, group_col])
    tmp = tmp[(tmp[item_col] != "") & (tmp[group_col] != "")]
    tmp = tmp.drop_duplicates(subset=[item_col, group_col])

    # Add nodes
    items = tmp[item_col].unique().tolist()
    groups = tmp[group_col].unique().tolist()
    for it in items:
        g.add_node(f"item::{it}", label=it, node_type="item")

    for gr in groups:
        g.add_node(f"group::{gr}", label=gr, node_type="group")

    if weight_col and weight_col in tmp.columns:
        w_raw = pd.to_numeric(tmp[weight_col], errors="coerce").fillna(1.0)

        for it, gr, raw in zip(tmp[item_col], tmp[group_col], w_raw):
            raw = float(raw)
            w_plot = edge_weight_from_adj_p(raw)

            g.add_edge(
                f"item::{it}",
                f"group::{gr}",
                weight_raw=raw,
                weight_plot=w_plot
            )
    else:
        for it, gr in zip(tmp[item_col], tmp[group_col]):
            g.add_edge(
                f"item::{it}",
                f"group::{gr}",
                weight_raw=1.0,
                weight_plot=1.0
            )

    return g


def layout_bipartite_two_column(g: nx.Graph) -> Dict[str, Tuple[float, float]]:
    """
    Deterministic bipartite layout:
    - groups on left (x=0), items on right (x=1)
    - y sorted by degree so hubs sit near center
    """
    groups = [n for n, d in g.nodes(data=True) if d.get("node_type") == "group"]
    items = [n for n, d in g.nodes(data=True) if d.get("node_type") == "item"]

    groups_sorted = sorted(groups, key=lambda n: g.degree(n), reverse=True)
    items_sorted = sorted(items, key=lambda n: g.degree(n), reverse=True)

    def y_positions(nodes_sorted):
        if len(nodes_sorted) == 1:
            return {nodes_sorted[0]: 0.0}
        step = 2.0 / (len(nodes_sorted) - 1)
        return {n: 1.0 - i * step for i, n in enumerate(nodes_sorted)}

    yg = y_positions(groups_sorted)
    yi = y_positions(items_sorted)

    pos = {}
    for n in groups_sorted:
        pos[n] = (0.0, yg[n])
    for n in items_sorted:
        pos[n] = (1.0, yi[n])

    return pos


def subgraph_filter(
    g: nx.Graph,
    search: str,
    min_degree: int,
    min_weight: float,
    max_groups: int,
    largest_component_only: bool,
) -> nx.Graph:
    """Return a filtered subgraph based on UI controls."""

    # ---- Helper: choose the right edge weight key (prefer weight_plot) ----
    def _edge_w(ed: dict) -> float:
        for k in ("weight_plot", "weight_raw", "weight"):
            if k in ed and ed[k] is not None:
                try:
                    return float(ed[k])
                except Exception:
                    pass
        return 1.0

    # Edge weight filter
    edges_keep = [(u, v) for u, v, ed in g.edges(data=True) if _edge_w(ed) >= min_weight]
    sg = g.edge_subgraph(edges_keep).copy()

    # Degree filter
    nodes_keep = [n for n in sg.nodes() if sg.degree(n) >= min_degree]
    sg = sg.subgraph(nodes_keep).copy()

    # Search (keep matches + their neighbors)
    s = (search or "").strip().lower()
    if s:
        hits = [n for n, d in sg.nodes(data=True) if s in str(d.get("label", "")).lower()]
        expanded = set(hits)
        for n in hits:
            expanded.update(list(sg.neighbors(n)))
        sg = sg.subgraph(list(expanded)).copy()

    # Limit number of groups shown
    if max_groups and max_groups > 0:
        groups = [(n, sg.degree(n)) for n, d in sg.nodes(data=True) if d.get("node_type") == "group"]
        groups_sorted = sorted(groups, key=lambda x: x[1], reverse=True)
        allowed_groups = set([n for n, _ in groups_sorted[:max_groups]])

        if allowed_groups:
            keep = []
            for n, d in sg.nodes(data=True):
                if d.get("node_type") == "group":
                    if n in allowed_groups:
                        keep.append(n)
                else:
                    if any(nei in allowed_groups for nei in sg.neighbors(n)):
                        keep.append(n)
            sg = sg.subgraph(keep).copy()

    # Largest connected component only
    if largest_component_only and sg.number_of_nodes() > 0:
        comps = list(nx.connected_components(sg))
        if comps:
            biggest = max(comps, key=len)
            sg = sg.subgraph(list(biggest)).copy()

    return sg


import math
from typing import Dict, Tuple
import networkx as nx
import plotly.graph_objects as go

import math
from typing import Dict, Tuple
import networkx as nx
import plotly.graph_objects as go


def make_plotly_network(
    g: nx.Graph,
    pos: Dict[str, Tuple[float, float]],
    show_labels: bool,
    thickness_by_weight: bool = False,
    edge_width_range: Tuple[float, float] = (1.5, 6.0),
    edge_weight_range: Tuple[float, float] = (0.0, 10.0),
    highlight_nodes: dict | None = None,
) -> go.Figure:
    """
    Draw a bipartite-ish network with:
      - edges colored by term (group node)
      - optional edge thickness scaling by weight
      - optional candidate highlighting from the Insights tab
      - automatic conversion of p-values to -log10(p) when needed

    Weight handling:
      Prefer edge attr: weight_plot -> weight -> weight_raw
      If the chosen value looks like a p-value (0<val<=1), convert to -log10(val).
    """

    # ----------------------------
    # Helper: pick best available weight
    # ----------------------------
    def get_edge_weight(ed: dict) -> float:
        """
        Prefer weight_plot (already -log10(padj) ideally).
        Otherwise fall back to weight, then weight_raw.
        If it looks like a p-value (0<val<=1), convert to -log10(val).
        """
        if "weight_plot" in ed and ed["weight_plot"] is not None:
            val = ed["weight_plot"]
        elif "weight" in ed and ed["weight"] is not None:
            val = ed["weight"]
        elif "weight_raw" in ed and ed["weight_raw"] is not None:
            val = ed["weight_raw"]
        else:
            return 1.0

        try:
            val = float(val)
        except Exception:
            return 1.0

        # Auto-convert p-values to -log10(p)
        if 0.0 < val <= 1.0:
            val = max(val, 1e-300)
            return max(0.0, -math.log10(val))

        return float(val)

    # ----------------------------
    # Group edges by term (group node), keep weights
    # ----------------------------
    term_edges: dict = {}
    all_w: list[float] = []

    for u, v, ed in g.edges(data=True):
        if g.nodes[u].get("node_type") == "group":
            term = u
            gene = v
        else:
            term = v
            gene = u

        w = get_edge_weight(ed)
        term_edges.setdefault(term, []).append((term, gene, w))
        all_w.append(w)

    term_colors = term_color_map(term_edges.keys())

    # ----------------------------
    # Thickness scaling (percentile-based so small differences still show)
    # ----------------------------
    min_px, max_px = edge_width_range

    # User slider range (kept, but we also do robust scaling from the data)
    min_w_user, max_w_user = edge_weight_range
    user_has_valid_range = (max_w_user > min_w_user)

    if thickness_by_weight and all_w:
        w_arr = np.asarray(all_w, dtype=float)

        # Robust range: percentiles prevent 1 extreme edge from dominating
        lo_data, hi_data = np.percentile(w_arr, [5, 95])
        if hi_data <= lo_data:
            lo_data, hi_data = float(w_arr.min()), float(w_arr.max())

        # Optional: incorporate user range if they gave a sane one
        # (We intersect user range with data range to keep behavior intuitive.)
        if user_has_valid_range:
            lo = max(float(lo_data), float(min_w_user))
            hi = min(float(hi_data), float(max_w_user))
            if hi <= lo:
                lo, hi = float(lo_data), float(hi_data)
        else:
            lo, hi = float(lo_data), float(hi_data)

        def width_from_weight(weight: float) -> float:
            w = float(weight)

            # Clamp to robust range
            if hi > lo:
                w = max(lo, min(hi, w))
                t = (w - lo) / (hi - lo)
            else:
                t = 0.5

            # Contrast boost:
            # gamma < 1 makes small differences *more* visible (what you want for p-values)
            gamma = 0.55
            t = t ** gamma

            return float(min_px + t * (max_px - min_px))

    else:
        def width_from_weight(weight: float) -> float:
            return float(min_px)

    # ----------------------------
    # Edges
    # ----------------------------
    edge_traces = []

    for term, edges in term_edges.items():
        if thickness_by_weight:
            # One trace per edge (Plotly can't vary width inside a single trace)
            for tnode, gnode, w in edges:
                x0, y0 = pos[tnode]
                x1, y1 = pos[gnode]
                edge_traces.append(
                    go.Scatter(
                        x=[x0, x1, None],
                        y=[y0, y1, None],
                        mode="lines",
                        hoverinfo="none",
                        line=dict(
                            width=width_from_weight(w),
                            color=term_colors[term],
                        ),
                        opacity=0.65,
                        showlegend=False,
                    )
                )
        else:
            # Fast mode: one trace per term, constant width
            ex, ey = [], []
            for tnode, gnode, _w in edges:
                x0, y0 = pos[tnode]
                x1, y1 = pos[gnode]
                ex += [x0, x1, None]
                ey += [y0, y1, None]

            edge_traces.append(
                go.Scatter(
                    x=ex,
                    y=ey,
                    mode="lines",
                    hoverinfo="none",
                    line=dict(
                        width=float(min_px),
                        color=term_colors[term],
                    ),
                    opacity=0.65,
                    showlegend=False,
                )
            )

    # ----------------------------
    # Nodes
    # ----------------------------
    highlight_nodes = highlight_nodes or {}
    node_x, node_y, node_text, node_hover, node_size, node_color = [], [], [], [], [], []

    for n, d in g.nodes(data=True):
        x, y = pos[n]
        node_x.append(x)
        node_y.append(y)

        label = str(d.get("label", n))
        ntype = d.get("node_type", "unknown")
        deg = g.degree(n)

        if n in highlight_nodes:
            h = highlight_nodes[n]
            rank = h.get("candidate_rank", "")
            score = h.get("followup_score", "")
            node_hover.append(
                f"⭐ Candidate #{rank}: {label}"
                f"<br>type={ntype}"
                f"<br>degree={deg}"
                f"<br>follow-up score={score}"
            )
            node_text.append(f"★ {rank}. {label}")
            node_size.append(22 + max(0, 12 - int(rank or 12)))
            node_color.append("rgba(255,193,7,0.98)")
        else:
            node_hover.append(f"{label}<br>type={ntype}<br>degree={deg}")
            node_text.append(label if show_labels else "")
            node_size.append(8 + min(24, deg * 2))

            if ntype == "group":
                node_color.append("rgba(0,119,182,0.95)")
            else:
                node_color.append("rgba(0,180,216,0.95)")

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text" if show_labels else "markers",
        text=node_text,
        textposition="top center",
        hovertext=node_hover,
        hoverinfo="text",
        marker=dict(
            size=node_size,
            color=node_color,
            line=dict(
                width=[3 if n in highlight_nodes else 1 for n in g.nodes()],
                color=["rgba(120,53,15,0.95)" if n in highlight_nodes else "rgba(0,0,0,0.25)" for n in g.nodes()],
            ),
        ),
        name="nodes",
    )

    fig = go.Figure(data=edge_traces + [node_trace])
    fig.update_layout(
        margin=dict(l=10, r=10, t=10, b=10),
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        dragmode="pan",
        hovermode="closest",
    )

    return fig


def graph_stats(g: nx.Graph) -> Dict[str, object]:
    if g.number_of_nodes() == 0:
        return dict(
            nodes=0, edges=0, items=0, groups=0, components=0, largest_component=0,
            top_groups=[], top_items=[],
            weight_key=None,
            weight_summary=None,
            top_weighted_edges=[]
        )

    items = [n for n, d in g.nodes(data=True) if d.get("node_type") == "item"]
    groups = [n for n, d in g.nodes(data=True) if d.get("node_type") == "group"]

    comps = list(nx.connected_components(g))
    largest = max((len(c) for c in comps), default=0)

    top_groups = sorted(
        [(g.nodes[n].get("label", n), g.degree(n)) for n in groups],
        key=lambda x: x[1], reverse=True
    )[:20]

    top_items = sorted(
        [(g.nodes[n].get("label", n), g.degree(n)) for n in items],
        key=lambda x: x[1], reverse=True
    )[:20]

    # ----------------------------
    # Edge weight stats (if present)
    # Prefer weight_plot > weight > weight_raw
    # ----------------------------
    def pick_weight(ed: dict):
        for k in ("weight_plot", "weight", "weight_raw"):
            if k in ed and ed[k] is not None:
                try:
                    return k, float(ed[k])
                except Exception:
                    pass
        return None, None

    weights = []
    chosen_key = None

    for u, v, ed in g.edges(data=True):
        k, w = pick_weight(ed)
        if w is None:
            continue
        if chosen_key is None:
            chosen_key = k
        # If some edges have multiple keys, keep using the preferred order:
        # Only accept values that match the chosen key.
        if k == chosen_key:
            weights.append(w)

    weight_summary = None
    top_weighted_edges = []

    if weights:
        arr = np.array(weights, dtype=float)

        weight_summary = dict(
            key=chosen_key,
            min=float(np.min(arr)),
            p25=float(np.percentile(arr, 25)),
            median=float(np.median(arr)),
            mean=float(np.mean(arr)),
            p75=float(np.percentile(arr, 75)),
            max=float(np.max(arr)),
        )

        # Top weighted edges (show endpoints + labels + weight)
        # Note: only include edges that actually have the chosen_key
        edges_with_w = []
        for u, v, ed in g.edges(data=True):
            if chosen_key in ed and ed[chosen_key] is not None:
                try:
                    w = float(ed[chosen_key])
                except Exception:
                    continue
                u_lab = g.nodes[u].get("label", u)
                v_lab = g.nodes[v].get("label", v)
                edges_with_w.append((u_lab, v_lab, w))

        edges_with_w.sort(key=lambda x: x[2], reverse=True)
        top_weighted_edges = edges_with_w[:20]

    return dict(
        nodes=g.number_of_nodes(),
        edges=g.number_of_edges(),
        items=len(items),
        groups=len(groups),
        components=len(comps),
        largest_component=largest,
        top_groups=top_groups,
        top_items=top_items,
        weight_key=chosen_key,
        weight_summary=weight_summary,
        top_weighted_edges=top_weighted_edges,
    )


# ----------------------------
# Large graph warning helpers
# ----------------------------
LARGE_GRAPH_NODE_WARNING = 1500
LARGE_GRAPH_EDGE_WARNING = 5000
VERY_LARGE_RENDER_NODE_LIMIT = 1200
VERY_LARGE_RENDER_EDGE_LIMIT = 4000


def graph_size_warning_component(node_count: int, edge_count: int, context: str = "network"):
    """Return a serializable warning dict for large graphs, or an empty dict for normal sizes."""
    if node_count < LARGE_GRAPH_NODE_WARNING and edge_count < LARGE_GRAPH_EDGE_WARNING:
        return {}
    return {
        "level": "warning",
        "title": "Large graph warning",
        "message": (
            f"The current {context} has {node_count:,} nodes and {edge_count:,} edges. Rendering may be slow. "
            "Increase minimum node degree, increase minimum edge weight, reduce maximum groups, "
            "use search to focus on a pathway family, or avoid force-directed layout."
        ),
    }


def large_graph_placeholder_figure(node_count: int, edge_count: int) -> go.Figure:
    """Avoid trying to render extremely large filtered graphs in the browser."""
    fig = go.Figure()
    fig.update_layout(
        annotations=[
            dict(
                text=(
                    "Large filtered graph skipped for browser performance.<br>"
                    f"Filtered graph: {node_count:,} nodes, {edge_count:,} edges.<br>"
                    "Use min degree, min edge weight, max groups, or search filters to reduce the graph."
                ),
                showarrow=False,
                x=0.5,
                y=0.5,
                xref="paper",
                yref="paper",
                font=dict(size=16),
                align="center",
            )
        ],
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    return fig



def rebuild_graph_from_store(graph_data: dict) -> nx.Graph:
    """Rebuild a NetworkX graph from the serialized Dash store."""
    g = nx.Graph()

    if not graph_data:
        return g

    for nd in graph_data.get("nodes", []):
        node_id = nd.get("id")
        if not node_id:
            continue
        attrs = nd.copy()
        attrs.pop("id", None)
        g.add_node(node_id, **attrs)

    for ed in graph_data.get("edges", []):
        u, v = ed.get("source"), ed.get("target")
        if not u or not v:
            continue
        attrs = ed.copy()
        attrs.pop("source", None)
        attrs.pop("target", None)
        g.add_edge(u, v, **attrs)

    return g


def _normalized_edge_weight(ed: dict) -> float:
    """
    Return a positive edge weight suitable for random-walk analysis.

    Prefer weight_plot because this app already stores adjusted p-values as
    -log10(adjusted p-value). Fall back to weight, then weight_raw.
    If a fallback looks like a p-value, convert it to -log10(p).
    """
    val = None
    for key in ("weight_plot", "weight", "weight_raw"):
        if key in ed and ed[key] is not None:
            val = ed[key]
            break

    if val is None:
        return 1.0

    try:
        val = float(val)
    except Exception:
        return 1.0

    if 0.0 < val <= 1.0:
        val = max(val, 1e-300)
        val = -math.log10(val)

    # PageRank needs non-negative weights. A tiny floor avoids all-zero issues.
    return max(float(val), 1e-12)


def prepare_markov_graph(g: nx.Graph, ranking_mode: str = "balanced") -> nx.Graph:
    """
    Copy the current graph and add a 'markov_weight' edge attribute.

    This does not alter the visualization graph. It only prepares a separate
    weighted graph for signal diffusion / random-walk analysis.

    Ranking modes:
      - evidence: edge weights are driven by -log10(adjusted p-value)
      - connectivity: edge weights are treated equally so topology dominates
      - balanced: dampens extreme p-value effects while preserving evidence
    """
    mode = (ranking_mode or "balanced").lower().strip()
    mg = g.copy()

    for u, v, ed in mg.edges(data=True):
        evidence_w = _normalized_edge_weight(ed)

        if mode == "evidence":
            markov_w = evidence_w
        elif mode == "connectivity":
            markov_w = 1.0
        else:
            # Balanced mode keeps statistical evidence but compresses extreme
            # single-edge p-value effects so connected pathway structure matters.
            markov_w = math.sqrt(evidence_w)

        ed["markov_weight"] = max(float(markov_w), 1e-12)

    return mg


def _degree_adjusted_score(g: nx.Graph, node_id: str, pagerank_score: float, ranking_mode: str) -> float:
    """
    Convert the raw random-walk score into a user-facing priority score.

    Raw PageRank can still reward a strong isolated edge. For enrichment
    interpretation, especially pathway prioritization, degree should matter:
    pathways supported by multiple genes/items should not be buried by a
    one-edge term with an extremely strong p-value.
    """
    mode = (ranking_mode or "balanced").lower().strip()
    degree = max(int(g.degree(node_id)), 1)

    if mode == "evidence":
        # Pure evidence mode stays closest to classic weighted PageRank.
        degree_factor = 1.0
    elif mode == "connectivity":
        # Connectivity mode strongly rewards shared-network structure.
        degree_factor = 1.0 + math.log1p(degree)
    else:
        # Balanced mode gives moderate support to connected pathways while
        # avoiding the "everything is just degree" problem.
        degree_factor = math.sqrt(1.0 + math.log1p(degree))

    return float(pagerank_score) * float(degree_factor)


def run_signal_diffusion(
    g: nx.Graph,
    seed_node: str | None = None,
    alpha: float = 0.85,
    ranking_mode: str = "balanced",
) -> dict:
    """
    Run weighted PageRank / personalized PageRank over the enrichment network.

    Unseeded mode asks: which nodes are globally influential in this network?
    Seeded mode asks: starting from one selected gene/pathway, which nodes are
    most reachable through weighted network diffusion?

    ranking_mode controls how strongly evidence vs connectivity affects the walk.
    """
    if g.number_of_nodes() == 0:
        return {}

    mg = prepare_markov_graph(g, ranking_mode=ranking_mode)
    personalization = None

    if seed_node and seed_node in mg.nodes:
        personalization = {n: 0.0 for n in mg.nodes}
        personalization[seed_node] = 1.0

    try:
        scores = nx.pagerank(
            mg,
            alpha=float(alpha),
            weight="markov_weight",
            personalization=personalization,
            max_iter=500,
            tol=1e-10,
        )
    except nx.PowerIterationFailedConvergence:
        scores = nx.pagerank(
            mg,
            alpha=float(alpha),
            weight="markov_weight",
            personalization=personalization,
            max_iter=2000,
            tol=1e-8,
        )

    return scores


def diffusion_rows(
    g: nx.Graph,
    scores: dict,
    top_n: int = 50,
    ranking_mode: str = "balanced",
) -> list[dict]:
    """Convert diffusion scores into table-friendly rows."""
    rows = []
    for node_id, score in scores.items():
        node_attrs = g.nodes[node_id]
        priority_score = _degree_adjusted_score(g, node_id, float(score), ranking_mode)
        rows.append(
            {
                "label": node_attrs.get("label", node_id),
                "node_type": node_attrs.get("node_type", "unknown"),
                "degree": int(g.degree(node_id)),
                "priority_score": round(float(priority_score), 8),
                "raw_diffusion_score": round(float(score), 8),
                "diffusion_score": round(float(priority_score), 8),
                "ranking_mode": ranking_mode,
                "node_id": node_id,
            }
        )

    rows = sorted(rows, key=lambda r: r["priority_score"], reverse=True)
    for i, row in enumerate(rows, start=1):
        row["rank"] = i

    return rows[: int(top_n)]


def split_diffusion_rows(rows: list[dict], top_n: int = 25) -> tuple[list[dict], list[dict]]:
    """Split diffusion results into pathway/group and gene/item result tables."""
    group_rows = [r for r in rows if r.get("node_type") == "group"][: int(top_n)]
    item_rows = [r for r in rows if r.get("node_type") == "item"][: int(top_n)]
    return group_rows, item_rows


def _minmax01(value: float, values: list[float]) -> float:
    """Scale a value to 0..1 using the current result set."""
    vals = [float(v) for v in values if v is not None and np.isfinite(float(v))]
    if not vals:
        return 0.0
    lo, hi = min(vals), max(vals)
    if hi == lo:
        return 1.0
    return float((float(value) - lo) / (hi - lo))


def _incident_evidence_summary(g: nx.Graph, node_id: str) -> dict:
    """Summarize direct statistical evidence around a node from incident edges."""
    weights = []
    raw_pvals = []

    for _u, _v, ed in g.edges(node_id, data=True):
        weights.append(_normalized_edge_weight(ed))
        if "weight_raw" in ed and ed["weight_raw"] is not None:
            try:
                raw = float(ed["weight_raw"])
                if raw > 0:
                    raw_pvals.append(raw)
            except Exception:
                pass

    if not weights:
        return {
            "max_evidence": 0.0,
            "mean_evidence": 0.0,
            "best_adj_p": None,
        }

    return {
        "max_evidence": float(max(weights)),
        "mean_evidence": float(np.mean(weights)),
        "best_adj_p": float(min(raw_pvals)) if raw_pvals else None,
    }


def candidate_rows(
    g: nx.Graph,
    diffusion_result_rows: list[dict],
    top_n: int = 30,
    node_type: str = "group",
) -> list[dict]:
    """
    Build a follow-up candidate table that combines statistical evidence,
    network diffusion, and support breadth.

    This is meant to answer: what should a scientist inspect first?
    It is not a new p-value. It is a pragmatic prioritization score.
    """
    rows = [r for r in diffusion_result_rows if r.get("node_type") == node_type]
    if not rows:
        return []

    enriched = []
    for r in rows:
        node_id = r.get("node_id")
        if node_id not in g.nodes:
            continue
        ev = _incident_evidence_summary(g, node_id)
        enriched.append({**r, **ev})

    if not enriched:
        return []

    priority_values = [float(r.get("priority_score", 0.0)) for r in enriched]
    evidence_values = [float(r.get("max_evidence", 0.0)) for r in enriched]
    degree_values = [float(r.get("degree", 0.0)) for r in enriched]

    out = []
    for r in enriched:
        diffusion_component = _minmax01(float(r.get("priority_score", 0.0)), priority_values)
        evidence_component = _minmax01(float(r.get("max_evidence", 0.0)), evidence_values)
        support_component = _minmax01(float(r.get("degree", 0.0)), degree_values)

        # Default product-facing score: diffusion is most important, direct evidence
        # keeps it grounded in enrichment strength, and support breadth prevents one-edge
        # terms from dominating when broader biology is present.
        followup_score = (
            0.50 * diffusion_component
            + 0.30 * evidence_component
            + 0.20 * support_component
        )

        best_adj_p = r.get("best_adj_p")
        out.append({
            "label": r.get("label"),
            "node_type": r.get("node_type"),
            "followup_score": round(float(followup_score), 8),
            "priority_score": r.get("priority_score"),
            "raw_diffusion_score": r.get("raw_diffusion_score"),
            "max_evidence": round(float(r.get("max_evidence", 0.0)), 6),
            "mean_evidence": round(float(r.get("mean_evidence", 0.0)), 6),
            "best_adj_p": None if best_adj_p is None else float(best_adj_p),
            "degree": r.get("degree"),
            "ranking_mode": r.get("ranking_mode"),
            "node_id": r.get("node_id"),
        })

    out = sorted(out, key=lambda r: r["followup_score"], reverse=True)[: int(top_n)]
    for i, row in enumerate(out, start=1):
        row["candidate_rank"] = i
    return out



# ----------------------------
# Pathway projection helpers
# ----------------------------

def build_pathway_projection_graph(g: nx.Graph, method: str = "jaccard") -> nx.Graph:
    """Build a pathway-only graph: pathways connect when they share genes/items."""
    method = method or "jaccard"
    pg = nx.Graph()
    groups = [n for n, d in g.nodes(data=True) if d.get("node_type") == "group"]
    items = [n for n, d in g.nodes(data=True) if d.get("node_type") == "item"]

    for gr in groups:
        attrs = g.nodes[gr].copy()
        pg.add_node(
            gr,
            label=attrs.get("label", gr),
            node_type="projected_group",
            source_degree=int(g.degree(gr)),
        )

    group_items = {
        gr: set(nei for nei in g.neighbors(gr) if g.nodes[nei].get("node_type") == "item")
        for gr in groups
    }

    item_groups = {}
    for item in items:
        connected_groups = [nei for nei in g.neighbors(item) if g.nodes[nei].get("node_type") == "group"]
        if len(connected_groups) >= 2:
            item_groups[item] = connected_groups

    pair_payload = {}
    for item, connected_groups in item_groups.items():
        for a, b in combinations(sorted(connected_groups), 2):
            key = (a, b)
            if key not in pair_payload:
                pair_payload[key] = {"shared_items": [], "weighted_support": 0.0}
            wa = _normalized_edge_weight(g.edges[a, item]) if g.has_edge(a, item) else 1.0
            wb = _normalized_edge_weight(g.edges[b, item]) if g.has_edge(b, item) else 1.0
            pair_payload[key]["shared_items"].append(g.nodes[item].get("label", item))
            pair_payload[key]["weighted_support"] += (wa + wb) / 2.0

    for (a, b), payload in pair_payload.items():
        shared_count = len(payload["shared_items"])
        union_count = len(group_items.get(a, set()) | group_items.get(b, set()))
        jaccard = (shared_count / union_count) if union_count else 0.0
        weighted_shared = float(payload["weighted_support"])

        if method == "shared_count":
            edge_weight = float(shared_count)
        elif method == "weighted_shared":
            edge_weight = weighted_shared
        else:
            edge_weight = float(jaccard)

        pg.add_edge(
            a,
            b,
            weight=float(edge_weight),
            shared_count=int(shared_count),
            jaccard=float(jaccard),
            weighted_shared=float(weighted_shared),
            shared_items="; ".join(map(str, payload["shared_items"][:50])),
        )

    isolated = [n for n in pg.nodes if pg.degree(n) == 0]
    pg.remove_nodes_from(isolated)
    return pg


def projection_graph_to_store(pg: nx.Graph) -> dict:
    nodes = [{"id": n, **d} for n, d in pg.nodes(data=True)]
    edges = [{"source": u, "target": v, **ed} for u, v, ed in pg.edges(data=True)]
    return {"nodes": nodes, "edges": edges}


def rebuild_projection_from_store(graph_data: dict) -> nx.Graph:
    pg = nx.Graph()
    if not graph_data:
        return pg
    for nd in graph_data.get("nodes", []):
        node_id = nd.get("id")
        if not node_id:
            continue
        attrs = nd.copy()
        attrs.pop("id", None)
        pg.add_node(node_id, **attrs)
    for ed in graph_data.get("edges", []):
        u, v = ed.get("source"), ed.get("target")
        if not u or not v:
            continue
        attrs = ed.copy()
        attrs.pop("source", None)
        attrs.pop("target", None)
        pg.add_edge(u, v, **attrs)
    return pg


def make_projection_figure(pg: nx.Graph, show_labels: bool = True) -> go.Figure:
    if pg.number_of_nodes() == 0:
        fig = go.Figure()
        fig.update_layout(
            annotations=[dict(text="Build a projection network after building the main graph.", showarrow=False)],
            margin=dict(l=10, r=10, t=10, b=10),
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
        )
        return fig

    pos = nx.spring_layout(pg, seed=11, k=1 / math.sqrt(max(1, pg.number_of_nodes())))
    weights = [float(ed.get("weight", 1.0)) for _, _, ed in pg.edges(data=True)]

    if weights:
        w_arr = np.asarray(weights, dtype=float)
        lo, hi = float(np.percentile(w_arr, 5)), float(np.percentile(w_arr, 95))
        if hi <= lo:
            lo, hi = float(w_arr.min()), float(w_arr.max())
    else:
        lo, hi = 0.0, 1.0

    def width_from_weight(w: float) -> float:
        if hi > lo:
            t = max(0.0, min(1.0, (float(w) - lo) / (hi - lo)))
        else:
            t = 0.5
        return 1.0 + t * 5.0

    edge_traces = []
    for u, v, ed in pg.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        hover = (
            f"{pg.nodes[u].get('label', u)} ↔ {pg.nodes[v].get('label', v)}"
            f"<br>weight={float(ed.get('weight', 0.0)):.4g}"
            f"<br>shared items={ed.get('shared_count', 0)}"
            f"<br>jaccard={float(ed.get('jaccard', 0.0)):.4g}"
        )
        edge_traces.append(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode="lines",
                hovertext=hover,
                hoverinfo="text",
                line=dict(width=width_from_weight(ed.get("weight", 1.0)), color="rgba(100,116,139,0.55)"),
                showlegend=False,
            )
        )

    pr = nx.pagerank(pg, weight="weight") if pg.number_of_edges() else {
        n: 1 / max(1, pg.number_of_nodes()) for n in pg.nodes
    }
    max_pr = max(pr.values()) if pr else 1.0

    node_x, node_y, node_text, node_hover, node_size = [], [], [], [], []
    for n, d in pg.nodes(data=True):
        x, y = pos[n]
        label = str(d.get("label", n))
        deg = pg.degree(n)
        node_x.append(x)
        node_y.append(y)
        node_text.append(label if show_labels else "")
        node_hover.append(
            f"{label}<br>projection degree={deg}"
            f"<br>source membership degree={d.get('source_degree', '—')}"
            f"<br>PageRank={pr.get(n, 0.0):.6f}"
        )
        node_size.append(12 + 26 * (pr.get(n, 0.0) / max_pr if max_pr else 0.0))

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text" if show_labels else "markers",
        text=node_text,
        textposition="top center",
        hovertext=node_hover,
        hoverinfo="text",
        marker=dict(
            size=node_size,
            color="rgba(0,119,182,0.92)",
            line=dict(width=1.5, color="rgba(0,0,0,0.35)"),
        ),
        name="pathways",
    )

    fig = go.Figure(data=edge_traces + [node_trace])
    fig.update_layout(
        margin=dict(l=10, r=10, t=10, b=10),
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        dragmode="pan",
        hovermode="closest",
    )
    return fig


def projection_stats_rows(pg: nx.Graph, top_n: int = 30) -> list[dict]:
    if pg.number_of_nodes() == 0:
        return []
    pr = nx.pagerank(pg, weight="weight") if pg.number_of_edges() else {n: 1 / pg.number_of_nodes() for n in pg.nodes}
    rows = []
    for n, score in sorted(pr.items(), key=lambda x: x[1], reverse=True):
        rows.append({
            "rank": len(rows) + 1,
            "pathway": pg.nodes[n].get("label", n),
            "projection_score": round(float(score), 8),
            "projection_degree": int(pg.degree(n)),
            "source_degree": int(pg.nodes[n].get("source_degree", 0)),
            "node_id": n,
        })
    return rows[:int(top_n)]


def prepare_projection_markov_graph(pg: nx.Graph, ranking_mode: str = "balanced") -> nx.Graph:
    mg = pg.copy()
    mode = (ranking_mode or "balanced").lower().strip()
    for _u, _v, ed in mg.edges(data=True):
        try:
            overlap_w = float(ed.get("weight", 1.0))
        except Exception:
            overlap_w = 1.0
        overlap_w = max(overlap_w, 1e-12)
        if mode == "connectivity":
            markov_w = 1.0
        elif mode == "evidence":
            markov_w = overlap_w
        else:
            markov_w = math.sqrt(overlap_w)
        ed["markov_weight"] = max(float(markov_w), 1e-12)
    return mg


def run_projection_diffusion(pg: nx.Graph, seed_node: str | None = None, alpha: float = 0.85, ranking_mode: str = "balanced") -> dict:
    if pg.number_of_nodes() == 0:
        return {}
    mg = prepare_projection_markov_graph(pg, ranking_mode=ranking_mode)
    personalization = None
    if seed_node and seed_node in mg.nodes:
        personalization = {n: 0.0 for n in mg.nodes}
        personalization[seed_node] = 1.0
    try:
        return nx.pagerank(mg, alpha=float(alpha), weight="markov_weight", personalization=personalization, max_iter=500, tol=1e-10)
    except nx.PowerIterationFailedConvergence:
        return nx.pagerank(mg, alpha=float(alpha), weight="markov_weight", personalization=personalization, max_iter=2000, tol=1e-8)


def _projection_degree_adjusted_score(pg: nx.Graph, node_id: str, pagerank_score: float, ranking_mode: str) -> float:
    mode = (ranking_mode or "balanced").lower().strip()
    degree = max(int(pg.degree(node_id)), 1)
    if mode == "evidence":
        degree_factor = 1.0
    elif mode == "connectivity":
        degree_factor = 1.0 + math.log1p(degree)
    else:
        degree_factor = math.sqrt(1.0 + math.log1p(degree))
    return float(pagerank_score) * float(degree_factor)


def projection_diffusion_rows(pg: nx.Graph, scores: dict, top_n: int = 50, ranking_mode: str = "balanced") -> list[dict]:
    rows = []
    for node_id, score in scores.items():
        if node_id not in pg.nodes:
            continue
        attrs = pg.nodes[node_id]
        priority_score = _projection_degree_adjusted_score(pg, node_id, float(score), ranking_mode)
        rows.append({
            "label": attrs.get("label", node_id),
            "node_type": "projected_group",
            "degree": int(pg.degree(node_id)),
            "source_degree": int(attrs.get("source_degree", 0)),
            "priority_score": round(float(priority_score), 8),
            "raw_diffusion_score": round(float(score), 8),
            "diffusion_score": round(float(priority_score), 8),
            "ranking_mode": ranking_mode,
            "node_id": node_id,
        })
    rows = sorted(rows, key=lambda r: r["priority_score"], reverse=True)
    for i, row in enumerate(rows, start=1):
        row["rank"] = i
    return rows[: int(top_n)]


def _projection_overlap_summary(pg: nx.Graph, node_id: str) -> dict:
    weights, shared_counts = [], []
    for _u, _v, ed in pg.edges(node_id, data=True):
        try:
            weights.append(float(ed.get("weight", 0.0)))
        except Exception:
            pass
        try:
            shared_counts.append(float(ed.get("shared_count", 0.0)))
        except Exception:
            pass
    if not weights:
        return {"mean_overlap_weight": 0.0, "max_overlap_weight": 0.0, "mean_shared_count": 0.0, "max_shared_count": 0.0}
    return {
        "mean_overlap_weight": float(np.mean(weights)),
        "max_overlap_weight": float(max(weights)),
        "mean_shared_count": float(np.mean(shared_counts)) if shared_counts else 0.0,
        "max_shared_count": float(max(shared_counts)) if shared_counts else 0.0,
    }


def projection_candidate_rows(pg: nx.Graph, diffusion_result_rows: list[dict], top_n: int = 30) -> list[dict]:
    if not diffusion_result_rows:
        return []
    enriched = []
    for r in diffusion_result_rows:
        node_id = r.get("node_id")
        if node_id not in pg.nodes:
            continue
        enriched.append({**r, **_projection_overlap_summary(pg, node_id)})
    if not enriched:
        return []
    priority_values = [float(r.get("priority_score", 0.0)) for r in enriched]
    degree_values = [float(r.get("degree", 0.0)) for r in enriched]
    overlap_values = [float(r.get("mean_overlap_weight", 0.0)) for r in enriched]
    source_degree_values = [float(r.get("source_degree", 0.0)) for r in enriched]
    out = []
    for r in enriched:
        followup_score = (0.50 * _minmax01(float(r.get("priority_score", 0.0)), priority_values)
                          + 0.25 * _minmax01(float(r.get("degree", 0.0)), degree_values)
                          + 0.15 * _minmax01(float(r.get("mean_overlap_weight", 0.0)), overlap_values)
                          + 0.10 * _minmax01(float(r.get("source_degree", 0.0)), source_degree_values))
        out.append({
            "label": r.get("label"),
            "node_type": r.get("node_type"),
            "followup_score": round(float(followup_score), 8),
            "priority_score": r.get("priority_score"),
            "raw_diffusion_score": r.get("raw_diffusion_score"),
            "degree": r.get("degree"),
            "source_degree": r.get("source_degree"),
            "mean_overlap_weight": round(float(r.get("mean_overlap_weight", 0.0)), 6),
            "max_overlap_weight": round(float(r.get("max_overlap_weight", 0.0)), 6),
            "mean_shared_count": round(float(r.get("mean_shared_count", 0.0)), 6),
            "max_shared_count": round(float(r.get("max_shared_count", 0.0)), 6),
            "ranking_mode": r.get("ranking_mode"),
            "node_id": r.get("node_id"),
        })
    out = sorted(out, key=lambda r: r["followup_score"], reverse=True)[: int(top_n)]
    for i, row in enumerate(out, start=1):
        row["candidate_rank"] = i
    return out


# ----------------------------
# Consensus candidate helpers
# ----------------------------

def consensus_candidate_rows(bipartite_store: dict | None, projection_store: dict | None, top_n: int = 30) -> list[dict]:
    """Combine bipartite and projection candidate scores into one final pathway list."""
    if not bipartite_store or not projection_store:
        return []

    bip_rows = bipartite_store.get("candidate_rows") or []
    proj_rows = projection_store.get("candidate_rows") or []

    bip_by_key = {}
    for r in bip_rows:
        key = r.get("node_id") or r.get("label")
        if key:
            bip_by_key[key] = r

    proj_by_key = {}
    for r in proj_rows:
        key = r.get("node_id") or r.get("label")
        if key:
            proj_by_key[key] = r

    keys = sorted(set(bip_by_key) | set(proj_by_key))
    if not keys:
        return []

    bip_scores = [float(r.get("followup_score", 0.0)) for r in bip_by_key.values()]
    proj_scores = [float(r.get("followup_score", 0.0)) for r in proj_by_key.values()]

    out = []
    for key in keys:
        b = bip_by_key.get(key, {})
        pr = proj_by_key.get(key, {})
        label = b.get("label") or pr.get("label") or str(key)

        b_score_raw = float(b.get("followup_score", 0.0)) if b else 0.0
        p_score_raw = float(pr.get("followup_score", 0.0)) if pr else 0.0

        b_score_norm = _minmax01(b_score_raw, bip_scores) if bip_scores else 0.0
        p_score_norm = _minmax01(p_score_raw, proj_scores) if proj_scores else 0.0
        consensus_score = 0.50 * b_score_norm + 0.50 * p_score_norm

        out.append({
            "label": label,
            "consensus_score": round(float(consensus_score), 8),
            "bipartite_followup_score": round(float(b_score_raw), 8),
            "projection_followup_score": round(float(p_score_raw), 8),
            "bipartite_candidate_rank": b.get("candidate_rank"),
            "projection_candidate_rank": pr.get("candidate_rank"),
            "bipartite_priority_score": b.get("priority_score"),
            "projection_priority_score": pr.get("priority_score"),
            "best_adj_p": b.get("best_adj_p"),
            "max_evidence": b.get("max_evidence"),
            "bipartite_degree": b.get("degree"),
            "projection_degree": pr.get("degree"),
            "source_degree": pr.get("source_degree"),
            "mean_overlap_weight": pr.get("mean_overlap_weight"),
            "max_shared_count": pr.get("max_shared_count"),
            "node_id": key,
        })

    out = sorted(out, key=lambda r: r["consensus_score"], reverse=True)[: int(top_n)]
    for i, row in enumerate(out, start=1):
        row["consensus_rank"] = i
    return out


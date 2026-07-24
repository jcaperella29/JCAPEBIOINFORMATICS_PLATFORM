
from __future__ import annotations

import io
import json
import math
import zipfile
from datetime import datetime, timezone
from typing import Any

import networkx as nx
import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import StreamingResponse

from . import core
from .schemas import (
    BuildGraphOptions,
    ConsensusRequest,
    DiffusionRequest,
    FilterGraphRequest,
    ProjectionBuildRequest,
    ProjectionDiffusionRequest,
)

app = FastAPI(
    title="Enrichment Network API",
    version="0.1.0",
    description="API module for building enrichment gene/pathway networks, running diffusion ranking, building pathway projections, and exporting dashboard-ready outputs.",
)


def _clean(value: Any) -> Any:
    """Make nested pandas/numpy/NaN values JSON-safe."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_clean(v) for v in value]
    if hasattr(value, "item"):
        try:
            return _clean(value.item())
        except Exception:
            pass
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def graph_to_store(g: nx.Graph) -> dict[str, Any]:
    return {
        "nodes": [{"id": n, **d} for n, d in g.nodes(data=True)],
        "edges": [{"source": u, "target": v, **ed} for u, v, ed in g.edges(data=True)],
    }


def store_to_csvs(store: dict[str, Any]) -> dict[str, str]:
    nodes = pd.DataFrame(store.get("nodes", []))
    edges = pd.DataFrame(store.get("edges", []))
    return {
        "nodes_csv": nodes.to_csv(index=False),
        "edges_csv": edges.to_csv(index=False),
    }


def stats_tables(g: nx.Graph) -> dict[str, Any]:
    st = core.graph_stats(g)
    return {
        "stats": st,
        "top_groups": [{"group": name, "degree": degree} for name, degree in st.get("top_groups", [])],
        "top_items": [{"item": name, "degree": degree} for name, degree in st.get("top_items", [])],
    }


def figure_json(fig) -> dict[str, Any]:
    return json.loads(fig.to_json())


async def upload_to_dataframe(file: UploadFile) -> pd.DataFrame:
    raw = await file.read()
    try:
        return pd.read_csv(io.BytesIO(raw))
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Could not parse CSV: {exc}") from exc


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "enrichment-network-api", "version": app.version}


@app.get("/presets")
def presets() -> dict[str, Any]:
    return {"column_mapping_presets": core.COLUMN_MAPPING_PRESETS}


@app.post("/network/build")
async def build_network(
    file: UploadFile = File(...),
    options_json: str = Form(default='{"preset":"custom_long"}'),
) -> dict[str, Any]:
    """
    Upload a CSV and build the main bipartite gene/item ↔ pathway/group graph.

    options_json example:
    {"preset":"enrichr", "apply_preset":true}
    or
    {"apply_preset":false, "item_col":"gene", "group_col":"term", "weight_col":"adjusted_pvalue"}
    """
    try:
        options = BuildGraphOptions.model_validate_json(options_json)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Invalid options_json: {exc}") from exc

    df = await upload_to_dataframe(file)
    preset_msg = None

    try:
        if options.apply_preset:
            df, guessed_item, guessed_group, guessed_weight, preset_msg = core.apply_column_mapping_preset(df, options.preset)
            item_col = options.item_col or guessed_item
            group_col = options.group_col or guessed_group
            weight_col = options.weight_col if options.weight_col is not None else guessed_weight
        else:
            item_col = options.item_col
            group_col = options.group_col
            weight_col = options.weight_col

        if not item_col or not group_col:
            raise HTTPException(status_code=400, detail="item_col and group_col are required or must be guessable from the preset.")
        if item_col not in df.columns or group_col not in df.columns:
            raise HTTPException(status_code=400, detail=f"Selected columns not present. Available columns: {list(df.columns)}")
        if weight_col == "":
            weight_col = None
        if weight_col and weight_col not in df.columns:
            raise HTTPException(status_code=400, detail=f"weight_col '{weight_col}' not present. Available columns: {list(df.columns)}")

        g = core.build_bipartite_graph(df, item_col=item_col, group_col=group_col, weight_col=weight_col)
        store = graph_to_store(g)
        out = {
            "message": f"Graph built with {g.number_of_nodes():,} nodes and {g.number_of_edges():,} edges.",
            "preset_message": preset_msg,
            "mapped_columns": {"item_col": item_col, "group_col": group_col, "weight_col": weight_col},
            "input_shape": {"rows": int(df.shape[0]), "columns": int(df.shape[1])},
            "graph": store,
            "exports": store_to_csvs(store),
            **stats_tables(g),
            "warning": core.graph_size_warning_component(g.number_of_nodes(), g.number_of_edges(), context="raw graph"),
        }
        return _clean(out)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.post("/network/filter")
def filter_network(req: FilterGraphRequest) -> dict[str, Any]:
    g = core.rebuild_graph_from_store(req.graph)
    opt = req.options
    sg = core.subgraph_filter(
        g,
        search=opt.search,
        min_degree=opt.min_degree,
        min_weight=opt.min_weight,
        max_groups=opt.max_groups,
        largest_component_only=opt.largest_component_only,
    )
    store = graph_to_store(sg)
    out: dict[str, Any] = {
        "message": f"Filtered graph has {sg.number_of_nodes():,} nodes and {sg.number_of_edges():,} edges.",
        "graph": store,
        "exports": store_to_csvs(store),
        **stats_tables(sg),
        "warning": core.graph_size_warning_component(sg.number_of_nodes(), sg.number_of_edges(), context="filtered graph"),
    }
    if opt.return_figure:
        if opt.layout_mode == "force":
            pos = nx.spring_layout(sg, seed=7, k=1 / math.sqrt(max(1, sg.number_of_nodes())))
        else:
            pos = core.layout_bipartite_two_column(sg)
        fig = core.make_plotly_network(
            sg,
            pos,
            show_labels=opt.show_labels,
            thickness_by_weight=opt.thickness_by_weight,
            highlight_nodes=req.highlight_nodes,
        )
        out["figure"] = figure_json(fig)
    return _clean(out)


@app.post("/diffusion/bipartite")
def bipartite_diffusion(req: DiffusionRequest) -> dict[str, Any]:
    g = core.rebuild_graph_from_store(req.graph)
    scores = core.run_signal_diffusion(g, seed_node=req.seed_node, alpha=req.alpha, ranking_mode=req.ranking_mode)
    rows = core.diffusion_rows(g, scores, top_n=req.top_n, ranking_mode=req.ranking_mode)
    group_rows, item_rows = core.split_diffusion_rows(rows, top_n=req.top_n)
    candidate_rows = core.candidate_rows(g, rows, top_n=req.candidate_top_n, node_type=req.candidate_node_type)
    return _clean({
        "seed_node": req.seed_node,
        "alpha": req.alpha,
        "ranking_mode": req.ranking_mode,
        "rows": rows,
        "group_rows": group_rows,
        "item_rows": item_rows,
        "candidate_rows": candidate_rows,
        "csv": pd.DataFrame(rows).to_csv(index=False),
        "candidate_csv": pd.DataFrame(candidate_rows).to_csv(index=False),
    })


@app.post("/projection/build")
def build_projection(req: ProjectionBuildRequest) -> dict[str, Any]:
    g = core.rebuild_graph_from_store(req.graph)
    pg = core.build_pathway_projection_graph(g, method=req.method)
    store = core.projection_graph_to_store(pg)
    nodes_df = pd.DataFrame(store.get("nodes", []))
    edges_df = pd.DataFrame(store.get("edges", []))
    out: dict[str, Any] = {
        "message": f"Projection built with {pg.number_of_nodes():,} pathway nodes and {pg.number_of_edges():,} overlap edges.",
        "projection_graph": store,
        "projection_rows": core.projection_stats_rows(pg),
        "exports": {"nodes_csv": nodes_df.to_csv(index=False), "edges_csv": edges_df.to_csv(index=False)},
    }
    if req.return_figure:
        out["figure"] = figure_json(core.make_projection_figure(pg, show_labels=req.show_labels))
    return _clean(out)


@app.post("/diffusion/projection")
def projection_diffusion(req: ProjectionDiffusionRequest) -> dict[str, Any]:
    pg = core.rebuild_projection_from_store(req.projection_graph)
    scores = core.run_projection_diffusion(pg, seed_node=req.seed_node, alpha=req.alpha, ranking_mode=req.ranking_mode)
    rows = core.projection_diffusion_rows(pg, scores, top_n=req.top_n, ranking_mode=req.ranking_mode)
    candidate_rows = core.projection_candidate_rows(pg, rows, top_n=req.candidate_top_n)
    return _clean({
        "seed_node": req.seed_node,
        "alpha": req.alpha,
        "ranking_mode": req.ranking_mode,
        "rows": rows,
        "candidate_rows": candidate_rows,
        "csv": pd.DataFrame(rows).to_csv(index=False),
        "candidate_csv": pd.DataFrame(candidate_rows).to_csv(index=False),
    })


@app.post("/consensus")
def consensus(req: ConsensusRequest) -> dict[str, Any]:
    rows = core.consensus_candidate_rows(req.bipartite_results, req.projection_results, top_n=req.top_n)
    return _clean({"candidate_rows": rows, "csv": pd.DataFrame(rows).to_csv(index=False)})


@app.post("/export/report-bundle")
def export_report_bundle(payload: dict[str, Any]) -> StreamingResponse:
    """
    Build a simple ZIP from already computed API outputs.

    Expected payload keys are flexible: graph, projection_graph, bipartite_results,
    projection_results, consensus_results, mapped_columns, settings.
    """
    graph = payload.get("graph") or {}
    projection_graph = payload.get("projection_graph") or {}
    bip = payload.get("bipartite_results") or {}
    proj = payload.get("projection_results") or {}
    con = payload.get("consensus_results") or {}

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("manifest.json", json.dumps({
            "created_at": datetime.now(timezone.utc).isoformat(),
            "product": "Enrichment Network API",
            "mapped_columns": payload.get("mapped_columns", {}),
            "settings": payload.get("settings", {}),
            "caveat": "Diffusion, follow-up, and consensus scores are prioritization scores, not statistical p-values.",
        }, indent=2))
        if graph:
            csvs = store_to_csvs(graph)
            z.writestr("main_nodes.csv", csvs["nodes_csv"])
            z.writestr("main_edges.csv", csvs["edges_csv"])
        if projection_graph:
            z.writestr("projection_nodes.csv", pd.DataFrame(projection_graph.get("nodes", [])).to_csv(index=False))
            z.writestr("projection_edges.csv", pd.DataFrame(projection_graph.get("edges", [])).to_csv(index=False))
        z.writestr("bipartite_diffusion_results.csv", bip.get("csv", ""))
        z.writestr("bipartite_top_candidates.csv", bip.get("candidate_csv", ""))
        z.writestr("projection_diffusion_results.csv", proj.get("csv", ""))
        z.writestr("projection_top_candidates.csv", proj.get("candidate_csv", ""))
        z.writestr("consensus_candidates.csv", con.get("csv", ""))
        z.writestr("interpretation_notes.md", """# Interpretation Notes\n\nDiffusion/follow-up/consensus scores are prioritization scores, not p-values. Use original adjusted p-values as statistical evidence.\n""")
    buf.seek(0)
    return StreamingResponse(
        buf,
        media_type="application/zip",
        headers={"Content-Disposition": "attachment; filename=enrichment_network_report_bundle.zip"},
    )

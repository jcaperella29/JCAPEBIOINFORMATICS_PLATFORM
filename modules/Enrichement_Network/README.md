# Enrichment Network API Module

FastAPI wrapper for the Dash enrichment-network logic. It turns the current Dash app into a reusable API module for a larger multi-omics dashboard.

## What this module exposes

- Build a bipartite gene/item ↔ pathway/group graph from uploaded CSVs
- Apply column presets for common enrichment outputs:
  - custom long-format CSV
  - Enrichr-style `term + genes`
  - g:Profiler/gprofiler2-style
  - clusterProfiler-style
  - GSEA/MSigDB-style
- Filter graph by search, minimum degree, edge weight, max groups, and component mode
- Return Plotly figure JSON for direct Dash rendering
- Run bipartite signal diffusion / PageRank-style prioritization
- Build pathway-only projection graphs from shared genes/items
- Run projection diffusion
- Build consensus pathway candidates
- Export a ZIP report bundle from computed outputs

## Install

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Run the API

```bash
uvicorn enrichment_network_api.main:app --reload --port 8000
```

Open:

```text
http://localhost:8000/docs
```

## Build a graph from CSV

```bash
curl -X POST "http://localhost:8000/network/build" \
  -F "file=@TEST.csv" \
  -F 'options_json={"preset":"enrichr","apply_preset":true}'
```

For an already-long table:

```bash
curl -X POST "http://localhost:8000/network/build" \
  -F "file=@edges.csv" \
  -F 'options_json={"apply_preset":false,"item_col":"gene","group_col":"term","weight_col":"adjusted_pvalue"}'
```

The response includes:

- `graph`: serializable nodes/edges store
- `stats`: graph statistics
- `exports.nodes_csv`, `exports.edges_csv`
- guessed/applied `mapped_columns`

## Dash integration pattern

In the dashboard, call `/network/build`, store the returned `graph` in a `dcc.Store`, and pass that graph to downstream endpoints.

For plotting, call `/network/filter` with `return_figure=true` and feed the returned `figure` directly to `dcc.Graph(figure=...)`.

Example callback sketch:

```python
import requests

API = "http://localhost:8000"

resp = requests.post(
    f"{API}/network/build",
    files={"file": (filename, csv_bytes, "text/csv")},
    data={"options_json": json.dumps({"preset": "enrichr", "apply_preset": True})},
)
data = resp.json()
graph_store = data["graph"]

plot_resp = requests.post(
    f"{API}/network/filter",
    json={
        "graph": graph_store,
        "options": {
            "search": "",
            "min_degree": 0,
            "min_weight": 0,
            "max_groups": 50,
            "layout_mode": "bipartite",
            "return_figure": True,
        },
    },
)
figure = plot_resp.json()["figure"]
```

## Suggested dashboard tab mapping

- Upload/preprocess tab → `/network/build`
- Main plot tab → `/network/filter`
- Stats tab → `stats`, `top_groups`, `top_items` from build/filter responses
- Bipartite Diffusion tab → `/diffusion/bipartite`
- Projection tab → `/projection/build`
- Projection Diffusion tab → `/diffusion/projection`
- Consensus tab → `/consensus`
- Report button → `/export/report-bundle`

## Notes

Diffusion, follow-up, and consensus scores are prioritization scores, not p-values. Use the original enrichment adjusted p-values as statistical evidence.

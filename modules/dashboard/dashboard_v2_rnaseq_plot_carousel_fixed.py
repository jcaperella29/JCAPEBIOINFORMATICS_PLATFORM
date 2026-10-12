from __future__ import annotations

from dash import Dash, Input, Output, State, dcc, html, dash_table, no_update, ctx
from pathlib import Path
import base64
import csv
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone


# -----------------------------------------------------------------------------
# App
# -----------------------------------------------------------------------------
app = Dash(__name__, suppress_callback_exceptions=True, title="JCAP Bioinformatics Platform")
server = app.server

# -----------------------------------------------------------------------------
# Integrated-pipeline execution configuration
# -----------------------------------------------------------------------------
RNASEQ_PIPELINE_DIR = Path(os.environ.get(
    "JCAP_RNASEQ_PIPELINE_DIR",
    "/mnt/c/Users/jcape/Downloads/omics_pipeline/modules/jcap_rnaseq_integrated_pipeline",
))
RNASEQ_PARAMS_FILE = Path(os.environ.get(
    "JCAP_RNASEQ_PARAMS_FILE",
    str(RNASEQ_PIPELINE_DIR / "params.json"),
))
RNASEQ_TEST_COUNTS = Path(os.environ.get(
    "JCAP_RNASEQ_TEST_COUNTS",
    "/mnt/c/Users/jcape/Downloads/omics_pipeline/modules/RNA_SEQ/sample_data(airways)/counts_data.csv",
))
RNASEQ_TEST_PHENOTYPE = Path(os.environ.get(
    "JCAP_RNASEQ_TEST_PHENOTYPE",
    "/mnt/c/Users/jcape/Downloads/omics_pipeline/modules/RNA_SEQ/sample_data(airways)/phenotype_data.csv",
))
DASHBOARD_RUN_ROOT = Path(os.environ.get(
    "JCAP_DASHBOARD_RUN_ROOT",
    str(Path.home() / "jcap_dashboard_runs"),
))
DASHBOARD_RUN_ROOT.mkdir(parents=True, exist_ok=True)
ACTIVE_PROCESSES: dict[str, subprocess.Popen] = {}

def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()

def decode_dash_upload(contents: str | None, destination: Path) -> Path:
    if not contents:
        raise ValueError("No uploaded file data received.")
    _header, encoded = contents.split(",", 1)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(base64.b64decode(encoded))
    return destination

def read_run_meta(run_dir: Path) -> dict:
    path = run_dir / "dashboard_run.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}

def write_run_meta(run_dir: Path, payload: dict) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "dashboard_run.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )

def stage_completed_rnaseq_results(run_dir: Path) -> Path | None:
    source = RNASEQ_PIPELINE_DIR / "results" / "jcap_integrated_results"
    if not source.exists():
        return None

    destination = run_dir / "results" / "jcap_integrated_results"
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Make staging safe if a refresh happens twice or the result tree contains
    # symlinked directories. Re-running this copy should not crash the dashboard.
    shutil.copytree(
        source,
        destination,
        dirs_exist_ok=True,
        symlinks=True,
    )

    zip_source = RNASEQ_PIPELINE_DIR / "results" / "jcap_integrated_results.zip"
    if zip_source.exists():
        shutil.copy2(zip_source, run_dir / "results" / zip_source.name)

    return destination


# -----------------------------------------------------------------------------
# Styling
# -----------------------------------------------------------------------------
COLORS = {
    "bg": "#f4f7fb",
    "card": "#ffffff",
    "navy": "#173b67",
    "blue": "#235a97",
    "blue2": "#eaf2fb",
    "border": "#d7e0eb",
    "text": "#172033",
    "muted": "#66758a",
    "good": "#eaf7ef",
    "warn": "#fff7df",
}


def card(children, **kwargs):
    style = {
        "background": COLORS["card"],
        "border": f"1px solid {COLORS['border']}",
        "borderRadius": "16px",
        "boxShadow": "0 6px 22px rgba(27, 50, 84, 0.08)",
        "padding": "18px",
        "marginBottom": "18px",
    }
    style.update(kwargs.pop("style", {}))
    return html.Div(children, style=style, **kwargs)


def upload_box(id_: str, label: str):
    return dcc.Upload(
        id=id_,
        multiple=False,
        children=html.Div(
            [
                html.Div("＋", style={"fontSize": "28px", "fontWeight": "700"}),
                html.Div(label, style={"fontWeight": "700"}),
                html.Div("Drag/drop or click to select", style={"fontSize": "12px", "color": COLORS["muted"]}),
            ],
            style={"textAlign": "center"},
        ),
        style={
            "border": f"2px dashed {COLORS['border']}",
            "borderRadius": "12px",
            "padding": "18px",
            "background": "#fbfdff",
            "cursor": "pointer",
        },
    )


def section_title(title: str, subtitle: str | None = None):
    children = [html.H2(title, style={"margin": "0 0 4px 0", "color": COLORS["navy"]})]
    if subtitle:
        children.append(html.Div(subtitle, style={"color": COLORS["muted"], "marginBottom": "12px"}))
    return html.Div(children)


def placeholder(text: str):
    return html.Div(
        text,
        style={
            "padding": "26px",
            "border": f"1px dashed {COLORS['border']}",
            "borderRadius": "12px",
            "background": "#fbfdff",
            "color": COLORS["muted"],
        },
    )


# -----------------------------------------------------------------------------
# Get Started
# -----------------------------------------------------------------------------
def get_started_layout():
    return card(
        [
            section_title(
                "Get Started",
                "Upload data, choose an integrated pipeline, run it, and track the job from one place.",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Pipeline", style={"fontWeight": "700"}),
                            dcc.Dropdown(
                                id="pipeline-picker",
                                value="rnaseq",
                                clearable=False,
                                options=[
                                    {"label": "Bulk RNA-seq", "value": "rnaseq"},
                                    {"label": "scRNA-seq", "value": "scrna"},
                                    {"label": "Perturb-seq / Mixscape", "value": "perturbseq"},
                                ],
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            html.Label("Run mode", style={"fontWeight": "700"}),
                            dcc.RadioItems(
                                id="run-mode",
                                value="uploaded",
                                options=[
                                    {"label": " Uploaded data", "value": "uploaded"},
                                    {"label": " Built-in test data", "value": "test"},
                                ],
                                inline=True,
                                style={"paddingTop": "8px"},
                            ),
                        ]
                    ),
                ],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "minmax(260px, 1fr) minmax(300px, 1fr)",
                    "gap": "16px",
                    "marginBottom": "16px",
                },
            ),
            html.Div(id="pipeline-upload-panel"),
            html.Div(
                [
                    html.Button(
                        "Run Integrated Pipeline",
                        id="run-integrated-btn",
                        n_clicks=0,
                        style={
                            "background": COLORS["blue"],
                            "color": "white",
                            "border": "none",
                            "borderRadius": "9px",
                            "padding": "11px 18px",
                            "fontWeight": "800",
                            "cursor": "pointer",
                            "marginRight": "10px",
                        },
                    ),
                    html.Button(
                        "Refresh Job Status",
                        id="refresh-job-btn",
                        n_clicks=0,
                        style={
                            "background": "white",
                            "color": COLORS["navy"],
                            "border": f"1px solid {COLORS['border']}",
                            "borderRadius": "9px",
                            "padding": "10px 18px",
                            "fontWeight": "700",
                            "cursor": "pointer",
                        },
                    ),
                ],
                style={"marginTop": "16px"},
            ),
            html.Div(
                id="job-status",
                children="No job started yet.",
                style={
                    "marginTop": "14px",
                    "padding": "12px",
                    "background": COLORS["blue2"],
                    "borderRadius": "10px",
                    "color": COLORS["navy"],
                },
            ),
            dcc.Store(id="active-run-store", storage_type="session"),
            dcc.Store(id="rna-omics-plot-index", data=0, storage_type="session"),
        ]
    )


@app.callback(
    Output("pipeline-upload-panel", "children"),
    Input("pipeline-picker", "value"),
    Input("run-mode", "value"),
)
def render_upload_panel(pipeline, run_mode):
    if run_mode == "test":
        labels = {
            "rnaseq": "Bulk RNA-seq test dataset",
            "scrna": "scRNA-seq test dataset",
            "perturbseq": "Perturb-seq / Mixscape test dataset",
        }
        return html.Div(
            [
                html.Div(
                    f"{labels[pipeline]} will be used. No upload required.",
                    style={
                        "padding": "16px",
                        "background": COLORS["good"],
                        "borderRadius": "10px",
                        "fontWeight": "700",
                    },
                ),
                # These stay mounted (hidden) because the launch callback declares
                # them as State inputs even when built-in test data are selected.
                html.Div(
                    [
                        dcc.Upload(id="run-counts-upload", multiple=False),
                        dcc.Upload(id="run-metadata-upload", multiple=False),
                    ],
                    style={"display": "none"},
                ),
            ]
        )

    if pipeline == "rnaseq":
        return html.Div(
            [
                html.Div([html.Label("Counts matrix", style={"fontWeight": "700"}), upload_box("run-counts-upload", "Upload counts CSV")]),
                html.Div([html.Label("Phenotype metadata", style={"fontWeight": "700"}), upload_box("run-metadata-upload", "Upload phenotype CSV")]),
            ],
            style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "16px"},
        )

    if pipeline == "scrna":
        return html.Div(
            [
                html.Div([html.Label("Counts matrix", style={"fontWeight": "700"}), upload_box("run-counts-upload", "Upload scRNA counts CSV")]),
                html.Div([html.Label("Cell metadata", style={"fontWeight": "700"}), upload_box("run-metadata-upload", "Upload cell metadata CSV")]),
            ],
            style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "16px"},
        )

    return html.Div(
        [
            html.Div([html.Label("Counts matrix", style={"fontWeight": "700"}), upload_box("run-counts-upload", "Upload Perturb-seq counts CSV")]),
            html.Div([html.Label("Cell / perturbation metadata", style={"fontWeight": "700"}), upload_box("run-metadata-upload", "Upload perturbation metadata CSV")]),
        ],
        style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "16px"},
    )


# -----------------------------------------------------------------------------
# Reusable output-family helpers
# -----------------------------------------------------------------------------
def export_bar(prefix: str):
    return html.Div(
        [
            html.Button("Download current table", id=f"{prefix}-download-table", n_clicks=0, className="export-btn"),
            html.Button("Download current plot", id=f"{prefix}-download-plot", n_clicks=0, className="export-btn"),
            html.Button("Download section files", id=f"{prefix}-download-section", n_clicks=0, className="export-btn"),
            html.Button("Download full pipeline bundle", id=f"{prefix}-download-bundle", n_clicks=0, className="export-btn"),
            dcc.Download(id=f"{prefix}-download"),
        ],
        style={"display": "flex", "gap": "8px", "flexWrap": "wrap", "marginBottom": "12px"},
    )


def nested_output_tabs(prefix: str, labels: list[tuple[str, str]], default: str):
    return dcc.Tabs(
        id=f"{prefix}-inner-tabs",
        value=default,
        children=[dcc.Tab(label=label, value=value) for label, value in labels],
    )


def output_panel(prefix: str, title: str, inner_tabs):
    return card(
        [
            html.H3(title, style={"color": COLORS["navy"], "marginTop": 0}),
            export_bar(prefix),
            inner_tabs,
            html.Div(
                id=f"{prefix}-inner-panel",
                children=placeholder("Run or load a pipeline result to populate this section."),
                style={"marginTop": "14px"},
            ),
        ]
    )



# -----------------------------------------------------------------------------
# Result rendering helpers
# -----------------------------------------------------------------------------
def _run_result_root(store: dict | None) -> Path | None:
    if not store:
        return None
    p = store.get("result_dir")
    if not p:
        return None
    root = Path(p)
    return root if root.exists() else None


def _find_first(root: Path, names: list[str]) -> Path | None:
    for name in names:
        hits = list(root.rglob(name))
        if hits:
            return hits[0]
    return None


def _find_many(root: Path, patterns: list[str]) -> list[Path]:
    found = []
    seen = set()
    for pat in patterns:
        for p in sorted(root.rglob(pat)):
            key = str(p.resolve()) if p.exists() else str(p)
            if key not in seen:
                seen.add(key)
                found.append(p)
    return found


def _file_download_link(path: Path, label: str | None = None):
    # Display the staged path for now; download callbacks are wired separately later.
    return html.Div(
        f"{label or path.name}: {path}",
        style={"fontSize": "12px", "color": COLORS["muted"], "marginTop": "6px"},
    )


def _render_csv(path: Path, title: str | None = None, max_rows: int = 250):
    try:
        import pandas as pd
        df = pd.read_csv(path)
        shown = df.head(max_rows)
        children = []
        if title:
            children.append(html.H4(title, style={"color": COLORS["navy"]}))
        children.append(
            html.Div(
                f"{len(df):,} rows × {len(df.columns):,} columns"
                + (f" — showing first {len(shown):,}" if len(df) > len(shown) else ""),
                style={"color": COLORS["muted"], "marginBottom": "8px"},
            )
        )
        children.append(
            dash_table.DataTable(
                data=shown.to_dict("records"),
                columns=[{"name": c, "id": c} for c in shown.columns],
                page_size=20,
                sort_action="native",
                filter_action="native",
                style_table={"overflowX": "auto", "maxHeight": "620px", "overflowY": "auto"},
                style_cell={
                    "textAlign": "left",
                    "fontSize": "12px",
                    "padding": "6px",
                    "minWidth": "90px",
                    "maxWidth": "300px",
                    "whiteSpace": "normal",
                },
                style_header={"fontWeight": "700", "background": "#eef4fb"},
            )
        )
        children.append(_file_download_link(path))
        return html.Div(children)
    except Exception as exc:
        return placeholder(f"Could not read {path.name}: {exc}")


def _render_json(path: Path, title: str | None = None):
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
        pretty = json.dumps(obj, indent=2)
    except Exception:
        pretty = path.read_text(encoding="utf-8", errors="replace")
    return html.Div([
        html.H4(title or path.name, style={"color": COLORS["navy"]}),
        html.Pre(
            pretty,
            style={
                "whiteSpace": "pre-wrap",
                "background": "#0f172a",
                "color": "#e5eef9",
                "padding": "12px",
                "borderRadius": "9px",
                "maxHeight": "650px",
                "overflowY": "auto",
                "fontSize": "12px",
            },
        ),
        _file_download_link(path),
    ])


def _render_markdown(path: Path, title: str | None = None):
    return html.Div([
        html.H4(title or path.name, style={"color": COLORS["navy"]}),
        dcc.Markdown(path.read_text(encoding="utf-8", errors="replace")),
        _file_download_link(path),
    ])


def _render_html_plot(path: Path, title: str | None = None):
    return html.Div([
        html.H4(title or path.stem, style={"color": COLORS["navy"]}),
        html.Iframe(
            srcDoc=path.read_text(encoding="utf-8", errors="replace"),
            style={"width": "100%", "height": "650px", "border": f"1px solid {COLORS['border']}", "borderRadius": "8px"},
        ),
        _file_download_link(path),
    ], style={"marginBottom": "18px"})


def _render_file_list(paths: list[Path], heading: str = "Files"):
    if not paths:
        return placeholder("No files found for this section.")
    return html.Div([
        html.H4(heading, style={"color": COLORS["navy"]}),
        html.Ul([html.Li(str(p)) for p in paths], style={"fontSize": "13px"}),
    ])


def _stack(children):
    kids = [c for c in children if c is not None]
    return html.Div(kids) if kids else placeholder("No matching result files were found for this run.")


def _stage_root(store):
    root = _run_result_root(store)
    if root:
        return root
    return None


def _section_dir(root: Path, folder: str) -> Path:
    p = root / folder
    return p if p.exists() else root


def _first_existing(root: Path, names: list[str]) -> Path | None:
    return _find_first(root, names)


def _ordered_rna_omics_plots(root: Path) -> list[Path]:
    """Return RNA omics HTML plots in a stable, user-friendly order."""
    priority = [
        "pca_plot.html",
        "umap_plot.html",
        "volcano_plot.html",
        "heatmap_plot.html",
    ]
    chosen: list[Path] = []
    seen: set[str] = set()

    for name in priority:
        p = _first_existing(root, [name])
        if p:
            key = str(p.resolve()) if p.exists() else str(p)
            if key not in seen:
                chosen.append(p)
                seen.add(key)

    for p in _find_many(root, ["*.html"]):
        key = str(p.resolve()) if p.exists() else str(p)
        if key not in seen:
            chosen.append(p)
            seen.add(key)

    return chosen


# -----------------------------------------------------------------------------
# RNA Outputs
# -----------------------------------------------------------------------------
def rna_outputs_layout():
    return html.Div(
        [
            section_title("RNA Outputs", "Bulk RNA-seq integrated analysis results."),
            dcc.Tabs(
                id="rna-family-tabs",
                value="rna-omics",
                children=[
                    dcc.Tab(label="Omics", value="rna-omics"),
                    dcc.Tab(label="Network", value="rna-network"),
                    dcc.Tab(label="Triage", value="rna-triage"),
                    dcc.Tab(label="Literature", value="rna-literature"),
                    dcc.Tab(label="Final Report", value="rna-final"),
                ],
            ),
            html.Div(id="rna-family-panel", style={"marginTop": "14px"}),
        ]
    )


@app.callback(Output("rna-family-panel", "children"), Input("rna-family-tabs", "value"))
def render_rna_family(tab):
    if tab == "rna-network":
        return output_panel(
            "rna-network",
            "Network",
            nested_output_tabs(
                "rna-network",
                [
                    ("Main Network", "main"),
                    ("Bipartite Diffusion", "bipartite"),
                    ("Pathway Projection", "projection"),
                    ("Projection Diffusion", "projection-diffusion"),
                    ("Consensus", "consensus"),
                    ("Tables / Files", "files"),
                ],
                "main",
            ),
        )
    if tab == "rna-triage":
        return output_panel(
            "rna-triage",
            "Evidence-Aware Triage",
            nested_output_tabs(
                "rna-triage",
                [
                    ("Overview", "overview"),
                    ("Triage Rows", "rows"),
                    ("Programs", "programs"),
                    ("Claims", "claims"),
                    ("LLM Report", "report"),
                    ("PubMed Context", "pubmed"),
                    ("Playbook Guidance", "playbook"),
                    ("Validation / Follow-up", "validation"),
                    ("Files", "files"),
                ],
                "overview",
            ),
        )
    if tab == "rna-literature":
        return output_panel(
            "rna-literature",
            "Literature Review",
            nested_output_tabs(
                "rna-literature",
                [
                    ("Evidence Summary", "summary"),
                    ("Findings", "findings"),
                    ("Assessments", "assessments"),
                    ("Direct / Mechanistic", "strong"),
                    ("Contextual / Background", "context"),
                    ("Excluded", "excluded"),
                    ("Mechanistic Paths", "paths"),
                    ("Raw JSON", "json"),
                ],
                "summary",
            ),
        )
    if tab == "rna-final":
        return output_panel(
            "rna-final",
            "Final Report",
            nested_output_tabs(
                "rna-final",
                [("Rendered Report", "report"), ("Report JSON", "json"), ("Files / Export", "files")],
                "report",
            ),
        )
    panel = output_panel(
        "rna-omics",
        "Omics",
        nested_output_tabs(
            "rna-omics",
            [
                ("Overview", "overview"),
                ("Differential Expression", "de"),
                ("Expression / QC", "qc"),
                ("Enrichment", "enrichment"),
                ("Plots", "plots"),
                ("Classifier", "classifier"),
                ("Power", "power"),
                ("Objects / Reproducibility", "objects"),
            ],
            "overview",
        ),
    )
    panel.children.insert(
        3,
        html.Div(
            [
                html.Button("← Previous plot", id="rna-omics-plot-prev", n_clicks=0, className="export-btn"),
                html.Button("Next plot →", id="rna-omics-plot-next", n_clicks=0, className="export-btn"),
            ],
            id="rna-omics-plot-nav",
            style={"display": "none", "gap": "8px", "marginBottom": "12px"},
        ),
    )
    return panel



# -----------------------------------------------------------------------------
# RNA result-population callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("rna-omics-inner-panel", "children"),
    Input("rna-omics-inner-tabs", "value"),
    Input("active-run-store", "data"),
    Input("rna-omics-plot-index", "data"),
)
def populate_rna_omics(tab, store, plot_index):
    root = _stage_root(store)
    if not root:
        return placeholder("Run or load a completed RNA-seq pipeline result to populate this section.")

    d = _section_dir(root, "01_rnaseq")

    if tab == "overview":
        files = []
        for name in ["run_summary.json", "run_status.json", "manifest.json"]:
            p = _first_existing(d, [name])
            if p:
                files.append(_render_json(p, name))
        return _stack(files)

    if tab == "de":
        candidates = _find_many(d, ["DE_results_*.csv", "limma_all_results_*.csv", "DE_results.csv", "limma_all_results.csv"])
        return _stack([_render_csv(p, p.name) for p in candidates])

    if tab == "qc":
        candidates = _find_many(d, ["normalized_logCPM_*.csv", "sample_qc_*.csv", "normalized_logCPM.csv", "sample_qc.csv"])
        return _stack([_render_csv(p, p.name) for p in candidates])

    if tab == "enrichment":
        candidates = _find_many(d, ["enrichment_all.csv", "enrichment_upregulated.csv", "enrichment_downregulated.csv", "enrichment_up.csv", "enrichment_down.csv"])
        plots = _find_many(d, ["*enrichment*.html"])
        return _stack([_render_csv(p, p.name) for p in candidates] + [_render_html_plot(p, p.name) for p in plots])

    if tab == "plots":
        plots = _ordered_rna_omics_plots(d)
        if not plots:
            return placeholder("No RNA omics HTML plots were found for this run.")

        idx = int(plot_index or 0)
        idx = idx % len(plots)
        current = plots[idx]

        return html.Div(
            [
                html.Div(
                    f"{idx + 1} of {len(plots)} — {current.name}",
                    style={
                        "fontWeight": "700",
                        "color": COLORS["navy"],
                        "marginBottom": "12px",
                    },
                ),
                _render_html_plot(current, current.name),
            ]
        )

    if tab == "classifier":
        candidates = _find_many(d, ["rf_*.csv", "classifier_*.csv"])
        roc = _find_many(d, ["rf_roc_plot.html"])
        return _stack([_render_csv(p, p.name) for p in candidates] + [_render_html_plot(p, p.name) for p in roc])

    if tab == "power":
        candidates = _find_many(d, ["*power*.csv"])
        plots = _find_many(d, ["*power*.html"])
        return _stack([_render_csv(p, p.name) for p in candidates] + [_render_html_plot(p, p.name) for p in plots])

    files = _find_many(d, ["*.rds", "sessionInfo.txt", "run_log.txt", "manifest.json", "run_status.json", "run_summary.json"])
    return _render_file_list(files, "RNA objects and reproducibility files")


def _render_network_json_summary(path: Path, title: str):
    """Render only lightweight network metadata.

    Network API JSON files can contain the entire graph and Plotly figure payload.
    Dumping those structures into html.Pre can create megabytes of DOM text and
    freeze the browser. This renderer intentionally extracts only small summaries.
    """
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return placeholder(f"Could not read {path.name}: {exc}")

    summary = {}
    for key in (
        "message",
        "preset_message",
        "mapped_columns",
        "input_shape",
        "stats",
        "top_groups",
        "top_items",
        "seed_node",
        "alpha",
        "ranking_mode",
        "warning",
    ):
        if key in obj:
            value = obj.get(key)
            # Keep summary lists compact.
            if isinstance(value, list):
                value = value[:20]
            summary[key] = value

    children = [
        html.H4(title, style={"color": COLORS["navy"]}),
        html.Pre(
            json.dumps(summary, indent=2, default=str),
            style={
                "whiteSpace": "pre-wrap",
                "background": "#0f172a",
                "color": "#e5eef9",
                "padding": "12px",
                "borderRadius": "9px",
                "maxHeight": "360px",
                "overflowY": "auto",
                "fontSize": "12px",
            },
        ),
    ]

    # Do not instantiate dcc.Graph here. On this Dash installation the Plotly
    # component's lazy-loaded JavaScript chunk times out when the Network tab is
    # first opened. The CSV result tables and compact JSON summary remain fully
    # usable; graph visualization can be reintroduced later as a prebuilt HTML
    # iframe without relying on Dash's async Plotly bundle.
    if isinstance(obj.get("figure"), dict):
        children.append(
            html.Div(
                "Interactive network figure is available in the API output, but "
                "browser rendering is temporarily disabled to keep this dashboard stable.",
                style={
                    "padding": "10px",
                    "background": COLORS["warn"],
                    "borderRadius": "8px",
                    "color": COLORS["navy"],
                    "marginBottom": "8px",
                },
            )
        )

    children.append(_file_download_link(path))
    return html.Div(children, style={"marginBottom": "16px"})


@app.callback(
    Output("rna-network-inner-panel", "children"),
    Input("rna-network-inner-tabs", "value"),
    Input("active-run-store", "data"),
)
def populate_rna_network(tab, store):
    root = _stage_root(store)
    if not root:
        return placeholder("Run or load a completed RNA-seq pipeline result to populate this section.")
    d = _section_dir(root, "02_network")

    mapping = {
        "main": (["main_nodes.csv", "main_edges.csv"], "network_build.json"),
        "bipartite": (["bipartite_diffusion_results.csv", "bipartite_top_candidates.csv"], "bipartite_diffusion.json"),
        "projection": (["projection_nodes.csv", "projection_edges.csv"], "projection_build.json"),
        "projection-diffusion": (["projection_diffusion_results.csv", "projection_top_candidates.csv"], "projection_diffusion.json"),
        "consensus": (["consensus_candidates.csv"], "consensus.json"),
    }

    if tab == "files":
        files = [p for p in _find_many(d, ["*"]) if p.is_file()]
        return _render_file_list(files[:100], "Network output files")

    csv_names, json_name = mapping.get(tab, ([], None))
    parts = []

    for name in csv_names:
        p = _first_existing(d, [name])
        if p:
            parts.append(_render_csv(p, name, max_rows=150))

    if json_name:
        p = _first_existing(d, [json_name])
        if p:
            parts.append(_render_network_json_summary(p, json_name))

    return _stack(parts)


@app.callback(
    Output("rna-triage-inner-panel", "children"),
    Input("rna-triage-inner-tabs", "value"),
    Input("active-run-store", "data"),
)
def populate_rna_triage(tab, store):
    root = _stage_root(store)
    if not root:
        return placeholder("Run or load a completed RNA-seq pipeline result to populate this section.")
    d = _section_dir(root, "03_llm_triage")

    if tab == "overview":
        p = _first_existing(d, ["analysis_result.json"])
        return _render_json(p, "Triage analysis result") if p else _render_file_list(_find_many(d, ["*.json"]), "Triage JSON files")
    if tab == "rows":
        p = _first_existing(d, ["triage_rows.csv"])
        return _render_csv(p, "Triage Rows") if p else placeholder("triage_rows.csv was not found.")
    if tab == "programs":
        p = _first_existing(d, ["programs.csv"])
        return _render_csv(p, "Programs") if p else placeholder("programs.csv was not found.")
    if tab == "claims":
        p = _first_existing(d, ["claims.csv"])
        return _render_csv(p, "Claims") if p else placeholder("claims.csv was not found.")
    if tab == "report":
        p = _first_existing(d, ["report.md"])
        return _render_markdown(p, "LLM Triage Report") if p else placeholder("report.md was not found.")
    if tab == "pubmed":
        p = _first_existing(d, ["analysis_result.json"])
        return _render_json(p, "PubMed / literature context") if p else placeholder("analysis_result.json was not found.")
    if tab == "playbook":
        p = _first_existing(d, ["analysis_result.json"])
        return _render_json(p, "Playbook retrieval context") if p else placeholder("analysis_result.json was not found.")
    if tab == "validation":
        p = _first_existing(d, ["claims.csv"])
        return _render_csv(p, "Claims and validation steps") if p else placeholder("claims.csv was not found.")
    return _render_file_list(_find_many(d, ["*"]), "Triage output files")


@app.callback(
    Output("rna-literature-inner-panel", "children"),
    Input("rna-literature-inner-tabs", "value"),
    Input("active-run-store", "data"),
)
def populate_rna_literature(tab, store):
    root = _stage_root(store)
    if not root:
        return placeholder("Run or load a completed RNA-seq pipeline result to populate this section.")
    d = _section_dir(root, "04_literature_review")
    req = _first_existing(d, ["literature_review_request.json"])
    resp = _first_existing(d, ["literature_review_response.json"])

    if tab == "json":
        return _stack([
            _render_json(req, "Literature Review Request") if req else None,
            _render_json(resp, "Literature Review Response") if resp else None,
        ])
    if not resp:
        return placeholder("literature_review_response.json was not found.")

    try:
        obj = json.loads(resp.read_text(encoding="utf-8"))
        assessments = obj.get("assessments", [])
        if tab == "summary":
            return html.Div([
                html.H4("Evidence Summary", style={"color": COLORS["navy"]}),
                html.Pre(json.dumps(obj.get("summary_counts", {}), indent=2)),
                html.Div(f"Assessments: {len(assessments):,}"),
            ])
        if tab == "findings":
            if req:
                return _render_json(req, "Literature findings / handoff")
            return placeholder("Literature request JSON was not found.")

        wanted = assessments
        if tab == "strong":
            wanted = [x for x in assessments if x.get("level") in {"direct", "mechanistic_support"}]
        elif tab == "context":
            wanted = [x for x in assessments if x.get("level") in {"contextual_support", "background"}]
        elif tab == "excluded":
            wanted = [x for x in assessments if x.get("level") == "exclude"]
        elif tab == "paths":
            wanted = [x for x in assessments if x.get("playbook_match")]

        import pandas as pd
        rows = []
        for a in wanted:
            paper = a.get("paper") or {}
            rows.append({
                "finding_id": a.get("finding_id"),
                "level": a.get("level"),
                "relevance_score": a.get("relevance_score"),
                "title": paper.get("title"),
                "year": paper.get("year"),
                "journal": paper.get("journal"),
                "pmid": paper.get("pmid"),
                "doi": paper.get("doi"),
                "assay_match": a.get("assay_match"),
                "context_match": a.get("context_match"),
                "playbook_distance": a.get("playbook_distance"),
                "reasons": " | ".join(a.get("reasons") or []),
                "caution": a.get("caution"),
            })
        if not rows:
            return placeholder("No assessments matched this evidence view.")
        temp = d / f"_dashboard_{tab}.csv"
        pd.DataFrame(rows).to_csv(temp, index=False)
        return _render_csv(temp, f"Literature — {tab}")
    except Exception as exc:
        return placeholder(f"Could not render literature response: {exc}")


@app.callback(
    Output("rna-final-inner-panel", "children"),
    Input("rna-final-inner-tabs", "value"),
    Input("active-run-store", "data"),
)
def populate_rna_final(tab, store):
    root = _stage_root(store)
    if not root:
        return placeholder("Run or load a completed RNA-seq pipeline result to populate this section.")

    md = _first_existing(root, ["rnaseq_final_report.md"])
    js = _first_existing(root, ["rnaseq_final_report.json"])

    if tab == "report":
        return _render_markdown(md, "Integrated RNA-seq Final Report") if md else placeholder("rnaseq_final_report.md was not found.")
    if tab == "json":
        return _render_json(js, "Final Report JSON") if js else placeholder("rnaseq_final_report.json was not found.")
    return _render_file_list(_find_many(root, ["*"]), "Integrated RNA pipeline files")



# -----------------------------------------------------------------------------
# RNA omics plot navigation
# -----------------------------------------------------------------------------
@app.callback(
    Output("rna-omics-plot-nav", "style"),
    Input("rna-omics-inner-tabs", "value"),
)
def toggle_rna_omics_plot_nav(tab):
    if tab == "plots":
        return {
            "display": "flex",
            "gap": "8px",
            "marginBottom": "12px",
            "flexWrap": "wrap",
        }
    return {"display": "none"}


@app.callback(
    Output("rna-omics-plot-index", "data"),
    Input("rna-omics-plot-prev", "n_clicks"),
    Input("rna-omics-plot-next", "n_clicks"),
    Input("rna-omics-inner-tabs", "value"),
    Input("active-run-store", "data"),
    State("rna-omics-plot-index", "data"),
    prevent_initial_call=True,
)
def update_rna_omics_plot_index(_prev, _next, tab, store, current_index):
    if tab != "plots":
        return 0

    root = _run_result_root(store)
    if not root:
        return 0

    d = _section_dir(root, "01_rnaseq")
    plots = _ordered_rna_omics_plots(d)
    if not plots:
        return 0

    triggered = ctx.triggered_id
    if triggered in ("rna-omics-inner-tabs", "active-run-store"):
        return 0

    idx = int(current_index or 0) % len(plots)
    if triggered == "rna-omics-plot-prev":
        return (idx - 1) % len(plots)
    if triggered == "rna-omics-plot-next":
        return (idx + 1) % len(plots)
    return idx

# -----------------------------------------------------------------------------
# RNA download callbacks
# -----------------------------------------------------------------------------
def _make_zip(source_dir: Path, zip_path: Path) -> Path:
    """Create/replace a ZIP containing one result section."""
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    if zip_path.exists():
        zip_path.unlink()
    base = str(zip_path.with_suffix(""))
    made = Path(shutil.make_archive(base, "zip", root_dir=str(source_dir)))
    return made


def _full_bundle_path(store: dict | None) -> Path | None:
    if not store:
        return None
    run_dir = Path(store.get("run_dir", ""))
    candidate = run_dir / "results" / "jcap_integrated_results.zip"
    if candidate.exists():
        return candidate

    root = _run_result_root(store)
    if root:
        # Some Nextflow bundles may not publish a sibling zip. Build one from the
        # staged result directory so the dashboard button always has a real file.
        return _make_zip(root, run_dir / "downloads" / "jcap_integrated_results.zip")
    return None


def _rna_current_file(prefix: str, tab: str, root: Path, kind: str, plot_index: int | None = None) -> Path | None:
    """Resolve the file represented by the currently visible RNA subtab."""
    if prefix == "rna-omics":
        d = _section_dir(root, "01_rnaseq")
        if kind == "table":
            mapping = {
                "de": ["DE_results_*.csv", "limma_all_results_*.csv"],
                "qc": ["sample_qc_*.csv", "normalized_logCPM_*.csv"],
                "enrichment": ["enrichment_all.csv", "enrichment_upregulated.csv", "enrichment_up.csv"],
                "classifier": ["rf_metrics_*.csv", "classifier_metrics*.csv", "rf_predictions_*.csv"],
                "power": ["*power*.csv"],
            }
        else:
            if tab == "plots":
                plots = _ordered_rna_omics_plots(d)
                if not plots:
                    return None
                idx = int(plot_index or 0) % len(plots)
                return plots[idx]

            mapping = {
                "enrichment": ["*enrichment*.html"],
                "classifier": ["rf_roc_plot.html"],
                "power": ["*power*.html"],
            }
    elif prefix == "rna-network":
        d = _section_dir(root, "02_network")
        mapping = {
            "main": ["main_nodes.csv", "main_edges.csv"],
            "bipartite": ["bipartite_diffusion_results.csv", "bipartite_top_candidates.csv"],
            "projection": ["projection_nodes.csv", "projection_edges.csv"],
            "projection-diffusion": ["projection_diffusion_results.csv", "projection_top_candidates.csv"],
            "consensus": ["consensus_candidates.csv"],
        } if kind == "table" else {}
    elif prefix == "rna-triage":
        d = _section_dir(root, "03_llm_triage")
        mapping = {
            "rows": ["triage_rows.csv"],
            "programs": ["programs.csv"],
            "claims": ["claims.csv"],
            "validation": ["claims.csv"],
        } if kind == "table" else {}
    elif prefix == "rna-literature":
        d = _section_dir(root, "04_literature_review")
        if kind == "table":
            mapping = {
                "assessments": ["_dashboard_assessments.csv"],
                "strong": ["_dashboard_strong.csv"],
                "context": ["_dashboard_context.csv"],
                "excluded": ["_dashboard_excluded.csv"],
                "paths": ["_dashboard_paths.csv"],
            }
        else:
            mapping = {}
    else:  # rna-final
        d = root
        mapping = {}

    patterns = mapping.get(tab, [])
    for pat in patterns:
        hits = sorted(d.rglob(pat))
        if hits:
            return hits[0]
    return None


def _rna_section_dir(prefix: str, root: Path) -> Path:
    return {
        "rna-omics": _section_dir(root, "01_rnaseq"),
        "rna-network": _section_dir(root, "02_network"),
        "rna-triage": _section_dir(root, "03_llm_triage"),
        "rna-literature": _section_dir(root, "04_literature_review"),
        "rna-final": root,
    }[prefix]


def _rna_download(prefix: str, tab: str, store: dict | None, plot_index: int | None = None):
    triggered = ctx.triggered_id
    if not triggered:
        return no_update

    root = _run_result_root(store)
    if not root:
        return no_update

    if triggered.endswith("-download-table"):
        path = _rna_current_file(prefix, tab, root, "table", plot_index)
        return dcc.send_file(str(path)) if path and path.exists() else no_update

    if triggered.endswith("-download-plot"):
        path = _rna_current_file(prefix, tab, root, "plot", plot_index)
        return dcc.send_file(str(path)) if path and path.exists() else no_update

    if triggered.endswith("-download-section"):
        section = _rna_section_dir(prefix, root)
        run_dir = Path(store.get("run_dir", root.parent))
        zip_path = run_dir / "downloads" / f"{prefix}_section.zip"
        return dcc.send_file(str(_make_zip(section, zip_path)))

    if triggered.endswith("-download-bundle"):
        path = _full_bundle_path(store)
        return dcc.send_file(str(path)) if path and path.exists() else no_update

    return no_update


def _register_rna_download_callback(prefix: str):
    @app.callback(
        Output(f"{prefix}-download", "data"),
        Input(f"{prefix}-download-table", "n_clicks"),
        Input(f"{prefix}-download-plot", "n_clicks"),
        Input(f"{prefix}-download-section", "n_clicks"),
        Input(f"{prefix}-download-bundle", "n_clicks"),
        State(f"{prefix}-inner-tabs", "value"),
        State("active-run-store", "data"),
        State("rna-omics-plot-index", "data"),
        prevent_initial_call=True,
    )
    def _download(_n1, _n2, _n3, _n4, tab, store, plot_index):
        return _rna_download(prefix, tab, store, plot_index)


for _prefix in ["rna-omics", "rna-network", "rna-triage", "rna-literature", "rna-final"]:
    _register_rna_download_callback(_prefix)


# -----------------------------------------------------------------------------
# scRNA Outputs
# -----------------------------------------------------------------------------
def scrna_outputs_layout():
    return html.Div(
        [
            section_title("scRNA Outputs", "Single-cell RNA-seq integrated analysis results."),
            dcc.Tabs(
                id="scrna-family-tabs",
                value="scrna-omics",
                children=[
                    dcc.Tab(label="Omics", value="scrna-omics"),
                    dcc.Tab(label="Cell State", value="scrna-cellstate"),
                    dcc.Tab(label="Network", value="scrna-network"),
                    dcc.Tab(label="Triage", value="scrna-triage"),
                    dcc.Tab(label="Literature", value="scrna-literature"),
                    dcc.Tab(label="Final Report", value="scrna-final"),
                ],
            ),
            html.Div(id="scrna-family-panel", style={"marginTop": "14px"}),
        ]
    )


@app.callback(Output("scrna-family-panel", "children"), Input("scrna-family-tabs", "value"))
def render_scrna_family(tab):
    if tab == "scrna-cellstate":
        return output_panel(
            "scrna-cellstate",
            "Cell State",
            nested_output_tabs(
                "scrna-cellstate",
                [
                    ("Overview", "overview"),
                    ("Overall Frequencies", "overall"),
                    ("State by Condition", "by-condition"),
                    ("Annotated Metadata", "metadata"),
                    ("Handoff JSON", "handoff"),
                ],
                "overview",
            ),
        )
    if tab == "scrna-network":
        return output_panel(
            "scrna-network",
            "Network",
            nested_output_tabs(
                "scrna-network",
                [
                    ("Main Network", "main"),
                    ("Bipartite Diffusion", "bipartite"),
                    ("Pathway Projection", "projection"),
                    ("Projection Diffusion", "projection-diffusion"),
                    ("Consensus", "consensus"),
                    ("Tables / Files", "files"),
                ],
                "main",
            ),
        )
    if tab == "scrna-triage":
        return output_panel(
            "scrna-triage",
            "Evidence-Aware Triage",
            nested_output_tabs(
                "scrna-triage",
                [
                    ("Overview", "overview"),
                    ("Triage Rows", "rows"),
                    ("Programs", "programs"),
                    ("Claims", "claims"),
                    ("LLM Report", "report"),
                    ("PubMed Context", "pubmed"),
                    ("Playbook Guidance", "playbook"),
                    ("Validation / Follow-up", "validation"),
                    ("Files", "files"),
                ],
                "overview",
            ),
        )
    if tab == "scrna-literature":
        return output_panel(
            "scrna-literature",
            "Literature Review",
            nested_output_tabs(
                "scrna-literature",
                [
                    ("Evidence Summary", "summary"),
                    ("Findings", "findings"),
                    ("Assessments", "assessments"),
                    ("Direct / Mechanistic", "strong"),
                    ("Contextual / Background", "context"),
                    ("Excluded", "excluded"),
                    ("Mechanistic Paths", "paths"),
                    ("Raw JSON", "json"),
                ],
                "summary",
            ),
        )
    if tab == "scrna-final":
        return output_panel(
            "scrna-final",
            "Final Report",
            nested_output_tabs(
                "scrna-final",
                [("Rendered Report", "report"), ("Report JSON", "json"), ("Files / Export", "files")],
                "report",
            ),
        )
    return output_panel(
        "scrna-omics",
        "Omics",
        nested_output_tabs(
            "scrna-omics",
            [
                ("Overview", "overview"),
                ("Differential Expression", "de"),
                ("PCA / UMAP", "dimred"),
                ("Enrichment", "enrichment"),
                ("Classifier", "classifier"),
                ("Power", "power"),
                ("Plots", "plots"),
                ("Objects / Reproducibility", "objects"),
            ],
            "overview",
        ),
    )


# -----------------------------------------------------------------------------
# Perturb-seq Outputs
# -----------------------------------------------------------------------------
def perturb_outputs_layout():
    return html.Div(
        [
            section_title("Perturb-seq Outputs", "Perturb-seq / Mixscape integrated analysis results."),
            dcc.Tabs(
                id="perturb-family-tabs",
                value="perturb-omics",
                children=[
                    dcc.Tab(label="Mixscape / Omics", value="perturb-omics"),
                    dcc.Tab(label="Cell State", value="perturb-cellstate"),
                    dcc.Tab(label="Network", value="perturb-network"),
                    dcc.Tab(label="Triage", value="perturb-triage"),
                    dcc.Tab(label="Literature", value="perturb-literature"),
                    dcc.Tab(label="Final Report", value="perturb-final"),
                ],
            ),
            html.Div(id="perturb-family-panel", style={"marginTop": "14px"}),
        ]
    )


@app.callback(Output("perturb-family-panel", "children"), Input("perturb-family-tabs", "value"))
def render_perturb_family(tab):
    if tab == "perturb-cellstate":
        return output_panel(
            "perturb-cellstate",
            "Cell State",
            nested_output_tabs(
                "perturb-cellstate",
                [
                    ("Overview", "overview"),
                    ("Overall Frequencies", "overall"),
                    ("State by Perturbation", "by-perturbation"),
                    ("Annotated Metadata", "metadata"),
                    ("Handoff JSON", "handoff"),
                ],
                "overview",
            ),
        )
    if tab == "perturb-network":
        return output_panel(
            "perturb-network",
            "Network",
            nested_output_tabs(
                "perturb-network",
                [
                    ("Main Network", "main"),
                    ("Bipartite Diffusion", "bipartite"),
                    ("Pathway Projection", "projection"),
                    ("Projection Diffusion", "projection-diffusion"),
                    ("Consensus", "consensus"),
                    ("Tables / Files", "files"),
                ],
                "main",
            ),
        )
    if tab == "perturb-triage":
        return output_panel(
            "perturb-triage",
            "Evidence-Aware Triage",
            nested_output_tabs(
                "perturb-triage",
                [
                    ("Overview", "overview"),
                    ("Triage Rows", "rows"),
                    ("Programs", "programs"),
                    ("Claims", "claims"),
                    ("LLM Report", "report"),
                    ("PubMed Context", "pubmed"),
                    ("Playbook Guidance", "playbook"),
                    ("Validation / Follow-up", "validation"),
                    ("Files", "files"),
                ],
                "overview",
            ),
        )
    if tab == "perturb-literature":
        return output_panel(
            "perturb-literature",
            "Literature Review",
            nested_output_tabs(
                "perturb-literature",
                [
                    ("Evidence Summary", "summary"),
                    ("Findings", "findings"),
                    ("Assessments", "assessments"),
                    ("Direct / Mechanistic", "strong"),
                    ("Contextual / Background", "context"),
                    ("Excluded", "excluded"),
                    ("Mechanistic Paths", "paths"),
                    ("Raw JSON", "json"),
                ],
                "summary",
            ),
        )
    if tab == "perturb-final":
        return output_panel(
            "perturb-final",
            "Final Report",
            nested_output_tabs(
                "perturb-final",
                [("Rendered Report", "report"), ("Report JSON", "json"), ("Files / Export", "files")],
                "report",
            ),
        )
    return output_panel(
        "perturb-omics",
        "Mixscape / Omics",
        nested_output_tabs(
            "perturb-omics",
            [
                ("Overview", "overview"),
                ("Mixscape Classification", "mixscape"),
                ("Perturbation DE", "de"),
                ("Summary Statistics", "summary"),
                ("Enrichment", "enrichment"),
                ("Plots", "plots"),
                ("Statistical Notes", "notes"),
                ("Objects / Reproducibility", "objects"),
            ],
            "overview",
        ),
    )


# -----------------------------------------------------------------------------
# Main page
# -----------------------------------------------------------------------------
app.layout = html.Div(
    [
        dcc.Markdown(
            """
            <style>
            body { margin: 0; background: #f4f7fb; font-family: Inter, Segoe UI, Arial, sans-serif; color: #172033; }
            .export-btn {
                background: white;
                border: 1px solid #d7e0eb;
                border-radius: 8px;
                padding: 8px 11px;
                font-weight: 650;
                color: #173b67;
                cursor: pointer;
            }
            .export-btn:hover { background: #f2f6fb; }
            </style>
            """,
            dangerously_allow_html=True,
        ),
        html.Div(
            [
                html.H1("JCAP Bioinformatics Platform", style={"margin": 0, "fontSize": "30px"}),
                html.Div(
                    "One workspace for bulk RNA-seq, scRNA-seq, Perturb-seq, network analysis, evidence-aware triage, literature review, and final reporting.",
                    style={"marginTop": "7px", "opacity": 0.9},
                ),
            ],
            style={
                "background": COLORS["navy"],
                "color": "white",
                "padding": "24px 30px",
                "boxShadow": "0 4px 16px rgba(0,0,0,0.12)",
            },
        ),
        html.Div(
            [
                get_started_layout(),
                card([rna_outputs_layout()]),
                card([scrna_outputs_layout()]),
                card([perturb_outputs_layout()]),
            ],
            style={"maxWidth": "1500px", "margin": "22px auto", "padding": "0 18px 50px 18px"},
        ),
    ]
)


# -----------------------------------------------------------------------------
# Integrated RNA-seq execution
# -----------------------------------------------------------------------------
@app.callback(
    Output("job-status", "children"),
    Output("active-run-store", "data"),
    Input("run-integrated-btn", "n_clicks"),
    State("pipeline-picker", "value"),
    State("run-mode", "value"),
    State("run-counts-upload", "contents"),
    State("run-counts-upload", "filename"),
    State("run-metadata-upload", "contents"),
    State("run-metadata-upload", "filename"),
    prevent_initial_call=True,
)
def launch_integrated_pipeline(n_clicks, pipeline, run_mode, counts_contents, counts_filename, metadata_contents, metadata_filename):
    if pipeline != "rnaseq":
        return "RNA-seq execution is wired first. scRNA-seq and Perturb-seq launch wiring come next.", no_update

    if not RNASEQ_PIPELINE_DIR.exists():
        return f"RNA pipeline directory was not found: {RNASEQ_PIPELINE_DIR}", no_update
    if not RNASEQ_PARAMS_FILE.exists():
        return f"RNA params file was not found: {RNASEQ_PARAMS_FILE}", no_update

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_id = f"rnaseq_{stamp}"
    run_dir = DASHBOARD_RUN_ROOT / run_id
    input_dir = run_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    try:
        if run_mode == "test":
            if not RNASEQ_TEST_COUNTS.exists() or not RNASEQ_TEST_PHENOTYPE.exists():
                return (
                    "Built-in RNA test files were not found. "
                    f"Counts={RNASEQ_TEST_COUNTS}; phenotype={RNASEQ_TEST_PHENOTYPE}",
                    no_update,
                )
            counts_path = RNASEQ_TEST_COUNTS
            phenotype_path = RNASEQ_TEST_PHENOTYPE
        else:
            if not counts_contents or not metadata_contents:
                return "Upload both the counts CSV and phenotype CSV before running RNA-seq.", no_update
            counts_name = Path(counts_filename or "counts_data.csv").name
            metadata_name = Path(metadata_filename or "phenotype_data.csv").name
            counts_path = decode_dash_upload(counts_contents, input_dir / counts_name)
            phenotype_path = decode_dash_upload(metadata_contents, input_dir / metadata_name)

        log_path = run_dir / "nextflow.log"
        cmd = [
            "nextflow", "run", "main.nf",
            "-params-file", str(RNASEQ_PARAMS_FILE),
            "--counts", str(counts_path),
            "--phenotype", str(phenotype_path),
        ]

        log_handle = log_path.open("w", encoding="utf-8")
        process = subprocess.Popen(
            cmd, cwd=str(RNASEQ_PIPELINE_DIR),
            stdout=log_handle, stderr=subprocess.STDOUT, text=True
        )
        log_handle.close()
        ACTIVE_PROCESSES[run_id] = process

        meta = {
            "run_id": run_id, "pipeline": "rnaseq", "run_mode": run_mode,
            "status": "running", "pid": process.pid,
            "started_at": utc_now_iso(), "finished_at": None,
            "return_code": None, "run_dir": str(run_dir),
            "log_path": str(log_path), "counts": str(counts_path),
            "phenotype": str(phenotype_path), "result_dir": None,
        }
        write_run_meta(run_dir, meta)
        return html.Div([
            html.Strong("RNA-seq pipeline launched."),
            html.Div(f"Run ID: {run_id}"),
            html.Div(f"PID: {process.pid}"),
            html.Div("Use Refresh Job Status to follow the Nextflow run."),
        ]), meta
    except Exception as exc:
        return f"Could not launch RNA-seq pipeline: {exc}", no_update

@app.callback(
    Output("job-status", "children", allow_duplicate=True),
    Output("active-run-store", "data", allow_duplicate=True),
    Input("refresh-job-btn", "n_clicks"),
    State("active-run-store", "data"),
    prevent_initial_call=True,
)
def refresh_integrated_job(n_clicks, store):
    if not store:
        return "No dashboard pipeline run has been started yet.", no_update

    run_id = store.get("run_id")
    run_dir = Path(store.get("run_dir", ""))
    log_path = Path(store.get("log_path", ""))
    if not run_id or not run_dir.exists():
        return "The saved dashboard run directory is no longer available.", no_update

    meta = read_run_meta(run_dir) or dict(store)
    process = ACTIVE_PROCESSES.get(run_id)
    if process is not None:
        rc = process.poll()
        if rc is None:
            meta["status"] = "running"
        else:
            meta["return_code"] = rc
            meta["finished_at"] = meta.get("finished_at") or utc_now_iso()
            if rc == 0:
                meta["status"] = "complete"
                staged = stage_completed_rnaseq_results(run_dir)
                if staged:
                    meta["result_dir"] = str(staged)
            else:
                meta["status"] = "failed"
            ACTIVE_PROCESSES.pop(run_id, None)
            write_run_meta(run_dir, meta)

    log_tail = ""
    if log_path.exists():
        try:
            lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
            log_tail = "\n".join(lines[-16:])
        except Exception:
            pass

    status_name = str(meta.get("status", "unknown")).upper()
    children = [html.Strong(f"{meta.get('pipeline', 'pipeline')} — {status_name}"), html.Div(f"Run ID: {run_id}")]
    if meta.get("return_code") is not None:
        children.append(html.Div(f"Nextflow return code: {meta['return_code']}"))
    if meta.get("result_dir"):
        children.append(html.Div(f"Results staged at: {meta['result_dir']}"))
    if log_tail:
        children.extend([
            html.Div("Latest Nextflow log:", style={"fontWeight": "700", "marginTop": "8px"}),
            html.Pre(log_tail, style={
                "whiteSpace": "pre-wrap", "maxHeight": "260px", "overflowY": "auto",
                "background": "#0f172a", "color": "#e5eef9", "padding": "10px",
                "borderRadius": "8px", "fontSize": "12px",
            }),
        ])
    return html.Div(children), meta


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)

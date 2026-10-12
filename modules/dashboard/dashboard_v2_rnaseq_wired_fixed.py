from __future__ import annotations

from dash import Dash, Input, Output, State, dcc, html, dash_table, no_update
from pathlib import Path
import base64
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
    if destination.exists():
        shutil.rmtree(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, destination)
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
                )
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
    return output_panel(
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

from __future__ import annotations

from dash import Dash, Input, Output, State, dcc, html, dash_table, no_update
from pathlib import Path


# -----------------------------------------------------------------------------
# App
# -----------------------------------------------------------------------------
app = Dash(__name__, suppress_callback_exceptions=True, title="JCAP Bioinformatics Platform")
server = app.server


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
            dcc.Store(id="active-run-store"),
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
        html.Style(
            """
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
            """
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
# Phase 1 placeholder callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("job-status", "children"),
    Input("run-integrated-btn", "n_clicks"),
    prevent_initial_call=True,
)
def phase1_run_placeholder(n_clicks):
    return (
        "Dashboard structure is ready. Pipeline execution wiring is the next step; "
        "this button is intentionally not launching Nextflow yet."
    )


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)

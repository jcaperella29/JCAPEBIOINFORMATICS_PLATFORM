from __future__ import annotations

from dash import Dash, Input, Output, State, dcc, html, dash_table, no_update, ctx
from pathlib import Path
import base64
import csv
import json
import html as html_std
import os
import shutil
import subprocess
import plotly.io as pio
import plotly.graph_objects as go
import networkx as nx
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



def _simple_markdown_to_html(md_text: str) -> str:
    """Lightweight Markdown renderer for the final report viewer.

    The full dcc.Markdown component can become sluggish on very large integrated
    reports. This renderer handles the structures used by our report (headings,
    bullets, paragraphs, fenced code, and simple pipe tables) and emits compact
    HTML for one section at a time.
    """
    import re as _re

    def esc(s):
        return html_std.escape(str(s), quote=False)

    lines = md_text.splitlines()
    out = []
    i = 0
    in_ul = False
    in_code = False
    code_lines = []

    def close_ul():
        nonlocal in_ul
        if in_ul:
            out.append("</ul>")
            in_ul = False

    while i < len(lines):
        line = lines[i]

        if line.strip().startswith("```"):
            close_ul()
            if not in_code:
                in_code = True
                code_lines = []
            else:
                out.append("<pre><code>" + esc("\n".join(code_lines)) + "</code></pre>")
                in_code = False
                code_lines = []
            i += 1
            continue

        if in_code:
            code_lines.append(line)
            i += 1
            continue

        # Simple pipe table: header row + separator row.
        if (
            "|" in line
            and i + 1 < len(lines)
            and _re.match(r"^\s*\|?\s*:?-{3,}", lines[i + 1])
        ):
            close_ul()
            headers = [x.strip() for x in line.strip().strip("|").split("|")]
            i += 2
            rows = []
            while i < len(lines) and "|" in lines[i] and lines[i].strip():
                rows.append([x.strip() for x in lines[i].strip().strip("|").split("|")])
                i += 1
            out.append("<div class='table-wrap'><table><thead><tr>")
            out.extend(f"<th>{esc(h)}</th>" for h in headers)
            out.append("</tr></thead><tbody>")
            for row in rows:
                out.append("<tr>")
                for cell in row:
                    out.append(f"<td>{esc(cell)}</td>")
                out.append("</tr>")
            out.append("</tbody></table></div>")
            continue

        m = _re.match(r"^(#{1,4})\s+(.*)$", line)
        if m:
            close_ul()
            level = len(m.group(1))
            out.append(f"<h{level}>{esc(m.group(2))}</h{level}>")
            i += 1
            continue

        if line.startswith("- "):
            if not in_ul:
                out.append("<ul>")
                in_ul = True
            out.append(f"<li>{esc(line[2:].strip())}</li>")
            i += 1
            continue

        if not line.strip():
            close_ul()
            i += 1
            continue

        close_ul()
        s = esc(line.strip())
        # Small amount of inline markdown support.
        s = _re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = _re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
        out.append(f"<p>{s}</p>")
        i += 1

    close_ul()
    if in_code and code_lines:
        out.append("<pre><code>" + esc("\n".join(code_lines)) + "</code></pre>")

    return "\n".join(out)


def _split_report_sections(md_text: str) -> list[dict]:
    """Split a large report at level-2 headings for light, paged rendering."""
    lines = md_text.splitlines()
    sections = []
    current_title = "Report overview"
    current_lines = []

    for line in lines:
        if line.startswith("## "):
            if current_lines:
                sections.append({
                    "title": current_title,
                    "html": _simple_markdown_to_html("\n".join(current_lines)),
                })
            current_title = line[3:].strip() or "Section"
            current_lines = [line]
        else:
            current_lines.append(line)

    if current_lines:
        sections.append({
            "title": current_title,
            "html": _simple_markdown_to_html("\n".join(current_lines)),
        })

    return [s for s in sections if s["html"].strip()]


def _render_final_report_carousel(path: Path):
    md_text = path.read_text(encoding="utf-8", errors="replace")
    sections = _split_report_sections(md_text)

    if not sections:
        return placeholder("The final report file was empty.")

    payload = json.dumps(sections).replace("</", "<\\/")

    page = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{
  margin: 0;
  font-family: Arial, sans-serif;
  color: #172033;
  background: white;
}}
.wrap {{
  padding: 14px;
}}
.toolbar {{
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 12px;
  flex-wrap: wrap;
  margin-bottom: 14px;
}}
.controls {{
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
}}
button {{
  background: white;
  border: 1px solid #d7e0eb;
  border-radius: 8px;
  padding: 8px 12px;
  font-weight: 700;
  color: #173b67;
  cursor: pointer;
}}
button:hover {{ background: #f2f6fb; }}
#counter {{
  font-weight: 700;
  color: #173b67;
}}
#section {{
  border: 1px solid #d7e0eb;
  border-radius: 12px;
  padding: 18px;
  background: #fbfdff;
  line-height: 1.55;
}}
h1, h2, h3, h4 {{
  color: #173b67;
  margin-top: 0.9em;
}}
p {{
  margin: 0.6em 0;
}}
ul {{
  margin-top: 0.4em;
}}
code {{
  background: #eef2f7;
  padding: 1px 4px;
  border-radius: 4px;
}}
pre {{
  background: #0f172a;
  color: #e5eef9;
  padding: 12px;
  border-radius: 8px;
  overflow-x: auto;
}}
.table-wrap {{
  overflow-x: auto;
  margin: 12px 0;
}}
table {{
  border-collapse: collapse;
  width: 100%;
  min-width: 600px;
}}
th, td {{
  border: 1px solid #d7e0eb;
  padding: 7px 9px;
  text-align: left;
  vertical-align: top;
}}
th {{
  background: #eef4fb;
}}
</style>
</head>
<body>
<div class="wrap">
  <div class="toolbar">
    <div class="controls">
      <button id="prev">← Previous section</button>
      <button id="next">Next section →</button>
    </div>
    <div id="counter"></div>
  </div>
  <div id="section"></div>
</div>
<script>
const sections = {payload};
let idx = 0;

function render() {{
  const s = sections[idx];
  document.getElementById("counter").textContent =
    (idx + 1) + " of " + sections.length + " — " + s.title;
  document.getElementById("section").innerHTML = s.html;
}}

document.getElementById("prev").addEventListener("click", () => {{
  idx = (idx - 1 + sections.length) % sections.length;
  render();
}});
document.getElementById("next").addEventListener("click", () => {{
  idx = (idx + 1) % sections.length;
  render();
}});
render();
</script>
</body>
</html>"""

    return html.Iframe(
        srcDoc=page,
        style={
            "width": "100%",
            "height": "900px",
            "border": f"1px solid {COLORS['border']}",
            "borderRadius": "10px",
            "background": "white",
        },
    )


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
    priority = [
        "pca_plot.html",
        "umap_plot.html",
        "volcano_plot.html",
        "heatmap_plot.html",
    ]
    chosen = []
    seen = set()

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


def _build_plot_carousel_html(plots: list[Path]) -> str:
    """Self-contained HTML/JS plot carousel.

    This avoids Dash callbacks for Previous/Next buttons, which were causing
    missing-component errors because the plot controls lived in dynamic tabs.
    """
    payload = []
    for p in plots:
        payload.append({
            "name": p.name,
            "html": p.read_text(encoding="utf-8", errors="replace"),
        })

    data_json = json.dumps(payload).replace("</", "<\\/")

    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{
  font-family: Arial, sans-serif;
  margin: 0;
  background: white;
  color: #172033;
}}
.toolbar {{
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 10px;
  padding: 10px 12px;
  border-bottom: 1px solid #d7e0eb;
  flex-wrap: wrap;
}}
.controls {{
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
}}
button {{
  background: white;
  border: 1px solid #d7e0eb;
  border-radius: 8px;
  padding: 8px 11px;
  font-weight: 650;
  color: #173b67;
  cursor: pointer;
}}
button:hover {{ background: #f2f6fb; }}
#label {{ font-weight: 700; color: #173b67; }}
#plotFrame {{
  width: 100%;
  height: 720px;
  border: 0;
}}
</style>
</head>
<body>
<div class="toolbar">
  <div class="controls">
    <button onclick="prevPlot()">← Previous</button>
    <button onclick="nextPlot()">Next →</button>
  </div>
  <div id="label"></div>
  <button onclick="downloadCurrent()">Download current plot</button>
</div>
<iframe id="plotFrame"></iframe>
<script>
const plots = {data_json};
let idx = 0;

function renderPlot() {{
  if (!plots.length) return;
  const p = plots[idx];
  document.getElementById("label").textContent =
    (idx + 1) + " of " + plots.length + " — " + p.name;
  document.getElementById("plotFrame").srcdoc = p.html;
}}

function prevPlot() {{
  idx = (idx - 1 + plots.length) % plots.length;
  renderPlot();
}}

function nextPlot() {{
  idx = (idx + 1) % plots.length;
  renderPlot();
}}

function downloadCurrent() {{
  const p = plots[idx];
  const blob = new Blob([p.html], {{type: "text/html;charset=utf-8"}});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = p.name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}}

renderPlot();
</script>
</body>
</html>"""


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
                    ("Graphs", "graphs"),
                    ("Nodes & Edges", "nodes-edges"),
                    ("Bipartite Diffusion", "bipartite"),
                    ("Projection Diffusion", "projection-diffusion"),
                    ("Consensus", "consensus"),
                    ("Tables / Files", "files"),
                ],
                "graphs",
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
        return card(
            [
                html.H3("Final Report", style={"color": COLORS["navy"], "marginTop": 0}),
                html.Div(
                    "Download the completed integrated report in Markdown or JSON format.",
                    style={"color": COLORS["muted"], "marginBottom": "14px"},
                ),
                html.Div(
                    [
                        html.Button("Download Markdown", id="rna-final-download-markdown", n_clicks=0, className="export-btn"),
                        html.Button("Download JSON", id="rna-final-download-json", n_clicks=0, className="export-btn"),
                        dcc.Download(id="rna-final-report-download"),
                    ],
                    style={"display": "flex", "gap": "10px", "flexWrap": "wrap"},
                ),
            ]
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
# RNA result-population callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("rna-omics-inner-panel", "children"),
    Input("rna-omics-inner-tabs", "value"),
    Input("active-run-store", "data"),
)
def populate_rna_omics(tab, store):
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
        return html.Iframe(
            srcDoc=_build_plot_carousel_html(plots),
            style={
                "width": "100%",
                "height": "790px",
                "border": f"1px solid {COLORS['border']}",
                "borderRadius": "10px",
                "background": "white",
            },
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



def _figure_from_graph_store(store: dict, title: str):
    """Build a standalone Plotly network from the API's saved graph store.

    /network/build intentionally saves graph data but no figure. This reconstructs
    the main bipartite graph from that store so the dashboard can display it next
    to the projection figure without another API call.
    """
    nodes = store.get("nodes") or []
    edges = store.get("edges") or []
    if not nodes or not edges:
        return None

    g = nx.Graph()

    for node in nodes:
        node_id = (
            node.get("id")
            or node.get("node_id")
            or node.get("name")
            or node.get("label")
        )
        if node_id is None:
            continue
        g.add_node(
            str(node_id),
            label=str(node.get("label") or node.get("name") or node_id),
            node_type=str(node.get("node_type") or node.get("type") or "unknown"),
        )

    for edge in edges:
        source = edge.get("source") or edge.get("from")
        target = edge.get("target") or edge.get("to")
        if source is None or target is None:
            continue
        source = str(source)
        target = str(target)
        if source not in g:
            g.add_node(source, label=source, node_type="unknown")
        if target not in g:
            g.add_node(target, label=target, node_type="unknown")
        weight = (
            edge.get("weight")
            or edge.get("weight_plot")
            or edge.get("adjusted_pvalue")
            or edge.get("value")
            or 1.0
        )
        try:
            weight = float(weight)
        except Exception:
            weight = 1.0
        g.add_edge(source, target, weight=weight)

    if g.number_of_nodes() == 0:
        return None

    # Prefer a true two-column bipartite arrangement: genes/items on the left,
    # pathways/groups on the right. Unknown nodes fall back to spring layout.
    left = [
        n for n, d in g.nodes(data=True)
        if str(d.get("node_type", "")).lower() in {"item", "gene"}
    ]
    right = [
        n for n, d in g.nodes(data=True)
        if str(d.get("node_type", "")).lower() in {"group", "term", "pathway"}
    ]

    if left and right and (len(left) + len(right) >= int(0.8 * g.number_of_nodes())):
        pos = {}
        left_sorted = sorted(left, key=lambda n: g.nodes[n].get("label", n))
        right_sorted = sorted(right, key=lambda n: g.nodes[n].get("label", n))

        def ys(count):
            if count <= 1:
                return [0.0]
            return [1.0 - 2.0 * i / (count - 1) for i in range(count)]

        for node, y in zip(left_sorted, ys(len(left_sorted))):
            pos[node] = (-1.0, y)
        for node, y in zip(right_sorted, ys(len(right_sorted))):
            pos[node] = (1.0, y)

        missing = [n for n in g.nodes if n not in pos]
        if missing:
            extra = nx.spring_layout(g.subgraph(missing), seed=7)
            for n, xy in extra.items():
                pos[n] = (float(xy[0]) * 0.4, float(xy[1]))
    else:
        pos = nx.spring_layout(g, seed=7)

    edge_x, edge_y = [], []
    for u, v in g.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        hoverinfo="none",
        line={"width": 0.7},
        opacity=0.45,
        showlegend=False,
    )

    node_x, node_y, labels, hover = [], [], [], []
    for n, data in g.nodes(data=True):
        x, y = pos[n]
        node_x.append(x)
        node_y.append(y)
        label = str(data.get("label") or n)
        labels.append(label)
        hover.append(
            f"{label}<br>Type: {data.get('node_type', 'unknown')}<br>Degree: {g.degree(n)}"
        )

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        text=hover,
        hoverinfo="text",
        marker={"size": 9},
        showlegend=False,
    )

    fig = go.Figure([edge_trace, node_trace])
    fig.update_layout(
        title=title,
        hovermode="closest",
        margin={"l": 20, "r": 20, "t": 55, "b": 20},
        xaxis={"visible": False},
        yaxis={"visible": False},
        height=720,
    )
    return fig


def _network_figure_html(path: Path, title: str) -> str | None:
    """Return a standalone HTML network figure from saved API output."""
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))

        # Projection/filter endpoints already save a Plotly figure.
        fig = obj.get("figure")
        if isinstance(fig, dict):
            return pio.to_html(
                fig,
                include_plotlyjs=True,
                full_html=True,
                config={"displaylogo": False, "responsive": True},
            )

        # /network/build saves the graph store but intentionally no figure.
        store = obj.get("graph")
        if isinstance(store, dict):
            rebuilt = _figure_from_graph_store(store, title)
            if rebuilt is not None:
                return pio.to_html(
                    rebuilt,
                    include_plotlyjs=True,
                    full_html=True,
                    config={"displaylogo": False, "responsive": True},
                )

        return None
    except Exception:
        return None


def _build_network_carousel_html(items: list[dict]) -> str:
    """Browser-side carousel for standalone network Plotly figures.

    Each item is {"label": ..., "html": ...}. This keeps graph switching entirely
    inside the iframe and avoids Dash's lazy-loaded dcc.Graph JavaScript chunk.
    """
    payload = json.dumps(items).replace("</", "<\\/")

    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{
  font-family: Arial, sans-serif;
  margin: 0;
  background: white;
  color: #172033;
}}
.toolbar {{
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 10px;
  padding: 10px 12px;
  border-bottom: 1px solid #d7e0eb;
  flex-wrap: wrap;
}}
.controls {{
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
}}
button {{
  background: white;
  border: 1px solid #d7e0eb;
  border-radius: 8px;
  padding: 8px 11px;
  font-weight: 650;
  color: #173b67;
  cursor: pointer;
}}
button:hover {{ background: #f2f6fb; }}
#label {{ font-weight: 700; color: #173b67; }}
#graphFrame {{
  width: 100%;
  height: 760px;
  border: 0;
}}
</style>
</head>
<body>
<div class="toolbar">
  <div class="controls">
    <button onclick="prevGraph()">← Previous graph</button>
    <button onclick="nextGraph()">Next graph →</button>
  </div>
  <div id="label"></div>
</div>
<iframe id="graphFrame"></iframe>
<script>
const items = {payload};
let idx = 0;

function currentFileName() {{
  if (!items.length) return "network_graph.html";
  const raw = items[idx].label || ("network_graph_" + (idx + 1));
  return raw.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_+|_+$/g, "") + ".html";
}}

function downloadCurrentGraph() {{
  if (!items.length) return;
  const item = items[idx];
  const blob = new Blob([item.html], {{type: "text/html;charset=utf-8"}});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = currentFileName();
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}}

function wireParentDownloadButton() {{
  try {{
    const btn = parent.document.getElementById("rna-network-download-plot");
    if (!btn || btn.dataset.jcapGraphDownloadWired === "1") return;

    btn.dataset.jcapGraphDownloadWired = "1";
    btn.addEventListener("click", function(ev) {{
      ev.preventDefault();
      ev.stopImmediatePropagation();
      downloadCurrentGraph();
    }}, true);
  }} catch (e) {{
    // Graph viewer still works if parent access is unavailable.
  }}
}}

function renderGraph() {{
  if (!items.length) return;
  const item = items[idx];
  document.getElementById("label").textContent =
    (idx + 1) + " of " + items.length + " — " + item.label;
  document.getElementById("graphFrame").srcdoc = item.html;
  wireParentDownloadButton();
}}

function prevGraph() {{
  idx = (idx - 1 + items.length) % items.length;
  renderGraph();
}}

function nextGraph() {{
  idx = (idx + 1) % items.length;
  renderGraph();
}}

renderGraph();
setTimeout(wireParentDownloadButton, 250);
</script>
</body>
</html>"""


def _render_network_graph_carousel(d: Path):
    candidates = [
        ("Main enrichment network", _first_existing(d, ["network_build.json"])),
        ("Pathway projection / diffusion network", _first_existing(d, ["projection_build.json"])),
    ]

    items = []
    for label, path in candidates:
        if not path:
            continue
        graph_html = _network_figure_html(path, label)
        if graph_html:
            items.append({"label": label, "html": graph_html})

    if not items:
        return placeholder(
            "No renderable network figure was found in network_build.json or projection_build.json."
        )

    return html.Iframe(
        srcDoc=_build_network_carousel_html(items),
        style={
            "width": "100%",
            "height": "825px",
            "border": f"1px solid {COLORS['border']}",
            "borderRadius": "10px",
            "background": "white",
        },
    )


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
        "nodes-edges": (
            ["main_nodes.csv", "main_edges.csv", "projection_nodes.csv", "projection_edges.csv"],
            None,
        ),
        "bipartite": (
            ["bipartite_diffusion_results.csv", "bipartite_top_candidates.csv"],
            "bipartite_diffusion.json",
        ),
        "projection-diffusion": (
            ["projection_diffusion_results.csv", "projection_top_candidates.csv"],
            "projection_diffusion.json",
        ),
        "consensus": (
            ["consensus_candidates.csv"],
            "consensus.json",
        ),
    }

    if tab == "graphs":
        return _render_network_graph_carousel(d)

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



def _load_analysis_result(d: Path) -> dict:
    p = _first_existing(d, ["analysis_result.json"])
    if not p:
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _render_triage_overview(obj: dict):
    metadata = obj.get("metadata") or {}
    validation = obj.get("input_validation") or {}
    pubmed = obj.get("pubmed") or {}
    playbook = obj.get("playbook_context") or {}
    claims = obj.get("claims") or []
    programs = (obj.get("programs") or {}).get("programs") or []
    triage_rows = (obj.get("triage") or {}).get("rows") or []

    cards = [
        ("Triage rows", len(triage_rows)),
        ("Programs", len(programs)),
        ("Claims", len(claims)),
        ("PubMed papers", len(pubmed.get("papers") or [])),
    ]

    return html.Div([
        html.Div(
            [
                html.Div(
                    [
                        html.Div(str(value), style={"fontSize": "28px", "fontWeight": "800", "color": COLORS["navy"]}),
                        html.Div(label, style={"color": COLORS["muted"]}),
                    ],
                    style={
                        "background": "#fbfdff",
                        "border": f"1px solid {COLORS['border']}",
                        "borderRadius": "10px",
                        "padding": "14px",
                    },
                )
                for label, value in cards
            ],
            style={"display": "grid", "gridTemplateColumns": "repeat(4, minmax(140px, 1fr))", "gap": "10px", "marginBottom": "16px"},
        ),
        html.H4("Run metadata", style={"color": COLORS["navy"]}),
        html.Pre(
            json.dumps(metadata, indent=2, default=str),
            style={"whiteSpace": "pre-wrap", "background": "#f8fafc", "padding": "10px", "borderRadius": "8px"},
        ),
        html.H4("Input validation", style={"color": COLORS["navy"]}),
        html.Pre(
            json.dumps(validation, indent=2, default=str),
            style={"whiteSpace": "pre-wrap", "background": "#f8fafc", "padding": "10px", "borderRadius": "8px"},
        ),
        html.Div(
            [
                html.Strong("Literature retrieval: "),
                html.Span(
                    f"{pubmed.get('retrieval_quality', 'unknown')} — {pubmed.get('retrieval_quality_reason', '')}"
                ),
            ],
            style={"marginTop": "10px"},
        ),
        html.Div(
            [
                html.Strong("Playbook retrieval mode: "),
                html.Span(playbook.get("mode", "unknown")),
            ],
            style={"marginTop": "6px"},
        ),
    ])


def _render_pubmed_context(obj: dict):
    pubmed = obj.get("pubmed") or {}
    papers = pubmed.get("papers") or []

    summary = html.Div(
        [
            html.Div([html.Strong("Status: "), html.Span(str(pubmed.get("status", "")))]),
            html.Div([html.Strong("Source: "), html.Span(str(pubmed.get("source", "")))]),
            html.Div([html.Strong("Query used: "), html.Code(str(pubmed.get("query_used") or pubmed.get("query") or ""))]),
            html.Div([html.Strong("Query strategy: "), html.Span(str(pubmed.get("query_strategy", "")))]),
            html.Div([html.Strong("Retrieval quality: "), html.Span(str(pubmed.get("retrieval_quality", "")))]),
            html.Div(str(pubmed.get("retrieval_quality_reason", "")), style={"color": COLORS["muted"], "marginTop": "4px"}),
            html.Div(
                [
                    html.Strong("Relevance counts: "),
                    html.Span(
                        ", ".join(
                            f"{k}: {v}" for k, v in (pubmed.get("relevance_counts") or {}).items()
                        ) or "none"
                    ),
                ],
                style={"marginTop": "6px"},
            ),
            html.Div(
                str(pubmed.get("literature_use_guidance", "")),
                style={
                    "marginTop": "10px",
                    "padding": "10px",
                    "background": COLORS["warn"],
                    "borderRadius": "8px",
                    "color": COLORS["navy"],
                },
            ),
        ],
        style={
            "padding": "14px",
            "background": "#fbfdff",
            "border": f"1px solid {COLORS['border']}",
            "borderRadius": "10px",
            "marginBottom": "16px",
        },
    )

    if not papers:
        return html.Div([summary, placeholder("No PubMed papers were retrieved for this run.")])

    rows = []
    for p in papers:
        rows.append({
            "Title": p.get("title", ""),
            "PMID": p.get("pmid", ""),
            "Date": p.get("pubdate", ""),
            "Source": p.get("source", ""),
            "DOI": p.get("doi", ""),
            "Relevance": p.get("relevance_label", ""),
            "Score": p.get("relevance_score", ""),
            "Why": p.get("relevance_reason", ""),
            "Use": p.get("literature_use", ""),
        })

    table = dash_table.DataTable(
        data=rows,
        columns=[{"name": c, "id": c} for c in rows[0].keys()],
        page_size=10,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={
            "textAlign": "left",
            "fontSize": "12px",
            "padding": "7px",
            "whiteSpace": "normal",
            "height": "auto",
            "minWidth": "90px",
            "maxWidth": "420px",
        },
        style_header={"fontWeight": "700", "background": "#eef4fb"},
    )

    abstract_blocks = []
    for i, p in enumerate(papers, start=1):
        abstract = (p.get("abstract") or "").strip()
        if not abstract:
            continue
        abstract_blocks.append(
            html.Details(
                [
                    html.Summary(f"{i}. {p.get('title') or 'Untitled'}"),
                    html.P(abstract, style={"lineHeight": "1.5"}),
                ],
                style={"marginTop": "8px"},
            )
        )

    return html.Div([
        summary,
        html.H4("Retrieved papers", style={"color": COLORS["navy"]}),
        table,
        html.H4("Abstracts", style={"color": COLORS["navy"], "marginTop": "18px"}),
        html.Div(abstract_blocks) if abstract_blocks else html.Div("No abstracts were available.", style={"color": COLORS["muted"]}),
    ])



def _render_claim_validation_carousel(obj: dict):
    claims = obj.get("claims") or []
    if not claims:
        return placeholder("No structured claims were available for validation review.")

    slides = []
    for i, claim in enumerate(claims, start=1):
        program = claim.get("program") or f"Claim {i}"
        role = claim.get("role") or ""
        strength = claim.get("evidence_strength") or ""
        rationale = claim.get("rationale") or claim.get("claim") or claim.get("summary") or ""
        confounders = claim.get("confounders") or []
        validations = claim.get("validation") or []
        evidence = claim.get("evidence") or {}

        ev_terms = evidence.get("terms") or []
        ev_genes = evidence.get("genes") or []
        ev_pmids = evidence.get("pmids") or []

        validation_html = []
        for j, v in enumerate(validations, start=1):
            validation_html.append(
                f"""
                <div class="val-card">
                  <div class="val-head">Validation {j}</div>
                  <div><b>Experiment:</b> {html_std.escape(str(v.get('experiment','')))}</div>
                  <div><b>Readout:</b> {html_std.escape(str(v.get('readout','')))}</div>
                  <div><b>Control:</b> {html_std.escape(str(v.get('control','')))}</div>
                  <div><b>Expected if causal:</b> {html_std.escape(str(v.get('expected_result_if_causal','')))}</div>
                  <div><b>Purpose:</b> {html_std.escape(str(v.get('purpose','')))}</div>
                </div>
                """
            )

        slides.append({
            "title": program,
            "role": role,
            "strength": strength,
            "rationale": rationale,
            "confounders": confounders,
            "terms": ev_terms,
            "genes": ev_genes,
            "pmids": ev_pmids,
            "validation_html": "".join(validation_html) or '<div class="muted">No validation experiment listed.</div>',
        })

    payload = json.dumps(slides).replace("</", "<\\/")

    page = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{
  margin: 0;
  font-family: Arial, sans-serif;
  color: #172033;
  background: white;
}}
.wrap {{
  padding: 14px;
}}
.toolbar {{
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 12px;
  flex-wrap: wrap;
  margin-bottom: 14px;
}}
.controls {{
  display: flex;
  gap: 8px;
}}
button {{
  background: white;
  border: 1px solid #d7e0eb;
  border-radius: 8px;
  padding: 8px 12px;
  font-weight: 700;
  color: #173b67;
  cursor: pointer;
}}
button:hover {{ background: #f2f6fb; }}
.counter {{
  font-weight: 700;
  color: #173b67;
}}
.hero {{
  border: 1px solid #d7e0eb;
  border-radius: 12px;
  padding: 16px;
  background: #fbfdff;
  margin-bottom: 14px;
}}
.title {{
  font-size: 22px;
  font-weight: 800;
  color: #173b67;
  margin-bottom: 8px;
}}
.badges {{
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
  margin-bottom: 12px;
}}
.badge {{
  border: 1px solid #d7e0eb;
  border-radius: 999px;
  padding: 5px 9px;
  background: white;
  font-size: 12px;
  font-weight: 700;
}}
.grid {{
  display: grid;
  grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
  gap: 14px;
}}
.panel {{
  border: 1px solid #d7e0eb;
  border-radius: 10px;
  padding: 14px;
  min-width: 0;
}}
.panel h3 {{
  margin: 0 0 10px 0;
  color: #173b67;
}}
.muted {{
  color: #66758a;
}}
ul {{
  margin: 6px 0 0 18px;
  padding: 0;
}}
.val-card {{
  border-top: 1px solid #e5eaf1;
  padding-top: 10px;
  margin-top: 10px;
  line-height: 1.45;
}}
.val-card:first-child {{
  border-top: 0;
  padding-top: 0;
  margin-top: 0;
}}
.val-head {{
  font-weight: 800;
  color: #173b67;
  margin-bottom: 6px;
}}
@media (max-width: 720px) {{
  .grid {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<div class="wrap">
  <div class="toolbar">
    <div class="controls">
      <button id="prev">← Previous claim</button>
      <button id="next">Next claim →</button>
    </div>
    <div class="counter" id="counter"></div>
  </div>

  <div class="hero">
    <div class="title" id="title"></div>
    <div class="badges">
      <span class="badge" id="role"></span>
      <span class="badge" id="strength"></span>
    </div>
    <div id="rationale"></div>
  </div>

  <div class="grid">
    <div class="panel">
      <h3>Evidence</h3>
      <div><b>Genes</b><ul id="genes"></ul></div>
      <div style="margin-top:10px"><b>Terms</b><ul id="terms"></ul></div>
      <div style="margin-top:10px"><b>PMIDs</b><ul id="pmids"></ul></div>
    </div>
    <div class="panel">
      <h3>Confounders</h3>
      <ul id="confounders"></ul>
    </div>
  </div>

  <div class="panel" style="margin-top:14px">
    <h3>Validation plan</h3>
    <div id="validation"></div>
  </div>
</div>

<script>
const slides = {payload};
let idx = 0;

function fillList(id, arr) {{
  const el = document.getElementById(id);
  el.innerHTML = "";
  if (!arr || !arr.length) {{
    const li = document.createElement("li");
    li.textContent = "None listed";
    li.className = "muted";
    el.appendChild(li);
    return;
  }}
  arr.forEach(v => {{
    const li = document.createElement("li");
    li.textContent = String(v);
    el.appendChild(li);
  }});
}}

function render() {{
  const s = slides[idx];
  document.getElementById("counter").textContent = `${{idx + 1}} of ${{slides.length}}`;
  document.getElementById("title").textContent = s.title;
  document.getElementById("role").textContent = s.role || "Role not specified";
  document.getElementById("strength").textContent = s.strength || "Evidence strength not specified";
  document.getElementById("rationale").textContent = s.rationale || "";
  fillList("genes", s.genes);
  fillList("terms", s.terms);
  fillList("pmids", s.pmids);
  fillList("confounders", s.confounders);
  document.getElementById("validation").innerHTML = s.validation_html;
}}

document.getElementById("prev").addEventListener("click", () => {{
  idx = (idx - 1 + slides.length) % slides.length;
  render();
}});
document.getElementById("next").addEventListener("click", () => {{
  idx = (idx + 1) % slides.length;
  render();
}});
render();
</script>
</body>
</html>"""

    return html.Iframe(
        srcDoc=page,
        style={
            "width": "100%",
            "height": "850px",
            "border": f"1px solid {COLORS['border']}",
            "borderRadius": "10px",
            "background": "white",
        },
    )


def _render_playbook_guidance(obj: dict):
    playbook = obj.get("playbook_context") or {}
    context_text = (playbook.get("context_text") or "").strip()

    pieces = [
        html.Div(
            [
                html.Div([html.Strong("Retrieval mode: "), html.Span(str(playbook.get("mode", "")))]),
                html.Div([html.Strong("Query: "), html.Span(str(playbook.get("query", "")))]),
            ],
            style={
                "padding": "14px",
                "background": "#fbfdff",
                "border": f"1px solid {COLORS['border']}",
                "borderRadius": "10px",
                "marginBottom": "12px",
            },
        )
    ]

    if playbook.get("error"):
        pieces.append(
            html.Div(
                str(playbook.get("error")),
                style={
                    "padding": "10px",
                    "background": COLORS["warn"],
                    "borderRadius": "8px",
                    "marginBottom": "12px",
                    "color": COLORS["navy"],
                },
            )
        )

    pieces.append(
        dcc.Markdown(context_text or "No playbook guidance was retrieved.")
    )
    return html.Div(pieces)


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
        obj = _load_analysis_result(d)
        return _render_triage_overview(obj) if obj else placeholder("analysis_result.json was not found or could not be read.")
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
        obj = _load_analysis_result(d)
        return _render_pubmed_context(obj) if obj else placeholder("analysis_result.json was not found or could not be read.")
    if tab == "playbook":
        obj = _load_analysis_result(d)
        return _render_playbook_guidance(obj) if obj else placeholder("analysis_result.json was not found or could not be read.")
    if tab == "validation":
        obj = _load_analysis_result(d)
        return _render_claim_validation_carousel(obj) if obj else placeholder("analysis_result.json was not found or could not be read.")
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


def _rna_current_file(prefix: str, tab: str, root: Path, kind: str) -> Path | None:
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
            mapping = {
                "plots": ["pca_plot.html", "umap_plot.html", "volcano_plot.html", "heatmap_plot.html", "*.html"],
                "enrichment": ["*enrichment*.html"],
                "classifier": ["rf_roc_plot.html"],
                "power": ["*power*.html"],
            }
    elif prefix == "rna-network":
        d = _section_dir(root, "02_network")
        mapping = {
            "nodes-edges": ["main_nodes.csv", "main_edges.csv", "projection_nodes.csv", "projection_edges.csv"],
            "bipartite": ["bipartite_diffusion_results.csv", "bipartite_top_candidates.csv"],
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


def _rna_download(prefix: str, tab: str, store: dict | None):
    triggered = ctx.triggered_id
    if not triggered:
        return no_update

    root = _run_result_root(store)
    if not root:
        return no_update

    if triggered.endswith("-download-table"):
        path = _rna_current_file(prefix, tab, root, "table")
        return dcc.send_file(str(path)) if path and path.exists() else no_update

    if triggered.endswith("-download-plot"):
        # The RNA Omics Plots view has an in-view browser-side download button
        # that follows the currently displayed carousel plot.
        if prefix == "rna-omics" and tab == "plots":
            return no_update
        path = _rna_current_file(prefix, tab, root, "plot")
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
        prevent_initial_call=True,
    )
    def _download(_n1, _n2, _n3, _n4, tab, store):
        return _rna_download(prefix, tab, store)


for _prefix in ["rna-omics", "rna-network", "rna-triage", "rna-literature"]:
    _register_rna_download_callback(_prefix)


# -----------------------------------------------------------------------------
# RNA final-report downloads
# -----------------------------------------------------------------------------
@app.callback(
    Output("rna-final-report-download", "data"),
    Input("rna-final-download-markdown", "n_clicks"),
    Input("rna-final-download-json", "n_clicks"),
    State("active-run-store", "data"),
    prevent_initial_call=True,
)
def download_rna_final_report(_md_clicks, _json_clicks, store):
    root = _run_result_root(store)
    if not root:
        return no_update

    triggered = ctx.triggered_id
    if triggered == "rna-final-download-markdown":
        path = _first_existing(root, ["rnaseq_final_report.md"])
    elif triggered == "rna-final-download-json":
        path = _first_existing(root, ["rnaseq_final_report.json"])
    else:
        return no_update

    return dcc.send_file(str(path)) if path and path.exists() else no_update


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

from __future__ import annotations

import io
import json
import os
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
from dotenv import load_dotenv

# Load .env BEFORE importing repo modules.
# reasoner.py constructs the OpenAI client at import time, so OPENAI_API_KEY must exist first.
load_dotenv(override=True)

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.encoders import jsonable_encoder
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.units import inch
from xml.sax.saxutils import escape

# Existing repo modules. Keep the rest of the repository as-is.
# Do NOT import pipeline.py here; pipeline imports reasoner/OpenAI.
# We lazy-import it only when mode=full is requested.
from input_validation import validate_enrichment_df
from program_summarizer import summarize_programs
from summarizer import build_triage_pdf
from triage import triage_enrichment_table

APP_VERSION = os.environ.get("APP_VERSION", "0.3.3-fastapi-readable-reports")
REPORTS_DIR = Path(os.environ.get("REPORTS_DIR") or Path(tempfile.gettempdir()) / "enrichment_llm_reports")
REPORTS_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(
    title="Enrichment Analysis LLM API",
    version=APP_VERSION,
    description=(
        "FastAPI module wrapper for the Enrichment Analysis LLM app. "
        "Keeps existing pipeline/triage/reasoner/summarizer files unchanged."
    ),
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=os.environ.get("CORS_ALLOW_ORIGINS", "*").split(","),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------------------
# Small compatibility helpers
# ---------------------------------------------------------------------


def _clean_for_json(value: Any) -> Any:
    """Make numpy/pandas values safe for JSON responses."""
    if isinstance(value, dict):
        return {str(k): _clean_for_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_clean_for_json(v) for v in value]
    if isinstance(value, tuple):
        return [_clean_for_json(v) for v in value]
    if isinstance(value, set):
        return sorted(_clean_for_json(v) for v in value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        v = float(value)
        return None if np.isnan(v) or np.isinf(v) else v
    if isinstance(value, float):
        return None if np.isnan(value) or np.isinf(value) else value
    if pd.isna(value) if not isinstance(value, (dict, list, tuple, set)) else False:
        return None
    return value


def _json_response(payload: Any, status_code: int = 200) -> JSONResponse:
    return JSONResponse(content=jsonable_encoder(_clean_for_json(payload)), status_code=status_code)


def _read_upload_csv(file: UploadFile) -> pd.DataFrame:
    try:
        raw = file.file.read()
        return pd.read_csv(io.BytesIO(raw))
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Could not read uploaded CSV: {exc}") from exc


def _context_from_form(
    *,
    organism: str = "",
    assay: str = "",
    tissue: str = "",
    cell_type: str = "",
    perturbation: str = "",
    timepoint: str = "",
    extra_context_json: str | None = None,
) -> Dict[str, Any]:
    context: Dict[str, Any] = {
        "organism": organism or "",
        "assay": assay or "",
        "tissue": tissue or "",
        "cell_type": cell_type or "",
        "perturbation": perturbation or "",
        "timepoint": timepoint or "",
    }
    if extra_context_json:
        try:
            extra = json.loads(extra_context_json)
        except Exception as exc:
            raise HTTPException(status_code=400, detail=f"extra_context_json is not valid JSON: {exc}") from exc
        if not isinstance(extra, dict):
            raise HTTPException(status_code=400, detail="extra_context_json must be a JSON object.")
        context.update(extra)
    return context


def _normalize_col_name(name: str) -> str:
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    normalized = {_normalize_col_name(c): c for c in df.columns}
    for candidate in candidates:
        hit = normalized.get(_normalize_col_name(candidate))
        if hit is not None:
            return hit
    return None


def _long_edges_to_enrichment_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert RNA-seq/Network API long edge files into the LLM app's expected shape.

    Input:
      gene,term,adjusted_pvalue

    Output:
      Term,Genes,Adjusted.P.value,Overlap
    """
    gene_col = _find_col(df, ["gene", "genes", "item", "symbol", "gene_symbol"])
    term_col = _find_col(df, ["term", "pathway", "description", "group", "name"])
    padj_col = _find_col(df, ["adjusted_pvalue", "adjusted p-value", "Adjusted.P.value", "padj", "fdr", "qvalue", "q_value"])

    if not gene_col or not term_col or not padj_col:
        raise ValueError(
            "Long-edge conversion requires gene, term, and adjusted p-value columns. "
            f"Found columns: {list(df.columns)}"
        )

    tmp = df[[gene_col, term_col, padj_col]].copy()
    tmp[gene_col] = tmp[gene_col].astype(str).str.strip()
    tmp[term_col] = tmp[term_col].astype(str).str.strip()
    tmp[padj_col] = pd.to_numeric(tmp[padj_col], errors="coerce")
    tmp = tmp[(tmp[gene_col] != "") & (tmp[term_col] != "")]

    rows = []
    for term, grp in tmp.groupby(term_col, sort=False):
        genes = sorted(set(g for g in grp[gene_col].astype(str) if g and g.lower() not in {"nan", "none", "null"}))
        if not genes:
            continue
        padj = grp[padj_col].dropna().min()
        rows.append(
            {
                "Term": term,
                "Genes": ";".join(genes),
                "Adjusted.P.value": None if pd.isna(padj) else float(padj),
                "Overlap": f"{len(genes)}/{len(genes)}",
            }
        )
    return pd.DataFrame(rows)


def _normalize_input_df(df: pd.DataFrame, input_format: str = "auto") -> pd.DataFrame:
    fmt = (input_format or "auto").strip().lower()
    if fmt not in {"auto", "enrichment", "long_edges"}:
        raise HTTPException(status_code=400, detail="input_format must be auto, enrichment, or long_edges.")

    if fmt == "enrichment":
        return df

    has_gene = _find_col(df, ["gene", "item", "symbol", "gene_symbol"]) is not None
    has_term = _find_col(df, ["term", "pathway", "description", "group", "name"]) is not None
    has_padj = _find_col(df, ["adjusted_pvalue", "Adjusted.P.value", "padj", "fdr", "qvalue", "q_value"]) is not None
    has_gene_list = _find_col(df, ["Genes", "genes", "overlap_genes", "gene_list", "core_enrichment"]) is not None

    if fmt == "long_edges" or (has_gene and has_term and has_padj and not has_gene_list):
        try:
            return _long_edges_to_enrichment_table(df)
        except Exception as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    return df


def _write_result_tables(result: dict, outdir: Path) -> dict[str, str]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    result_json = outdir / "analysis_result.json"
    result_json.write_text(json.dumps(_clean_for_json(result), indent=2, default=str), encoding="utf-8")
    paths["analysis_result_json"] = str(result_json)

    triage_rows = result.get("triage", {}).get("rows", []) or []
    if triage_rows:
        p = outdir / "triage_rows.csv"
        pd.DataFrame(triage_rows).to_csv(p, index=False)
        paths["triage_rows_csv"] = str(p)

    programs = result.get("programs", {}).get("programs", []) or []
    if programs:
        p = outdir / "programs.csv"
        pd.DataFrame(programs).to_csv(p, index=False)
        paths["programs_csv"] = str(p)

    claims = result.get("claims", []) or []
    if claims:
        p = outdir / "claims.csv"
        pd.json_normalize(claims).to_csv(p, index=False)
        paths["claims_csv"] = str(p)

    markdown = result.get("markdown_report") or _make_markdown_report(result)
    if markdown:
        p = outdir / "report.md"
        p.write_text(markdown, encoding="utf-8")
        paths["markdown_report"] = str(p)

    return paths


def _markdown_to_pdf(markdown_text: str, out_pdf_path: Path, *, title: str = "Evidence-Aware Enrichment Interpretation") -> None:
    """Small, dependency-light Markdown-ish PDF renderer for API bundles."""
    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(
        str(out_pdf_path),
        pagesize=letter,
        leftMargin=0.65 * inch,
        rightMargin=0.65 * inch,
        topMargin=0.65 * inch,
        bottomMargin=0.65 * inch,
        title=title,
    )
    story = []
    story.append(Paragraph(escape(title), styles["Title"]))
    story.append(Paragraph(escape(f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"), styles["BodyText"]))
    story.append(Spacer(1, 0.18 * inch))

    def clean_inline(line: str) -> str:
        # Basic markdown cleanup before ReportLab paragraph rendering.
        line = line.replace("**", "")
        line = line.replace("__", "")
        return escape(line)

    for raw_line in (markdown_text or "").splitlines():
        line = raw_line.rstrip()
        if not line.strip():
            story.append(Spacer(1, 0.08 * inch))
            continue
        if line.startswith("### "):
            story.append(Paragraph(clean_inline(line[4:]), styles["Heading3"]))
        elif line.startswith("## "):
            story.append(Paragraph(clean_inline(line[3:]), styles["Heading2"]))
        elif line.startswith("# "):
            story.append(Paragraph(clean_inline(line[2:]), styles["Heading1"]))
        elif line.startswith("- ") or line.startswith("* "):
            story.append(Paragraph("• " + clean_inline(line[2:]), styles["BodyText"]))
        else:
            story.append(Paragraph(clean_inline(line), styles["BodyText"]))
    doc.build(story)


def _build_pdf(result: dict, out_pdf_path: Path) -> None:
    context = result.get("context") or {}
    assay = context.get("assay", "") if isinstance(context, dict) else ""
    title = f"{assay} Evidence-Aware Enrichment Interpretation" if assay else "Evidence-Aware Enrichment Interpretation"
    markdown = result.get("markdown_report") or _make_markdown_report(result)
    if not markdown or not str(markdown).strip():
        markdown = "# Evidence-Aware Enrichment Interpretation\n\nNo report text was available in the analysis result."
    _markdown_to_pdf(str(markdown), out_pdf_path, title=title)


def _make_markdown_report(result: dict) -> str:
    """Return full GPT markdown if present; otherwise build deterministic pretriage markdown."""
    gpt = result.get("gpt") or {}
    parsed = gpt.get("parsed") or {}
    if isinstance(parsed, dict) and parsed.get("markdown_report"):
        return str(parsed["markdown_report"])

    raw_text = gpt.get("raw_text") or ""
    if isinstance(raw_text, str) and raw_text.strip():
        return raw_text.strip()

    display = gpt.get("display") or {}
    display_bits = []
    for key in ["headline", "executive_summary", "most_plausible_biology", "likely_artifacts_confounders", "follow_up_experiments"]:
        val = display.get(key)
        if val:
            display_bits.append(f"## {key.replace('_', ' ').title()}\n\n{val}")
    if display_bits:
        return "# Enrichment interpretation report\n\n" + "\n\n".join(display_bits)

    phenotype = result.get("phenotype", "")
    context = result.get("context") or {}
    programs = (result.get("programs") or {}).get("programs", []) or []
    triage_rows = (result.get("triage") or {}).get("rows", []) or []

    lines = [
        "# Enrichment pretriage report",
        "",
        "This deterministic pretriage report was generated without PubMed or LLM reasoning.",
        "It summarizes validation, pathway/program grouping, and ranked enrichment signals.",
        "",
        f"**Phenotype:** {phenotype or 'Not provided'}",
        "",
        "## Context",
        "",
    ]
    if isinstance(context, dict) and context:
        for k, v in context.items():
            if v:
                lines.append(f"- **{k}:** {v}")
    else:
        lines.append("- No context provided.")

    lines += ["", "## Top biological programs", ""]
    if programs:
        for i, p in enumerate(programs[:10], start=1):
            program = p.get("program") or p.get("label") or p.get("name") or "Unknown"
            score = p.get("program_score", "")
            member_count = p.get("member_count", "")
            top_genes = p.get("top_genes", []) or []
            gene_text = ", ".join(map(str, top_genes[:12]))
            lines.append(f"{i}. **{program}** — score: {score}; members: {member_count}; top genes: {gene_text}")
    else:
        lines.append("No program summaries returned.")

    lines += ["", "## Top enrichment rows", ""]
    if triage_rows:
        for i, r in enumerate(triage_rows[:15], start=1):
            term = r.get("term", "")
            score = r.get("combined_pre_gpt_score", r.get("triage_score", ""))
            padj = r.get("adjusted_p_value", "")
            overlap = r.get("overlap_k", "")
            genes = r.get("genes_list", []) or []
            flags = r.get("flags", []) or []
            lines.append(f"{i}. **{term}** — score: {score}; adj p: {padj}; overlap: {overlap}; genes: {', '.join(map(str, genes[:10]))}; flags: {', '.join(map(str, flags))}")
    else:
        lines.append("No triage rows returned.")

    lines += [
        "",
        "## Caveat",
        "",
        "This pretriage output is a prioritization and organization layer, not a causal biological conclusion. Use full mode for PubMed-aware LLM interpretation and follow-up experiment recommendations.",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------
# API routes
# ---------------------------------------------------------------------


@app.get("/health")
def health() -> dict[str, Any]:
    return {"ok": True, "service": "enrichment-llm-api", "version": APP_VERSION}


@app.get("/reports/{filename}")
def get_report(filename: str) -> FileResponse:
    path = REPORTS_DIR / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail="Report not found.")
    return FileResponse(path)


@app.post("/validate")
def validate_upload(
    file: UploadFile = File(...),
    input_format: str = Form("auto"),
) -> JSONResponse:
    df = _normalize_input_df(_read_upload_csv(file), input_format=input_format)
    validation = validate_enrichment_df(df)
    validation["normalized_columns"] = list(df.columns)
    validation["normalized_n_rows"] = int(len(df))
    return _json_response(validation, status_code=200 if validation.get("ok") else 400)


@app.post("/pretriage")
def pretriage(
    file: UploadFile = File(...),
    phenotype: str = Form(...),
    organism: str = Form(""),
    assay: str = Form(""),
    tissue: str = Form(""),
    cell_type: str = Form(""),
    perturbation: str = Form(""),
    timepoint: str = Form(""),
    extra_context_json: str | None = Form(None),
    input_format: str = Form("auto"),
) -> JSONResponse:
    """Deterministic module mode: validation + triage + program grouping. No PubMed/OpenAI call."""
    df = _normalize_input_df(_read_upload_csv(file), input_format=input_format)
    validation = validate_enrichment_df(df)
    if not validation.get("ok"):
        return _json_response({"error": validation.get("error"), "input_validation": validation}, status_code=400)

    context = _context_from_form(
        organism=organism,
        assay=assay,
        tissue=tissue,
        cell_type=cell_type,
        perturbation=perturbation,
        timepoint=timepoint,
        extra_context_json=extra_context_json,
    )

    try:
        triage = triage_enrichment_table(df, phenotype=phenotype, context=context)
        programs = summarize_programs(triage["rows"], phenotype=phenotype)
        result = {
            "mode": "pretriage_no_llm",
            "phenotype": phenotype,
            "context": context,
            "input_validation": validation,
            "triage": triage,
            "programs": programs,
        }
        result["markdown_report"] = _make_markdown_report(result)
        return _json_response(result)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.post("/analyze")
def analyze(
    file: UploadFile = File(...),
    phenotype: str = Form(...),
    organism: str = Form(""),
    assay: str = Form(""),
    tissue: str = Form(""),
    cell_type: str = Form(""),
    perturbation: str = Form(""),
    timepoint: str = Form(""),
    extra_context_json: str | None = Form(None),
    input_format: str = Form("auto"),
    mode: str = Form("full"),
) -> JSONResponse:
    """
    Full module mode by default: validation -> triage -> programs -> PubMed -> GPT -> claims.
    Use mode=pretriage for no PubMed/OpenAI call.
    """
    if (mode or "full").lower() in {"pretriage", "no_llm", "deterministic"}:
        return pretriage(
            file=file,
            phenotype=phenotype,
            organism=organism,
            assay=assay,
            tissue=tissue,
            cell_type=cell_type,
            perturbation=perturbation,
            timepoint=timepoint,
            extra_context_json=extra_context_json,
            input_format=input_format,
        )

    df = _normalize_input_df(_read_upload_csv(file), input_format=input_format)
    validation = validate_enrichment_df(df)
    if not validation.get("ok"):
        return _json_response({"error": validation.get("error"), "input_validation": validation}, status_code=400)

    context = _context_from_form(
        organism=organism,
        assay=assay,
        tissue=tissue,
        cell_type=cell_type,
        perturbation=perturbation,
        timepoint=timepoint,
        extra_context_json=extra_context_json,
    )

    try:
        from pipeline import run_enrichment_pipeline
        result = run_enrichment_pipeline(df, phenotype=phenotype, context=context)
        result["phenotype"] = phenotype
        result["context"] = context
        result["mode"] = "full_pubmed_llm"
        result["markdown_report"] = _make_markdown_report(result)
        return _json_response(result)
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail={
                "error": str(exc),
                "support_hint": "Check CSV fields, OPENAI_API_KEY, NCBI_EMAIL, and network access for PubMed/OpenAI.",
            },
        ) from exc


@app.post("/summarize")
def summarize_to_pdf(payload: Dict[str, Any]) -> dict[str, str]:
    try:
        REPORTS_DIR.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        pdf_filename = f"triage_report_{ts}.pdf"
        pdf_path = REPORTS_DIR / pdf_filename
        _build_pdf(payload, pdf_path)
        return {"pdf_url": f"/reports/{pdf_filename}", "pdf_path": str(pdf_path)}
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.post("/analyze-bundle")
def analyze_bundle(
    file: UploadFile = File(...),
    phenotype: str = Form(...),
    organism: str = Form(""),
    assay: str = Form(""),
    tissue: str = Form(""),
    cell_type: str = Form(""),
    perturbation: str = Form(""),
    timepoint: str = Form(""),
    extra_context_json: str | None = Form(None),
    input_format: str = Form("auto"),
    mode: str = Form("full"),
    make_pdf: bool = Form(True),
) -> StreamingResponse:
    """Run analysis and return a zip with JSON/CSV outputs plus optional PDF."""
    # Use the same underlying functions directly instead of calling the route functions.
    df = _normalize_input_df(_read_upload_csv(file), input_format=input_format)
    validation = validate_enrichment_df(df)
    if not validation.get("ok"):
        raise HTTPException(status_code=400, detail={"error": validation.get("error"), "input_validation": validation})

    context = _context_from_form(
        organism=organism,
        assay=assay,
        tissue=tissue,
        cell_type=cell_type,
        perturbation=perturbation,
        timepoint=timepoint,
        extra_context_json=extra_context_json,
    )

    try:
        if (mode or "full").lower() in {"pretriage", "no_llm", "deterministic"}:
            triage = triage_enrichment_table(df, phenotype=phenotype, context=context)
            programs = summarize_programs(triage["rows"], phenotype=phenotype)
            result = {
                "mode": "pretriage_no_llm",
                "phenotype": phenotype,
                "context": context,
                "input_validation": validation,
                "triage": triage,
                "programs": programs,
            }
            result["markdown_report"] = _make_markdown_report(result)
        else:
            from pipeline import run_enrichment_pipeline
            result = run_enrichment_pipeline(df, phenotype=phenotype, context=context)
            result["phenotype"] = phenotype
            result["context"] = context
            result["mode"] = "full_pubmed_llm"
            result["markdown_report"] = _make_markdown_report(result)

        tmpdir = Path(tempfile.mkdtemp(prefix="enrichment_llm_bundle_"))
        outdir = tmpdir / "bundle"
        paths = _write_result_tables(result, outdir)
        if make_pdf:
            pdf_path = outdir / "enrichment_llm_report.pdf"
            _build_pdf(result, pdf_path)
            paths["pdf_report"] = str(pdf_path)

        manifest = {
            "created_at": datetime.now().isoformat(),
            "service": "enrichment-llm-api",
            "version": APP_VERSION,
            "mode": result.get("mode"),
            "phenotype": phenotype,
            "context": context,
            "outputs": paths,
        }
        (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")

        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for p in outdir.rglob("*"):
                if p.is_file():
                    z.write(p, arcname=p.relative_to(outdir))
        buf.seek(0)
        return StreamingResponse(
            buf,
            media_type="application/zip",
            headers={"Content-Disposition": "attachment; filename=enrichment_llm_bundle.zip"},
        )
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="127.0.0.1", port=int(os.environ.get("PORT", "8002")), reload=True)

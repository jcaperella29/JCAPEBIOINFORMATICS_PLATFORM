from __future__ import annotations

import os
import re
from datetime import datetime
from typing import Any, Dict, List, Tuple, Optional
from xml.sax.saxutils import escape

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib import colors


def _safe(x: Any) -> str:
    if x is None:
        return ""
    return escape(str(x))


def _slugify(s: str) -> str:
    s = (s or "").strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "report"


def _get(d: Dict[str, Any], path: str, default=None):
    cur: Any = d
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def _join(values: List[Any], limit: int = 8) -> str:
    vals = [str(v) for v in (values or []) if v]
    return ", ".join(vals[:limit])


def _para(story: List[Any], text: str, style) -> None:
    story.append(Paragraph(_safe(text) if text else "—", style))


def _bullet_list(story: List[Any], items: List[str], body, limit: int = 8) -> None:
    if not items:
        story.append(Paragraph("—", body))
        return
    for item in items[:limit]:
        story.append(Paragraph("• " + _safe(item), body))


def build_triage_pdf(triage_json: Dict[str, Any], out_pdf_path: str, title: Optional[str] = None, subtitle: Optional[str] = None) -> None:
    _build_pdf(triage_json, out_pdf_path, title or "Evidence-Aware Enrichment Interpretation", subtitle)


def generate_pdf_from_triage_json(triage_json: Dict[str, Any], out_dir: str = "static/reports", filename_prefix: str = "triage_report") -> Tuple[str, str]:
    os.makedirs(out_dir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    phenotype = _get(triage_json, "programs.meta.phenotype") or _get(triage_json, "gpt.phenotype") or "enrichment triage"
    base = f"{filename_prefix}_{stamp}_{_slugify(phenotype)[:40]}.pdf"
    pdf_path = os.path.join(out_dir, base)
    pdf_url = f"/static/reports/{base}"
    build_triage_pdf(triage_json, pdf_path)
    return pdf_path, pdf_url



def _norm_program_name(s: str) -> str:
    """Normalize program labels so deterministic buckets and GPT-renamed claims can match."""
    s = (s or "").strip().upper()
    s = re.sub(r"[^A-Z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def _program_aliases(name: str) -> set[str]:
    """
    Return equivalent program names used by the deterministic program summarizer
    and by the GPT/role-normalizer layer.

    Example:
      program_summarizer.py may emit INFLAMMATION_NFKB.
      reasoner.py may relabel the claim as ANTI_INFLAMMATORY_NFKB_ATTENUATION
      in a dexamethasone/glucocorticoid context.
    """
    n = _norm_program_name(name)
    aliases = {n}

    nfkb_aliases = {
        "INFLAMMATION_NFKB",
        "NFKB_INFLAMMATION",
        "ANTI_INFLAMMATORY_NFKB_ATTENUATION",
        "INFLAMMATION_NFKB_ATTENUATION",
        "NFKB_ATTENUATION",
        "NF_KB_ATTENUATION",
        "ANTI_INFLAMMATORY_NF_KB_ATTENUATION",
    }
    if n in nfkb_aliases:
        aliases.update(nfkb_aliases)

    gr_aliases = {
        "GLUCOCORTICOID_GR_RESPONSE",
        "GLUCOCORTICOID_RECEPTOR_RESPONSE",
        "GR_RESPONSE",
        "NR3C1_RESPONSE",
    }
    if n in gr_aliases:
        aliases.update(gr_aliases)

    mapk_aliases = {
        "MAPK_STRESS_KINASE_ATTENUATION",
        "MAPK_STRESS_ATTENUATION",
        "STRESS_KINASE_ATTENUATION",
    }
    if n in mapk_aliases:
        aliases.update(mapk_aliases)

    ecm_aliases = {
        "ECM_FIBROSIS",
        "ECM_REMODELING",
        "EXTRACELLULAR_MATRIX",
        "EXTRACELLULAR_MATRIX_REMODELING",
    }
    if n in ecm_aliases:
        aliases.update(ecm_aliases)

    return aliases


def _find_claim_for_program(program: Dict[str, Any], claims: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Find the GPT claim that best corresponds to a deterministic program row.

    First use exact/alias label matching. If the GPT layer renamed a program, fall
    back to overlap in evidence genes/terms so the ranked table can still inherit
    the normalized role/evidence from the claim matrix.
    """
    p_name = program.get("program", "")
    p_aliases = _program_aliases(p_name)

    # 1) Exact/alias program name match.
    for c in claims or []:
        if _norm_program_name(c.get("program", "")) in p_aliases:
            return c

    # 2) Evidence overlap fallback.
    p_genes = {str(g).upper() for g in program.get("top_genes", []) or [] if g}
    p_terms = {
        str(t.get("term", "")).lower()
        for t in program.get("representative_terms", []) or []
        if isinstance(t, dict) and t.get("term")
    }

    best_claim: Dict[str, Any] = {}
    best_score = 0

    for c in claims or []:
        ev = c.get("evidence", {}) or {}
        c_genes = {str(g).upper() for g in ev.get("genes", []) or [] if g}
        c_terms = {str(t).lower() for t in ev.get("terms", []) or [] if t}

        gene_overlap = len(p_genes & c_genes)
        term_overlap = sum(
            1
            for pt in p_terms
            for ct in c_terms
            if pt and ct and (pt in ct or ct in pt)
        )

        # Gene overlap is usually stronger than term substring overlap.
        score = (2 * gene_overlap) + term_overlap
        if score > best_score:
            best_score = score
            best_claim = c

    return best_claim if best_score > 0 else {}


def _role_priority(role: str) -> int:
    x = (role or "").lower()
    if "primary driver" in x:
        return 6
    if "downstream mechanism" in x:
        return 5
    if "downstream response" in x:
        return 4
    if "reactive" in x:
        return 3
    if "uncertain" in x:
        return 2
    if "artifact" in x or "confounded" in x:
        return 1
    if "driver" in x:
        return 5  # backward compatibility with older reports
    return 0


def _evidence_priority(evidence_strength: str) -> int:
    x = (evidence_strength or "").lower()
    if "strong" in x:
        return 3
    if "moderate" in x:
        return 2
    if "weak" in x:
        return 1
    return 0


def _claim_priority(c: Dict[str, Any]) -> int:
    return 10 * _role_priority(c.get("role", "")) + _evidence_priority(c.get("evidence_strength", ""))


def _program_sort_key(program: Dict[str, Any], claims: List[Dict[str, Any]]) -> Tuple[int, float]:
    claim_for_prog = _find_claim_for_program(program, claims)
    return (_claim_priority(claim_for_prog), float(program.get("program_score", 0.0)))


def _sorted_programs_for_report(programs: List[Dict[str, Any]], claims: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return sorted(programs or [], key=lambda p: _program_sort_key(p, claims), reverse=True)


def _top_primary_driver(claims: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Pick the best top-driver KPI under the new role taxonomy."""
    if not claims:
        return None

    primary = [c for c in claims if "primary driver" in (c.get("role") or "").lower()]
    if primary:
        return sorted(primary, key=_claim_priority, reverse=True)[0]

    # Backward compatibility for older saved outputs.
    legacy = [c for c in claims if (c.get("role") or "").lower().startswith("likely driver")]
    if legacy:
        return sorted(legacy, key=_claim_priority, reverse=True)[0]

    return sorted(claims, key=_claim_priority, reverse=True)[0]


def _is_synthetic_positive_control(triage_json: Dict[str, Any]) -> bool:
    rows = ((triage_json.get("triage", {}) or {}).get("rows", []) or [])
    for r in rows:
        if "synthetic positive control" in str(r.get("term", "")).lower():
            return True
    # Also inspect claim evidence terms, in case triage rows are absent from a saved result.
    for c in triage_json.get("claims", []) or []:
        ev = c.get("evidence", {}) or {}
        for term in ev.get("terms", []) or []:
            if "synthetic positive control" in str(term).lower():
                return True
    return False


def _claim_heading(c: Dict[str, Any]) -> str:
    return str(c.get("program") or "Unassigned")


def _md_text(x: Any) -> str:
    """Plain markdown-safe-ish string. Keep content readable for machines/LLMs."""
    if x is None:
        return ""
    return str(x).replace("\r\n", "\n").replace("\r", "\n").strip()


def _md_join(values: List[Any], limit: int = 8) -> str:
    vals = [_md_text(v) for v in (values or []) if _md_text(v)]
    return ", ".join(vals[:limit]) if vals else "—"


def _md_table_cell(x: Any) -> str:
    s = _md_text(x) or "—"
    # Keep markdown tables valid while preserving meaning.
    return s.replace("|", "\\|").replace("\n", " ")


def _md_bullets(items: List[Any], limit: int = 12) -> str:
    vals = [_md_text(x) for x in (items or []) if _md_text(x)]
    if not vals:
        return "- —"
    return "\n".join(f"- {v}" for v in vals[:limit])


def build_triage_markdown(
    triage_json: Dict[str, Any],
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
) -> str:
    """
    Build a machine-friendly Markdown report from the same triage_json used by the PDF.

    This intentionally reuses the same claim/program matching helpers as the PDF so
    the Markdown report stays synchronized with the ranked program table and claim
    matrix, including GPT-normalized labels such as ANTI_INFLAMMATORY_NFKB_ATTENUATION.
    """
    report_title = title or "Evidence-Aware Enrichment Interpretation"

    gpt = triage_json.get("gpt", {}) or {}
    parsed = gpt.get("parsed", {}) or {}
    display = gpt.get("display", {}) or triage_json.get("gpt_display", {}) or {}
    claims = triage_json.get("claims") or []
    programs = (triage_json.get("programs", {}) or {}).get("programs", []) or []
    pubmed = triage_json.get("pubmed", {}) or {}
    metadata = triage_json.get("metadata", {}) or {}
    context = triage_json.get("context") or gpt.get("experiment_context") or {}
    synthetic_positive_control = _is_synthetic_positive_control(triage_json)
    programs_for_report = _sorted_programs_for_report(programs, claims)
    claims_for_report = sorted(claims, key=_claim_priority, reverse=True)

    executive_summary = (
        parsed.get("executive_summary")
        or display.get("executive_summary")
        or display.get("gpt_summary")
        or display.get("headline")
        or gpt.get("raw_text", "")[:900]
        or "—"
    )

    top_driver = _top_primary_driver(claims)
    biggest_confounder = ""
    for c in claims:
        if c.get("confounders"):
            biggest_confounder = c["confounders"][0]
            break

    top_exp = parsed.get("next_best_experiment") or ""
    if not top_exp:
        for c in claims:
            vals = c.get("validation") or []
            if vals:
                v = vals[0]
                top_exp = f"{v.get('experiment', '')}; readout: {v.get('readout', '')}; control: {v.get('control', '')}"
                break

    lines: List[str] = []
    lines.append(f"# {_md_text(report_title)}")
    lines.append("")
    lines.append(_md_text(subtitle or f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"))
    lines.append("")

    lines.append("## Executive Summary")
    lines.append("")
    lines.append(_md_text(executive_summary))
    lines.append("")
    if synthetic_positive_control:
        lines.append("> **Demo input notice:** This run contains terms labeled synthetic positive control. Use this report to test product behavior and presentation, not as biological evidence.")
        lines.append("")

    lines.append("## Key Takeaways")
    lines.append("")
    lines.append(f"- **Top interpretation:** {_md_text(top_driver.get('claim', '—') if top_driver else '—')}")
    lines.append(f"- **Top driver/program:** {_md_text(top_driver.get('program', '—') if top_driver else '—')}")
    lines.append(f"- **Biggest confounder:** {_md_text(biggest_confounder or '—')}")
    lines.append(f"- **Recommended next experiment:** {_md_text(top_exp or '—')}")
    lines.append("")

    lines.append("## Run Metadata")
    lines.append("")
    meta_lines = [
        f"App version: {metadata.get('app_version', '—')}",
        f"Prompt version: {metadata.get('prompt_version', '—')}",
        f"Playbook version: {metadata.get('playbook_version', '—')}",
        f"Model: {metadata.get('model_name', '—')}",
        f"Input hash: {metadata.get('input_hash', '—')}",
    ]
    lines.append(_md_bullets(meta_lines, limit=10))
    lines.append("")

    lines.append("## Ranked Program Table")
    lines.append("")
    if programs_for_report:
        lines.append("| Program | Score | Role | Evidence | Key genes | Supporting terms |")
        lines.append("|---|---:|---|---|---|---|")
        for p in programs_for_report[:10]:
            claim_for_prog = _find_claim_for_program(p, claims)
            display_program = claim_for_prog.get("program") or p.get("program", "")
            terms = [x.get("term", "") for x in p.get("representative_terms", [])[:3]]
            lines.append(
                "| "
                + " | ".join([
                    _md_table_cell(display_program),
                    _md_table_cell(f"{float(p.get('program_score', 0.0)):.1f}"),
                    _md_table_cell(claim_for_prog.get("role", "—")),
                    _md_table_cell(claim_for_prog.get("evidence_strength", "—")),
                    _md_table_cell(_md_join(p.get("top_genes", []), 8)),
                    _md_table_cell(_md_join(terms, 3)),
                ])
                + " |"
            )
    else:
        lines.append("No program summaries available.")
    lines.append("")

    lines.append("## Claim-Evidence Matrix")
    lines.append("")
    if claims_for_report:
        for idx, c in enumerate(claims_for_report[:8], 1):
            heading = _claim_heading(c)
            lines.append(f"### {idx}. {_md_text(heading)}: {_md_text(c.get('role', 'Uncertain'))} / {_md_text(c.get('evidence_strength', 'Weak'))}")
            lines.append("")
            lines.append(_md_text(c.get("claim", "—")))
            lines.append("")
            if c.get("rationale"):
                lines.append("**Rationale**")
                lines.append("")
                lines.append(_md_text(c.get("rationale", "")))
                lines.append("")
            ev = c.get("evidence", {}) or {}
            lines.append("**Evidence**")
            lines.append("")
            lines.append(f"- **Terms:** {_md_join(ev.get('terms', []), 6)}")
            lines.append(f"- **Genes:** {_md_join(ev.get('genes', []), 10)}")
            lines.append(f"- **PMIDs:** {_md_join(ev.get('pmids', []), 6)}")
            lines.append(f"- **Literature status:** {_md_text(ev.get('literature_status', 'not_assessed'))}")
            lines.append("")
    else:
        lines.append(_md_text(gpt.get("raw_text", "No structured claims available.")))
        lines.append("")

    lines.append("## Confounders and Assay Limitations")
    lines.append("")
    lines.append("### Assay/context")
    lines.append("")
    ctx_bits = [f"{k}: {v}" for k, v in (context or {}).items() if v]
    lines.append(_md_bullets(ctx_bits, limit=10))
    lines.append("")

    lines.append("### Main confounders")
    lines.append("")
    confounders = parsed.get("main_confounders") or []
    if not confounders:
        for c in claims:
            for x in c.get("confounders", []) or []:
                if x not in confounders:
                    confounders.append(x)
    lines.append(_md_bullets(confounders, limit=12))
    lines.append("")

    lines.append("### Assay limitations")
    lines.append("")
    lines.append(_md_bullets(parsed.get("assay_limitations", []), limit=12))
    lines.append("")

    lines.append("## Validation Plan")
    lines.append("")
    seen_validation_headings = set()
    wrote_validation = False
    for c in claims_for_report[:8]:
        vals = c.get("validation", []) or []
        if not vals:
            continue
        heading = _claim_heading(c)
        if heading in seen_validation_headings:
            continue
        seen_validation_headings.add(heading)
        wrote_validation = True
        lines.append(f"### {_md_text(heading)}")
        lines.append("")
        for v in vals[:3]:
            lines.append(f"- **Experiment:** {_md_text(v.get('experiment', ''))}")
            lines.append(f"  - **Readout:** {_md_text(v.get('readout', ''))}")
            lines.append(f"  - **Control:** {_md_text(v.get('control', ''))}")
            lines.append(f"  - **Expected if causal:** {_md_text(v.get('expected_result_if_causal', ''))}")
        lines.append("")
    if not wrote_validation:
        lines.append("—")
        lines.append("")

    lines.append("## Literature Context (PubMed)")
    lines.append("")
    papers = pubmed.get("papers", []) or []
    if not papers:
        lines.append("No PubMed papers available for this run.")
    else:
        for i, p in enumerate(papers[:5], start=1):
            pmid = _md_text(p.get("pmid", ""))
            title_p = _md_text(p.get("title", "Untitled")) or "Untitled"
            source = _md_text(p.get("source", ""))
            pubdate = _md_text(p.get("pubdate", ""))
            lines.append(f"{i}. **{title_p}**")
            lines.append(f"   - PMID: {pmid or '—'} | {source or '—'} | {pubdate or '—'}")
    lines.append("")

    return "\n".join(lines).strip() + "\n"


def write_triage_markdown(
    triage_json: Dict[str, Any],
    out_md_path: str,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
) -> None:
    """Write the Markdown report to disk."""
    md = build_triage_markdown(triage_json, title=title, subtitle=subtitle)
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write(md)


def _build_pdf(triage_json: Dict[str, Any], pdf_path: str, title: str, subtitle: Optional[str] = None) -> None:
    styles = getSampleStyleSheet()
    h1 = styles["Heading1"]
    h2 = styles["Heading2"]
    h3 = styles["Heading3"]
    body = styles["BodyText"]
    small = ParagraphStyle("Small", parent=body, fontSize=8, leading=10)

    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=letter,
        leftMargin=0.65 * inch,
        rightMargin=0.65 * inch,
        topMargin=0.65 * inch,
        bottomMargin=0.65 * inch,
        title=title,
    )

    story: List[Any] = []
    gpt = triage_json.get("gpt", {}) or {}
    parsed = gpt.get("parsed", {}) or {}
    display = gpt.get("display", {}) or triage_json.get("gpt_display", {}) or {}
    claims = triage_json.get("claims") or []
    programs = (triage_json.get("programs", {}) or {}).get("programs", []) or []
    pubmed = triage_json.get("pubmed", {}) or {}
    metadata = triage_json.get("metadata", {}) or {}
    context = triage_json.get("context") or gpt.get("experiment_context") or {}
    synthetic_positive_control = _is_synthetic_positive_control(triage_json)
    programs_for_report = _sorted_programs_for_report(programs, claims)

    # Page 1: executive summary
    story.append(Paragraph(_safe(title), h1))
    story.append(Paragraph(_safe(subtitle or datetime.now().strftime("%Y-%m-%d %H:%M:%S")), body))
    story.append(Spacer(1, 0.18 * inch))

    story.append(Paragraph("Executive summary", h2))
    _para(story, parsed.get("executive_summary") or display.get("executive_summary") or display.get("gpt_summary") or display.get("headline") or gpt.get("raw_text", "")[:900], body)
    if synthetic_positive_control:
        story.append(Spacer(1, 0.08 * inch))
        story.append(Paragraph("<b>Demo input notice:</b> This run contains terms labeled synthetic positive control. Use this report to test product behavior and presentation, not as biological evidence.", body))
    story.append(Spacer(1, 0.15 * inch))

    top_driver = _top_primary_driver(claims)
    biggest_confounder = ""
    for c in claims:
        if c.get("confounders"):
            biggest_confounder = c["confounders"][0]
            break
    top_exp = parsed.get("next_best_experiment") or ""
    if not top_exp:
        for c in claims:
            vals = c.get("validation") or []
            if vals:
                v = vals[0]
                top_exp = f"{v.get('experiment', '')}; readout: {v.get('readout', '')}; control: {v.get('control', '')}"
                break

    kpi_data = [
        ["Top interpretation", top_driver.get("claim", "—") if top_driver else "—"],
        ["Top driver/program", top_driver.get("program", "—") if top_driver else "—"],
        ["Biggest confounder", biggest_confounder or "—"],
        ["Recommended next experiment", top_exp or "—"],
    ]
    t = Table([[Paragraph(_safe(a), body), Paragraph(_safe(b), body)] for a, b in kpi_data], colWidths=[1.8 * inch, 4.8 * inch])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#eef2ff")),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
    ]))
    story.append(t)

    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph("Run metadata", h3))
    meta_lines = [
        f"App version: {metadata.get('app_version', '—')}",
        f"Prompt version: {metadata.get('prompt_version', '—')}",
        f"Playbook version: {metadata.get('playbook_version', '—')}",
        f"Model: {metadata.get('model_name', '—')}",
        f"Input hash: {metadata.get('input_hash', '—')}",
    ]
    _bullet_list(story, meta_lines, small, limit=10)
    story.append(PageBreak())

    # Page 2: ranked programs
    story.append(Paragraph("Ranked program table", h2))
    if programs_for_report:
        rows = [["Program", "Score", "Role", "Evidence", "Key genes", "Supporting terms"]]
        for p in programs_for_report[:10]:
            claim_for_prog = _find_claim_for_program(p, claims)
            # Prefer the GPT-normalized claim label when available. This keeps
            # the ranked table synchronized with the claim matrix.
            display_program = claim_for_prog.get("program") or p.get("program", "")
            terms = [x.get("term", "") for x in p.get("representative_terms", [])[:3]]
            rows.append([
                display_program,
                f"{float(p.get('program_score', 0.0)):.1f}",
                claim_for_prog.get("role", "—"),
                claim_for_prog.get("evidence_strength", "—"),
                _join(p.get("top_genes", []), 8),
                _join(terms, 3),
            ])
        table = Table([[Paragraph(_safe(str(cell)), small if i else body) for cell in row] for i, row in enumerate(rows)], colWidths=[1.35*inch, .55*inch, 1.0*inch, .75*inch, 1.5*inch, 1.6*inch])
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dbeafe")),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ]))
        story.append(table)
    else:
        story.append(Paragraph("No program summaries available.", body))
    story.append(PageBreak())

    # Page 3: claim-evidence matrix
    story.append(Paragraph("Claim-evidence matrix", h2))
    claims_for_report = sorted(claims, key=_claim_priority, reverse=True)
    if claims_for_report:
        for idx, c in enumerate(claims_for_report[:8], 1):
            story.append(Paragraph(f"{idx}. {_safe(_claim_heading(c))}: {_safe(c.get('role', 'Uncertain'))} / {_safe(c.get('evidence_strength', 'Weak'))}", h3))
            _para(story, c.get("claim", ""), body)
            if c.get("rationale"):
                story.append(Paragraph("<b>Rationale</b>", body))
                _para(story, c.get("rationale", ""), body)
            ev = c.get("evidence", {}) or {}
            ev_rows = [
                ["Terms", _join(ev.get("terms", []), 6)],
                ["Genes", _join(ev.get("genes", []), 10)],
                ["PMIDs", _join(ev.get("pmids", []), 6)],
                ["Literature status", ev.get("literature_status", "not_assessed")],
            ]
            ev_table = Table([[Paragraph(_safe(a), small), Paragraph(_safe(b), small)] for a, b in ev_rows], colWidths=[1.3*inch, 5.3*inch])
            ev_table.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
            story.append(ev_table)
            story.append(Spacer(1, 0.12 * inch))
    else:
        _para(story, gpt.get("raw_text", "No structured claims available."), body)
    story.append(PageBreak())

    # Page 4: confounders and assay limits
    story.append(Paragraph("Confounders and assay limitations", h2))
    story.append(Paragraph("Assay/context", h3))
    ctx_bits = [f"{k}: {v}" for k, v in context.items() if v]
    _bullet_list(story, ctx_bits, body, limit=10)
    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph("Main confounders", h3))
    confounders = parsed.get("main_confounders") or []
    if not confounders:
        for c in claims:
            for x in c.get("confounders", []) or []:
                if x not in confounders:
                    confounders.append(x)
    _bullet_list(story, confounders, body, limit=12)
    story.append(Paragraph("Assay limitations", h3))
    _bullet_list(story, parsed.get("assay_limitations", []), body, limit=12)
    story.append(PageBreak())

    # Page 5: validation plan + PubMed
    story.append(Paragraph("Validation plan", h2))
    seen_validation_headings = set()
    for c in claims_for_report[:8]:
        vals = c.get("validation", []) or []
        if not vals:
            continue
        heading = _claim_heading(c)
        # Avoid duplicate validation sections when GPT repeats a program. Alternative explanations
        # keep their explicit alternative label, so they remain visible without seeming duplicated.
        if heading in seen_validation_headings:
            continue
        seen_validation_headings.add(heading)
        story.append(Paragraph(_safe(heading), h3))
        for v in vals[:3]:
            _para(story, f"Experiment: {v.get('experiment', '')}; Readout: {v.get('readout', '')}; Control: {v.get('control', '')}; Expected if causal: {v.get('expected_result_if_causal', '')}", body)
    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph("Literature context (PubMed)", h2))
    papers = pubmed.get("papers", []) or []
    if not papers:
        story.append(Paragraph("No PubMed papers available for this run.", body))
    else:
        for i, p in enumerate(papers[:5], start=1):
            pmid = p.get("pmid", "") or ""
            title_p = p.get("title", "") or "Untitled"
            source = p.get("source", "") or ""
            pubdate = p.get("pubdate", "") or ""
            story.append(Paragraph(f"<b>{i}. {_safe(title_p)}</b>", body))
            story.append(Paragraph(_safe(f"PMID: {pmid} | {source} | {pubdate}"), small))
            story.append(Spacer(1, 0.08 * inch))

    doc.build(story)

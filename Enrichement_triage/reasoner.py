from __future__ import annotations

import os
import json
import re
from pathlib import Path
from typing import Dict, Any, Optional, List

from openai import OpenAI
from schemas import claim_schema_for_prompt
from playbook_retriever import retrieve_playbook_context

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

SYSTEM = """You are a senior computational biologist.
Your job is to convert enrichment results into cautious, auditable biological interpretation claims.
Be strict about causal vs reactive vs artifact/confounders.
Propose follow-up experiments with concrete readouts and controls.
Use cautious, evidence-weighted language.
Never claim causality from enrichment alone.
When external literature evidence is provided, use it carefully:
- treat it as supporting context, not automatic proof
- note when literature aligns with the enrichment results
- note when literature and the current data do not align
- do not overclaim causality from literature alone
- when you refer to a retrieved paper, cite it inline as (PMID: XXXXXXXX)
- respect literature relevance labels: weak_background papers are background only, general_support papers are context, and direct_support_candidate papers are still not causal proof
"""

# Works when reasoner.py lives at project root; also tolerates older nested layout.
_HERE = Path(__file__).resolve().parent
PLAYBOOK_DIR = _HERE / "playbook"
if not PLAYBOOK_DIR.exists():
    PLAYBOOK_DIR = _HERE.parent / "playbook"

GLOBAL_PLAYBOOKS = [
    "16_evidence_weighting_and_translation.md",   # preferred production name
    "6_evidence_weighting_and_translation.md",    # tolerated older/uploaded name
    "03_growth_axis_and_overlap_rules.md",
    "08_tissue_celltype_prior_map.md",
    "09_followup_experiment_menu.md",
]

ASSAY_TO_PLAYBOOK = {
    "bulk_rnaseq": ["01_assay_confounders_rnaseq.md"],
    "scrnaseq": ["02_assay_confounders_scrna.md"],
    "perturbseq": ["10_assay_confounders_perturbseq.md"],
    "atacseq": ["11_assay_confounders_atacseq.md"],
    "mirnaseq": ["12_assay_confounders_mirnaseq.md"],
    "gwas": ["13_assay_confounders_gwas.md"],
    "dna_methylation": [
        "14_epigenetic_vs_transcriptional_priors.md",
        "15_assay_confounders_dna_methylation.md",
    ],
}

SECTION_ALIASES = {
    "headline": "headline",
    "executive summary": "executive_summary",
    "experimental context": "experimental_context",
    "most plausible biology": "most_plausible_biology",
    "likely reactive programs": "likely_reactive_programs",
    "likely artifacts / confounders": "likely_artifacts_confounders",
    "likely artifacts/confounders": "likely_artifacts_confounders",
    "evidence strength and rationale": "evidence_strength_rationale",
    "follow-up experiments": "follow_up_experiments",
    "main uncertainties": "main_uncertainties",
    "literature context": "literature_context",
}


def _norm_assay(s: str) -> str:
    s = (s or "").strip().lower()
    s = s.replace("-", "").replace(" ", "").replace("_", "")
    if s in ("bulkrnaseq", "rnaseq"):
        return "bulk_rnaseq"
    if s in ("scrnaseq", "singlecellrnaseq"):
        return "scrnaseq"
    if s in ("perturbseq",):
        return "perturbseq"
    if s in ("atacseq",):
        return "atacseq"
    if s in ("mirnaseq",):
        return "mirnaseq"
    if s in ("gwas",):
        return "gwas"
    if s in ("dnamethylation", "methylation"):
        return "dna_methylation"
    return s


def _load_md_files(files: List[str]) -> str:
    chunks: List[str] = []
    seen_paths = set()
    for fn in files:
        p = PLAYBOOK_DIR / fn
        if p.exists() and p not in seen_paths:
            chunks.append(f"# {fn}\n\n" + p.read_text(encoding="utf-8"))
            seen_paths.add(p)
    return "\n\n---\n\n".join(chunks).strip()


def _load_playbook_md(assay: str) -> str:
    assay_key = _norm_assay(assay)
    assay_files = ASSAY_TO_PLAYBOOK.get(assay_key, [])
    all_files = GLOBAL_PLAYBOOKS + assay_files
    return _load_md_files(all_files)


def parse_gpt_markdown_sections(text: str) -> Dict[str, str]:
    if not text or not text.strip():
        return {}
    pattern = re.compile(r"(?ms)^##\s+(.+?)\s*$(.*?)(?=^##\s+.+?$|\Z)")
    sections: Dict[str, str] = {}
    for match in pattern.finditer(text):
        raw_heading = match.group(1).strip()
        body = match.group(2).strip()
        norm = raw_heading.lower()
        key = SECTION_ALIASES.get(norm, re.sub(r"[^a-z0-9]+", "_", norm).strip("_"))
        sections[key] = body
    return sections


def build_gpt_display_fields(gpt: Dict[str, Any]) -> Dict[str, str]:
    parsed = gpt.get("parsed", {}) or {}
    sections = gpt.get("sections", {}) or {}
    claims = parsed.get("claims", []) or []

    ranked = []
    confs = []
    validations = []
    top_primary_driver = ""
    top_downstream_mechanism = ""
    for c in claims:
        role_l = (c.get('role') or '').lower()
        if not top_primary_driver and role_l.startswith('likely primary driver'):
            top_primary_driver = c.get('program') or c.get('claim') or ''
        if not top_downstream_mechanism and role_l.startswith('likely downstream'):
            top_downstream_mechanism = c.get('program') or c.get('claim') or ''
        ranked.append(f"{c.get('program', '')}: {c.get('role', 'Uncertain')} / {c.get('evidence_strength', 'Weak')} — {c.get('claim', '')}")
        for x in c.get("confounders", []) or []:
            if x not in confs:
                confs.append(x)
        for v in c.get("validation", []) or []:
            exp = v.get("experiment", "")
            readout = v.get("readout", "")
            control = v.get("control", "")
            if exp:
                validations.append(f"{exp}; readout: {readout}; control: {control}")

    return {
        "headline": parsed.get("headline") or sections.get("headline", ""),
        "experimental_context": sections.get("experimental_context", ""),
        "executive_summary": parsed.get("executive_summary", ""),
        "most_plausible_biology": "\n".join(ranked) or sections.get("most_plausible_biology", ""),
        "likely_reactive_programs": sections.get("likely_reactive_programs", ""),
        "likely_artifacts_confounders": "\n".join(confs) or sections.get("likely_artifacts_confounders", ""),
        "evidence_strength_rationale": sections.get("evidence_strength_rationale", ""),
        "follow_up_experiments": "\n".join(validations) or sections.get("follow_up_experiments", ""),
        "main_uncertainties": "\n".join(parsed.get("assay_limitations", []) or []) or sections.get("main_uncertainties", ""),
        "top_primary_driver": top_primary_driver,
        "top_downstream_mechanism": top_downstream_mechanism,
        "raw_text": gpt.get("raw_text", ""),
    }


def _compact_pubmed_context(pubmed_context: Optional[Dict[str, Any]], max_papers: int = 5) -> Dict[str, Any]:
    if not isinstance(pubmed_context, dict):
        return {}
    papers = pubmed_context.get("papers", []) or []
    compact_papers = []
    for p in papers[:max_papers]:
        if not isinstance(p, dict):
            continue
        compact_papers.append({
            "pmid": str(p.get("pmid", "")),
            "title": p.get("title", ""),
            "pubdate": p.get("pubdate", ""),
            "source": p.get("source", ""),
            "authors": (p.get("authors", []) or [])[:5],
            "abstract": (p.get("abstract", "") or "")[:2500],
            "url": p.get("url", ""),
            "relevance_label": p.get("relevance_label", "not_assessed"),
            "relevance_score": p.get("relevance_score", 0),
            "relevance_reason": p.get("relevance_reason", ""),
            "literature_use": p.get("literature_use", ""),
        })
    return {
        "status": pubmed_context.get("status", ""),
        "source": pubmed_context.get("source", "PubMed via NCBI E-utilities"),
        "query": pubmed_context.get("query", ""),
        "query_used": pubmed_context.get("query_used", ""),
        "query_strategy": pubmed_context.get("query_strategy", ""),
        "retrieval_quality": pubmed_context.get("retrieval_quality", "not_assessed"),
        "retrieval_quality_reason": pubmed_context.get("retrieval_quality_reason", ""),
        "relevance_counts": pubmed_context.get("relevance_counts", {}),
        "literature_use_guidance": pubmed_context.get("literature_use_guidance", ""),
        "top_terms_used": pubmed_context.get("top_terms_used", []),
        "top_genes_used": pubmed_context.get("top_genes_used", []),
        "papers": compact_papers,
        "error": pubmed_context.get("error", ""),
    }


def _extract_json_object(text: str) -> Dict[str, Any]:
    try:
        return json.loads(text)
    except Exception:
        pass
    m = re.search(r"\{.*\}", text or "", flags=re.S)
    if not m:
        raise ValueError("No JSON object found in model response.")
    return json.loads(m.group(0))


def _markdown_from_parsed(parsed: Dict[str, Any]) -> str:
    if parsed.get("markdown_report"):
        return parsed["markdown_report"]
    lines = ["## Headline", parsed.get("headline", ""), "", "## Executive Summary", parsed.get("executive_summary", "")]
    lines.append("\n## Claims")
    for c in parsed.get("claims", []) or []:
        lines.append(f"- **{c.get('program', '')}** [{c.get('role', 'Uncertain')} / {c.get('evidence_strength', 'Weak')}]: {c.get('claim', '')}")
    lines.append("\n## Main Confounders")
    for x in parsed.get("main_confounders", []) or []:
        lines.append(f"- {x}")
    lines.append("\n## Next Best Experiment")
    lines.append(parsed.get("next_best_experiment", ""))
    return "\n".join(lines).strip()


ROLE_VOCABULARY = [
    "Likely primary driver",
    "Likely downstream mechanism",
    "Likely downstream response",
    "Likely reactive",
    "Likely artifact / confounded",
    "Uncertain",
]


def _text_has_any(text: str, terms: List[str]) -> bool:
    text = (text or "").lower()
    return any(t.lower() in text for t in terms)


def normalize_program_label(program: str, phenotype: str, context: Dict[str, Any]) -> str:
    """
    Make program labels clearer for reporting without changing the underlying biology.
    This mainly prevents generic labels like INFLAMMATION_NFKB from reading as
    "inflammation is the primary driver" in glucocorticoid/dexamethasone runs.
    """
    p = (program or "").strip()
    p_l = p.lower().replace("-", "_")
    perturbation_l = str((context or {}).get("perturbation", "")).lower()
    phenotype_l = (phenotype or "").lower()
    combined = " ".join([p_l, perturbation_l, phenotype_l])

    is_glucocorticoid_context = _text_has_any(
        combined,
        ["dexamethasone", "glucocorticoid", "corticosteroid", "steroid"],
    )

    if is_glucocorticoid_context:
        if _text_has_any(p_l, ["inflammation", "nfkb", "nf_kappa", "nfκb", "cytokine", "chemokine"]):
            return "ANTI_INFLAMMATORY_NFKB_ATTENUATION"
        if _text_has_any(p_l, ["mapk", "stress_kinase", "stress kinase"]):
            return "MAPK_STRESS_KINASE_ATTENUATION"

    return p


def normalize_interpretation_role(
    role: str,
    program: str,
    claim: str,
    phenotype: str,
    context: Dict[str, Any],
) -> str:
    """
    Normalize model-produced role labels so downstream effects are not overcalled
    as primary drivers. This is intentionally conservative and mostly affects
    display/reporting. It is safe even if the JSON schema still uses the older
    role labels, because it runs after the model response is parsed.
    """
    role = (role or "").strip()
    role_l = role.lower()
    program_l = (program or "").lower().replace("-", "_")
    claim_l = (claim or "").lower()
    phenotype_l = (phenotype or "").lower()
    perturbation_l = str((context or {}).get("perturbation", "")).lower()
    text = " ".join([program_l, claim_l, phenotype_l, perturbation_l])

    # Preserve explicit negative/uncertain calls.
    if "artifact" in role_l or "confounded" in role_l:
        return "Likely artifact / confounded"
    if "uncertain" in role_l:
        return "Uncertain"

    # Already-normalized labels pass through.
    if role in ROLE_VOCABULARY:
        return role

    has_glucocorticoid_perturbation = _text_has_any(
        perturbation_l,
        ["dexamethasone", "glucocorticoid", "corticosteroid", "steroid"],
    )
    is_gr_program = _text_has_any(
        text,
        ["glucocorticoid", "nr3c1", "gr_response", "gr response", "steroid hormone", "fkbp5", "tsc22d3"],
    )
    is_nfkb_or_inflammation = _text_has_any(
        text,
        ["nfkb", "nf-kappa", "nf_kappa", "nfκb", "inflammation", "cytokine", "chemokine", "tnf", "il6", "cxcl8", "ccl2"],
    )
    is_mapk_stress = _text_has_any(
        text,
        ["mapk", "stress kinase", "dusp1", "dusp5", "dusp10", "p38", "jnk", "erk"],
    )

    # Dexamethasone/glucocorticoid runs: GR is upstream; NF-kB/MAPK/cytokine
    # outputs are generally downstream unless independent evidence says otherwise.
    if has_glucocorticoid_perturbation:
        if is_gr_program:
            return "Likely primary driver"
        if is_nfkb_or_inflammation or is_mapk_stress:
            if "reactive" in role_l or _text_has_any(claim_l, ["late", "reactive", "cytokine output", "chemokine output"]):
                return "Likely downstream response"
            return "Likely downstream mechanism"

    # General language-based correction.
    if _text_has_any(claim_l, ["downstream of", "secondary to", "consequence of", "attenuation downstream"]):
        return "Likely downstream mechanism"
    if _text_has_any(claim_l, ["late timepoint", "reactive consequence", "secondary response", "endpoint"]):
        return "Likely downstream response"

    # Backward compatibility with older prompt/schema labels.
    if role_l in {"likely driver", "driver"}:
        return "Likely primary driver"
    if "reactive" in role_l:
        return "Likely reactive"

    return role or "Uncertain"


def normalize_claim_roles_and_labels(parsed: Dict[str, Any], phenotype: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """Apply role and label normalization in-place and return parsed."""
    if not isinstance(parsed, dict):
        return parsed

    for c in parsed.get("claims", []) or []:
        if not isinstance(c, dict):
            continue
        old_program = c.get("program", "")
        new_program = normalize_program_label(old_program, phenotype, context)
        c["program"] = new_program
        c["role"] = normalize_interpretation_role(
            role=c.get("role", ""),
            program=new_program or old_program,
            claim=c.get("claim", ""),
            phenotype=phenotype,
            context=context,
        )

    # Keep top-level summaries aligned if they exist.
    if parsed.get("top_driver"):
        parsed["top_driver"] = normalize_program_label(str(parsed.get("top_driver", "")), phenotype, context)

    return parsed


def gpt5_reason_simple(
    *,
    phenotype: str,
    context: Dict[str, Any],
    triage: Dict[str, Any],
    programs: Dict[str, Any],
    pubmed_context: Optional[Dict[str, Any]] = None,
    playbook_context: Optional[Dict[str, Any]] = None,
    vector_store_id: Optional[str] = None,
    model: str = "gpt-5",
) -> Dict[str, Any]:
    vs_id = vector_store_id or os.environ.get("VECTOR_STORE_ID")

    payload = {
        "phenotype": phenotype,
        "experiment_context": context,
        "top_programs": (programs.get("programs") or [])[:12],
        "top_terms": (triage.get("rows") or [])[:50],
    }

    # Retrieval-first playbook grounding.
    # This replaces the earlier direct injection of full selected markdown playbooks.
    # If no vector store is configured, retrieve_playbook_context uses a compact local fallback.
    if playbook_context is None:
        playbook_context = retrieve_playbook_context(
            phenotype=phenotype,
            context=context,
            triage=triage,
            programs=programs,
            vector_store_id=vs_id,
        )

    playbook_guidance = (playbook_context or {}).get("context_text", "").strip()
    pubmed_payload = _compact_pubmed_context(pubmed_context)
    schema = claim_schema_for_prompt()

    prompt = f"""
Return ONLY valid JSON matching this schema. Do not wrap it in markdown.

JSON schema:
{json.dumps(schema, indent=2)}

Experiment context:
{json.dumps(context, indent=2)}

Phenotype:
{phenotype}

RETRIEVED PLAYBOOK GUIDANCE (authoritative; follow these when relevant):
{playbook_guidance if playbook_guidance else "(none retrieved; rely on required interpretation behavior)"}

External biomedical literature context (PubMed / NCBI):
{json.dumps(pubmed_payload, indent=2)}

Enrichment summary:
{json.dumps(payload, indent=2)}

Required interpretation behavior:
- Create 3 to 8 InterpretationClaim-style claims.
- Tie every claim to supporting terms, row IDs, genes, and PMIDs when available.
- Use paper-level relevance_label values: weak_background papers must not be used as direct support; general_support papers can support context; direct_support_candidate papers can be cited as direct literature support but still do not prove causality.
- If retrieval_quality is low or low_to_moderate, say the PubMed section is background retrieval and keep literature_status as background/general_support unless a specific paper clearly overlaps the claim.
- Classify each program using one of: Likely primary driver, Likely downstream mechanism, Likely downstream response, Likely reactive, Likely artifact / confounded, or Uncertain.
- Use Likely primary driver only for the most upstream perturbation-linked or phenotype-linked program.
- Use Likely downstream mechanism for pathways mechanistically downstream of the primary driver but still biologically important.
- Use Likely downstream response for late transcriptional outputs, cytokine/chemokine changes, compensation, or endpoint biology.
- Do not label both an upstream regulator and its downstream consequence as primary drivers unless independent evidence supports both.
- Assign evidence_strength as Stronger, Moderate, or Weak.
- Every claim needs at least one validation experiment with a readout and control.
- Include assay limitations and top confounders.
- Use cautious language: consistent with, suggestive of, cannot distinguish from, requires validation.
""".strip()

    messages = [
        {"role": "system", "content": SYSTEM},
        {"role": "user", "content": prompt},
    ]

    tools = []
    if vs_id:
        tools = [{"type": "file_search", "vector_store_ids": [vs_id]}]

    # Prefer JSON mode when available; fall back to text for older SDK/model combinations.
    try:
        resp = client.responses.create(
            model=model,
            input=messages,
            tools=tools,
            text={"format": {"type": "json_schema", "name": "enrichment_interpretation", "schema": schema, "strict": False}},
        )
    except Exception:
        resp = client.responses.create(
            model=model,
            input=messages,
            tools=tools,
            text={"format": {"type": "text"}},
        )

    out = getattr(resp, "output_text", None)
    if not out:
        raise RuntimeError(f"No output_text returned. Raw response: {resp}")

    try:
        parsed = _extract_json_object(out)
        parsed = normalize_claim_roles_and_labels(parsed, phenotype=phenotype, context=context)
        raw_text = _markdown_from_parsed(parsed)
        sections = parse_gpt_markdown_sections(raw_text)
        gpt_result = {
            "model": model,
            "raw_text": raw_text,
            "raw_json_text": out,
            "parsed": parsed,
            "sections": sections,
            "phenotype": phenotype,
            "experiment_context": context,
            "pubmed_context": pubmed_payload,
            "playbook_context": playbook_context or {},
            "playbook_retrieval_mode": (playbook_context or {}).get("mode", "none"),
        }
    except Exception:
        sections = parse_gpt_markdown_sections(out)
        gpt_result = {
            "model": model,
            "raw_text": out,
            "sections": sections,
            "parsed": {},
            "phenotype": phenotype,
            "experiment_context": context,
            "pubmed_context": pubmed_payload,
            "playbook_context": playbook_context or {},
            "playbook_retrieval_mode": (playbook_context or {}).get("mode", "none"),
        }

    gpt_result["display"] = build_gpt_display_fields(gpt_result)
    return gpt_result

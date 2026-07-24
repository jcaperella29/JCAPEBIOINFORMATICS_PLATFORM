from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from openai import OpenAI

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

_HERE = Path(__file__).resolve().parent
PLAYBOOK_DIR = _HERE / "playbook"
if not PLAYBOOK_DIR.exists():
    PLAYBOOK_DIR = _HERE.parent / "playbook"

GLOBAL_PLAYBOOKS = [
    "16_evidence_weighting_and_translation.md",
    "6_evidence_weighting_and_translation.md",
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


def _norm_text(x: Any) -> str:
    return "" if x is None else str(x).strip()


def _norm_assay(s: str) -> str:
    s = (s or "").strip().lower()
    s = s.replace("-", "").replace(" ", "").replace("_", "")
    if s in ("bulkrnaseq", "rnaseq"):
        return "bulk_rnaseq"
    if s in ("scrnaseq", "singlecellrnaseq"):
        return "scrnaseq"
    if s == "perturbseq":
        return "perturbseq"
    if s == "atacseq":
        return "atacseq"
    if s == "mirnaseq":
        return "mirnaseq"
    if s == "gwas":
        return "gwas"
    if s in ("dnamethylation", "methylation"):
        return "dna_methylation"
    return s


def _safe_list(x: Any) -> List[Any]:
    return x if isinstance(x, list) else []


def _extract_program_names(programs: Dict[str, Any], max_items: int = 8) -> List[str]:
    out: List[str] = []
    for p in _safe_list((programs or {}).get("programs"))[:max_items]:
        if not isinstance(p, dict):
            continue
        for key in ("program", "label", "name"):
            val = _norm_text(p.get(key))
            if val and val.upper() != "OTHER":
                out.append(val)
                break
        for term_obj in _safe_list(p.get("representative_terms"))[:2]:
            if isinstance(term_obj, dict):
                term = _norm_text(term_obj.get("term"))
                if term:
                    out.append(term)
    return _dedupe(out)[:max_items]


def _extract_top_terms(triage: Dict[str, Any], max_items: int = 10) -> List[str]:
    out: List[str] = []
    for row in _safe_list((triage or {}).get("rows"))[:max_items]:
        if not isinstance(row, dict):
            continue
        for key in ("term", "pathway", "name", "description"):
            val = _norm_text(row.get(key))
            if val:
                out.append(val)
                break
    return _dedupe(out)[:max_items]


def _extract_top_genes(programs: Dict[str, Any], triage: Dict[str, Any], max_items: int = 12) -> List[str]:
    genes: List[str] = []
    for p in _safe_list((programs or {}).get("programs"))[:5]:
        if isinstance(p, dict):
            genes.extend(_norm_text(g) for g in _safe_list(p.get("top_genes"))[:8])
    if len(genes) < max_items:
        for row in _safe_list((triage or {}).get("rows"))[:10]:
            if not isinstance(row, dict):
                continue
            for key in ("genes", "gene_symbols", "leading_edge", "overlap_genes"):
                val = row.get(key)
                if isinstance(val, list):
                    genes.extend(_norm_text(g) for g in val[:8])
                elif isinstance(val, str):
                    genes.extend([x.strip() for x in re.split(r"[,;\s]+", val) if x.strip()][:8])
    return _dedupe([g for g in genes if g])[:max_items]


def _dedupe(items: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for item in items:
        key = item.lower().strip()
        if key and key not in seen:
            out.append(item.strip())
            seen.add(key)
    return out


def build_playbook_query(
    *,
    phenotype: str,
    context: Dict[str, Any],
    triage: Dict[str, Any],
    programs: Dict[str, Any],
) -> str:
    """Build a compact semantic query for playbook retrieval."""
    context = context or {}
    assay = _norm_text(context.get("assay"))
    assay_norm = _norm_assay(assay)
    ctx_bits = [
        f"assay: {assay_norm or assay}" if (assay_norm or assay) else "",
        f"phenotype: {phenotype}" if phenotype else "",
        f"tissue: {_norm_text(context.get('tissue'))}" if _norm_text(context.get("tissue")) else "",
        f"cell type: {_norm_text(context.get('cell_type'))}" if _norm_text(context.get("cell_type")) else "",
        f"perturbation: {_norm_text(context.get('perturbation'))}" if _norm_text(context.get("perturbation")) else "",
        f"organism: {_norm_text(context.get('organism'))}" if _norm_text(context.get("organism")) else "",
    ]
    program_bits = _extract_program_names(programs)
    term_bits = _extract_top_terms(triage)
    gene_bits = _extract_top_genes(programs, triage)

    return "\n".join([
        "Retrieve interpretation rules, assay confounders, evidence weighting guidance, tissue/cell-type priors, and follow-up experiments relevant to this enrichment analysis.",
        "Experimental context: " + "; ".join(x for x in ctx_bits if x),
        "Top biological programs/terms: " + "; ".join(_dedupe(program_bits + term_bits)[:12]),
        "Top genes: " + ", ".join(gene_bits[:12]),
        "Prefer rules that help classify likely driver vs reactive vs artifact/confounded vs uncertain.",
    ]).strip()


def _local_playbook_fallback(query: str, assay: str, max_chars: int = 6000) -> str:
    """
    Small fallback when no vector store is configured.
    This does a local keyword-scored extract rather than pasting whole playbooks.
    """
    assay_key = _norm_assay(assay)
    files = GLOBAL_PLAYBOOKS + ASSAY_TO_PLAYBOOK.get(assay_key, [])
    query_terms = {
        t.lower()
        for t in re.findall(r"[A-Za-z0-9_/-]{4,}", query)
        if t.lower() not in {"assay", "phenotype", "context", "genes", "terms", "programs", "retrieve"}
    }

    scored_chunks: List[tuple[int, str]] = []
    for fn in files:
        p = PLAYBOOK_DIR / fn
        if not p.exists():
            continue
        text = p.read_text(encoding="utf-8", errors="ignore")
        # Split by markdown headings so fallback still returns focused sections.
        sections = re.split(r"(?m)(?=^#{1,3}\s+)", text)
        for section in sections:
            sec = section.strip()
            if not sec:
                continue
            low = sec.lower()
            score = sum(1 for t in query_terms if t in low)
            # Always keep first-principles / evidence sections in the candidate pool.
            if any(x in low for x in ["evidence", "confounder", "validation", "artifact", "causal", "reactive"]):
                score += 2
            if score > 0:
                scored_chunks.append((score, f"## Source: {fn}\n{sec}"))

    if not scored_chunks:
        return "No local playbook guidance retrieved. Use the required interpretation behavior and be cautious."

    scored_chunks.sort(key=lambda x: x[0], reverse=True)
    out: List[str] = []
    total = 0
    for _, chunk in scored_chunks[:8]:
        remaining = max_chars - total
        if remaining <= 0:
            break
        clipped = chunk[:remaining]
        out.append(clipped)
        total += len(clipped)
    return "\n\n---\n\n".join(out).strip()


def retrieve_playbook_context(
    *,
    phenotype: str,
    context: Dict[str, Any],
    triage: Dict[str, Any],
    programs: Dict[str, Any],
    vector_store_id: Optional[str] = None,
    model: str = "gpt-5-mini",
    max_chunks: int = 8,
    max_chars: int = 9000,
) -> Dict[str, Any]:
    """
    Retrieval-first playbook grounding.

    Returns a compact dict with:
      - mode: vector_store | local_fallback | none/error
      - query: semantic retrieval query used
      - context_text: compact guidance to place in the main reasoning prompt
    """
    vs_id = (
        vector_store_id
        or os.environ.get("OPENAI_VECTOR_STORE_ID")
        or os.environ.get("VECTOR_STORE_ID")
    )
    assay = _norm_text((context or {}).get("assay"))
    query = build_playbook_query(
        phenotype=phenotype,
        context=context,
        triage=triage,
        programs=programs,
    )

    if not vs_id:
        fallback = _local_playbook_fallback(query=query, assay=assay, max_chars=max_chars)
        return {
            "mode": "local_fallback",
            "vector_store_id": "",
            "query": query,
            "context_text": fallback,
            "error": "OPENAI_VECTOR_STORE_ID/VECTOR_STORE_ID not set; used local focused fallback.",
        }

    retrieval_prompt = f"""
Use the file_search tool to retrieve only the most relevant playbook guidance for this enrichment analysis.
Do not summarize unrelated playbook material.
Return concise bullets grouped under:
- Assay limitations / confounders
- Evidence weighting / driver-vs-reactive-vs-artifact rules
- Tissue or cell-type priors
- Follow-up experiments and controls

Keep it compact. Include source file names when available.

Retrieval query:
{query}
""".strip()

    try:
        resp = client.responses.create(
            model=model,
            input=[{"role": "user", "content": retrieval_prompt}],
            tools=[{
                "type": "file_search",
                "vector_store_ids": [vs_id],
                "max_num_results": max_chunks,
            }],
            text={"format": {"type": "text"}},
        )
        text = (getattr(resp, "output_text", None) or "").strip()
        if not text:
            raise RuntimeError("No output_text returned from playbook retrieval response.")
        return {
            "mode": "vector_store",
            "vector_store_id": vs_id,
            "query": query,
            "context_text": text[:max_chars],
            "error": "",
        }
    except Exception as e:
        fallback = _local_playbook_fallback(query=query, assay=assay, max_chars=max_chars)
        return {
            "mode": "local_fallback_after_error",
            "vector_store_id": vs_id,
            "query": query,
            "context_text": fallback,
            "error": str(e),
        }

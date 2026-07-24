from __future__ import annotations

import os
import time
from typing import Any, Dict, List, Optional
import xml.etree.ElementTree as ET

import requests


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _safe_list(x: Any) -> List[Any]:
    return x if isinstance(x, list) else []


def _norm_text(s: Any) -> str:
    if s is None:
        return ""
    return str(s).strip()


def _pick_top_program_terms(programs: Dict[str, Any], max_programs: int = 3) -> List[str]:
    out: List[str] = []
    prog_list = _safe_list(programs.get("programs"))

    for prog in prog_list[:max_programs]:
        if not isinstance(prog, dict):
            continue

        prog_name = _norm_text(prog.get("program"))
        if prog_name and prog_name != "OTHER":
            out.append(prog_name)

        rep_terms = _safe_list(prog.get("representative_terms"))
        for term_obj in rep_terms[:2]:
            if isinstance(term_obj, dict):
                term = _norm_text(term_obj.get("term"))
                if term:
                    out.append(term)

    seen = set()
    deduped = []
    for x in out:
        if x not in seen:
            deduped.append(x)
            seen.add(x)
    return deduped[:8]


def _pick_top_genes(programs: Dict[str, Any], max_genes: int = 8) -> List[str]:
    genes: List[str] = []
    prog_list = _safe_list(programs.get("programs"))

    for prog in prog_list[:3]:
        if not isinstance(prog, dict):
            continue
        for g in _safe_list(prog.get("top_genes"))[:5]:
            g = _norm_text(g)
            if g:
                genes.append(g)

    seen = set()
    deduped = []
    for g in genes:
        if g not in seen:
            deduped.append(g)
            seen.add(g)
    return deduped[:max_genes]


def _build_pubmed_query(
    *,
    phenotype: str,
    context: Dict[str, Any],
    programs: Dict[str, Any],
) -> str:
    assay = _norm_text(context.get("assay"))
    tissue = _norm_text(context.get("tissue"))
    cell_type = _norm_text(context.get("cell_type"))
    perturbation = _norm_text(context.get("perturbation"))
    organism = _norm_text(context.get("organism"))

    top_terms = _pick_top_program_terms(programs)
    top_genes = _pick_top_genes(programs)

    query_parts: List[str] = []

    if phenotype:
        query_parts.append(f'("{phenotype}"[Title/Abstract])')

    ctx_bits = []
    for x in [cell_type, tissue, perturbation, organism]:
        if x:
            ctx_bits.append(f'"{x}"[Title/Abstract]')
    if ctx_bits:
        query_parts.append("(" + " AND ".join(ctx_bits) + ")")

    if top_terms:
        term_clause = " OR ".join(f'"{t}"[Title/Abstract]' for t in top_terms[:3])
        query_parts.append(f"({term_clause})")

    if top_genes:
        gene_clause = " OR ".join(f"{g}[Title/Abstract]" for g in top_genes[:4])
        query_parts.append(f"({gene_clause})")

    if assay:
        query_parts.append(f'("{assay}"[Title/Abstract])')

    if not query_parts:
        return "biomedical literature"

    return " AND ".join(query_parts)


def _build_fallback_queries(
    *,
    phenotype: str,
    context: Dict[str, Any],
    programs: Dict[str, Any],
) -> List[str]:
    tissue = _norm_text(context.get("tissue"))
    cell_type = _norm_text(context.get("cell_type"))
    perturbation = _norm_text(context.get("perturbation"))
    organism = _norm_text(context.get("organism"))

    top_terms = _pick_top_program_terms(programs)
    top_genes = _pick_top_genes(programs)

    queries: List[str] = []

    queries.append(_build_pubmed_query(
        phenotype=phenotype,
        context=context,
        programs=programs,
    ))

    if phenotype and perturbation and (cell_type or tissue):
        queries.append(
            " ".join(x for x in [phenotype, perturbation, cell_type or tissue, organism] if x)
        )

    if perturbation and (cell_type or tissue):
        queries.append(
            " ".join(x for x in [perturbation, cell_type or tissue, organism] if x)
        )

    if perturbation and top_terms:
        queries.append(
            " ".join([perturbation] + top_terms[:3] + ([cell_type or tissue] if (cell_type or tissue) else []))
        )

    if perturbation and top_genes:
        queries.append(
            " ".join([perturbation] + top_genes[:4] + ([cell_type or tissue] if (cell_type or tissue) else []))
        )

    if perturbation:
        broad_ctx = cell_type or tissue or "airway"
        queries.append(f"{perturbation} {broad_ctx}")

        # Productization patch 2: domain-aware fallback queries for the common
        # Bioconductor airway/dexamethasone use case and related steroid-response
        # datasets. These reduce drift toward generic airway inflammation papers.
        ctx_low = " ".join([phenotype, tissue, cell_type, perturbation]).lower()
        looks_like_dex_airway = (
            any(x in ctx_low for x in ["dexamethasone", "glucocorticoid", "corticosteroid"])
            and "airway" in ctx_low
        )
        if looks_like_dex_airway:
            queries.extend([
                "dexamethasone airway smooth muscle RNA-seq glucocorticoid",
                "dexamethasone airway smooth muscle FKBP5 TSC22D3 DUSP1 SGK1",
                "glucocorticoid receptor airway smooth muscle transcriptome",
                "dexamethasone human airway smooth muscle cytokine function CRISPLD2",
                "dexamethasone TNF alpha human airway smooth muscle NF-kappaB",
            ])

    if perturbation:
        # Generic fallbacks, ordered from more specific to broader.
        if cell_type:
            queries.append(f"{perturbation} glucocorticoid {cell_type}")
        if tissue:
            queries.append(f"{perturbation} glucocorticoid {tissue}")
        queries.append(f"{perturbation} airway smooth muscle")
        queries.append(f"{perturbation} airway epithelium")
        queries.append(f"{perturbation} epithelial cells")
        queries.append(f"{perturbation}")

    seen = set()
    deduped = []
    for q in queries:
        q = q.strip()
        if q and q not in seen:
            deduped.append(q)
            seen.add(q)

    return deduped


def _ncbi_params(extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    params: Dict[str, Any] = {
        "tool": "enrichment_llm_app",
        "retmode": "json",
    }

    email = os.environ.get("NCBI_EMAIL")
    api_key = os.environ.get("NCBI_API_KEY")

    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    if extra:
        params.update(extra)

    return params


def _esearch(query: str, retmax: int = 5) -> List[str]:
    params = _ncbi_params({
        "db": "pubmed",
        "term": query,
        "retmax": retmax,
        "sort": "relevance",
    })

    r = requests.get(ESEARCH_URL, params=params, timeout=20)
    r.raise_for_status()
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])


def _esummary(pmids: List[str]) -> List[Dict[str, Any]]:
    if not pmids:
        return []

    params = _ncbi_params({
        "db": "pubmed",
        "id": ",".join(pmids),
    })

    r = requests.get(ESUMMARY_URL, params=params, timeout=20)
    r.raise_for_status()
    data = r.json()

    result = data.get("result", {})
    out: List[Dict[str, Any]] = []

    for pmid in pmids:
        item = result.get(pmid, {})
        if not isinstance(item, dict):
            continue

        out.append({
            "pmid": pmid,
            "title": _norm_text(item.get("title")),
            "pubdate": _norm_text(item.get("pubdate")),
            "source": _norm_text(item.get("source")),
            "authors": [a.get("name", "") for a in item.get("authors", []) if isinstance(a, dict)],
            "doi": _norm_text(item.get("elocationid")),
        })

    return out


def _efetch_abstracts(pmids: List[str]) -> Dict[str, str]:
    if not pmids:
        return {}

    params = _ncbi_params({
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
    })

    r = requests.get(EFETCH_URL, params=params, timeout=25)
    r.raise_for_status()

    abstracts: Dict[str, str] = {}
    root = ET.fromstring(r.text)

    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        if pmid_el is None or not (pmid_el.text or "").strip():
            continue

        pmid = pmid_el.text.strip()
        abstract_nodes = article.findall(".//Abstract/AbstractText")

        parts: List[str] = []
        for node in abstract_nodes:
            text = "".join(node.itertext()).strip()
            label = node.attrib.get("Label", "").strip()
            if text:
                if label:
                    parts.append(f"{label}: {text}")
                else:
                    parts.append(text)

        abstracts[pmid] = " ".join(parts).strip()

    return abstracts


def _with_retry_esummary(pmids: List[str], retries: int = 1, delay_s: float = 1.2) -> List[Dict[str, Any]]:
    last_exc: Exception | None = None
    for attempt in range(retries + 1):
        try:
            return _esummary(pmids)
        except requests.HTTPError as e:
            last_exc = e
            if "429" in str(e) and attempt < retries:
                time.sleep(delay_s)
                continue
            raise
    if last_exc:
        raise last_exc
    return []


def _with_retry_abstracts(pmids: List[str], retries: int = 1, delay_s: float = 1.2) -> Dict[str, str]:
    last_exc: Exception | None = None
    for attempt in range(retries + 1):
        try:
            return _efetch_abstracts(pmids)
        except requests.HTTPError as e:
            last_exc = e
            if "429" in str(e) and attempt < retries:
                time.sleep(delay_s)
                continue
            raise
    if last_exc:
        raise last_exc
    return {}



# -------------------------
# Literature relevance labeling
# -------------------------

def _tokenize_for_relevance(text: Any) -> List[str]:
    """Small, dependency-free tokenizer for rough PubMed relevance labeling."""
    s = _norm_text(text).lower()
    toks = []
    for t in s.replace("-", " ").replace("/", " ").split():
        t = "".join(ch for ch in t if ch.isalnum())
        if len(t) >= 4 and t not in {"with", "from", "that", "this", "into", "cell", "cells", "gene", "genes", "study", "paper", "response"}:
            toks.append(t)
    return toks


def _paper_relevance(
    paper: Dict[str, Any],
    *,
    phenotype: str,
    context: Dict[str, Any],
    top_terms: List[str],
    top_genes: List[str],
) -> Dict[str, Any]:
    """
    Label PubMed hits without filtering them out.

    This is intentionally permissive: weak hits stay visible, but the app tells
    the user whether they should be treated as direct support, general support,
    or only background retrieval.
    """
    text = " ".join([
        _norm_text(paper.get("title")),
        _norm_text(paper.get("source")),
        _norm_text(paper.get("abstract")),
    ]).lower()

    score = 0.0
    reasons: List[str] = []

    # Context terms: phenotype and user-provided experimental context.
    context_values = [phenotype]
    for key in ["perturbation", "cell_type", "tissue", "organism", "assay"]:
        val = _norm_text((context or {}).get(key))
        if val and val.lower() not in {"not specified", "unknown", "na", "n/a"}:
            context_values.append(val)

    context_hits = []
    for val in context_values:
        val_norm = val.lower().strip()
        if len(val_norm) >= 4 and val_norm in text:
            context_hits.append(val)
    if context_hits:
        score += min(5.0, 2.0 * len(context_hits))
        reasons.append("matches experimental context: " + ", ".join(context_hits[:3]))

    # Exact gene hits are strong, but capped so one gene list does not dominate.
    gene_hits = []
    for g in top_genes[:10]:
        g_norm = _norm_text(g)
        if g_norm and len(g_norm) >= 3 and g_norm.lower() in text:
            gene_hits.append(g_norm)
    if gene_hits:
        score += min(6.0, 2.0 * len(gene_hits))
        reasons.append("mentions top genes: " + ", ".join(gene_hits[:5]))

    # Program/term phrase hits and token overlap.
    term_phrase_hits = []
    term_token_hits = []
    text_tokens = set(_tokenize_for_relevance(text))
    for term in top_terms[:8]:
        term_norm = _norm_text(term).lower()
        if not term_norm:
            continue
        if len(term_norm) >= 6 and term_norm in text:
            term_phrase_hits.append(term)
        else:
            toks = set(_tokenize_for_relevance(term))
            overlap = sorted(toks & text_tokens)
            if len(overlap) >= 2:
                term_token_hits.append(term)

    if term_phrase_hits:
        score += min(6.0, 3.0 * len(term_phrase_hits))
        reasons.append("matches top enriched terms: " + ", ".join(term_phrase_hits[:3]))
    elif term_token_hits:
        score += min(3.0, 1.0 * len(term_token_hits))
        reasons.append("partially overlaps top enriched terms: " + ", ".join(term_token_hits[:3]))

    # Use labels as guidance, not hard evidence calls.
    has_context = bool(context_hits)
    has_biology = bool(gene_hits or term_phrase_hits or term_token_hits)
    if score >= 8 and has_context and has_biology:
        label = "direct_support_candidate"
        explanation = "Potentially relevant to both the experiment context and the enriched biology. Still not causal proof."
    elif score >= 4 and has_biology:
        label = "general_support"
        explanation = "Relevant to at least part of the enriched biology, but not enough to treat as direct evidence."
    else:
        label = "weak_background"
        explanation = "Broad background retrieval only; do not use as direct support for an interpretation claim."

    return {
        "relevance_score": round(float(score), 2),
        "relevance_label": label,
        "relevance_reason": "; ".join(reasons) if reasons else explanation,
        "literature_use": explanation,
    }


def _annotate_pubmed_relevance(
    out: Dict[str, Any],
    *,
    phenotype: str,
    context: Dict[str, Any],
    programs: Dict[str, Any],
) -> Dict[str, Any]:
    top_terms = _pick_top_program_terms(programs)
    top_genes = _pick_top_genes(programs)
    papers = out.get("papers", []) or []

    annotated = []
    counts = {"direct_support_candidate": 0, "general_support": 0, "weak_background": 0}
    for p in papers:
        if not isinstance(p, dict):
            continue
        rel = _paper_relevance(p, phenotype=phenotype, context=context, top_terms=top_terms, top_genes=top_genes)
        p = {**p, **rel}
        counts[p["relevance_label"]] = counts.get(p["relevance_label"], 0) + 1
        annotated.append(p)

    # Sparse context means the search can still be useful, but should be framed as background.
    context_values = [
        _norm_text((context or {}).get("tissue")),
        _norm_text((context or {}).get("cell_type")),
        _norm_text((context or {}).get("perturbation")),
        _norm_text((context or {}).get("timepoint")),
        _norm_text((context or {}).get("organism")),
    ]
    informative_context = [x for x in context_values if x and x.lower() not in {"not specified", "unknown", "na", "n/a"}]
    sparse_context = len(informative_context) < 2

    if not annotated:
        quality = "none"
        reason = "No PubMed hits were retrieved."
    elif counts.get("direct_support_candidate", 0) > 0 and not sparse_context:
        quality = "higher"
        reason = "At least one retrieved paper overlaps both experimental context and enriched genes/terms."
    elif counts.get("general_support", 0) > 0:
        quality = "moderate" if not sparse_context else "low_to_moderate"
        reason = "Retrieved papers overlap parts of the enriched biology, but context is limited or incomplete."
    else:
        quality = "low"
        reason = "Retrieved papers should be treated as broad background because overlap with the input context/enriched biology is weak."

    if sparse_context and annotated:
        reason += " Tissue, cell type, perturbation, organism, or timepoint context is sparse."

    out["papers"] = annotated
    out["retrieval_quality"] = quality
    out["retrieval_quality_reason"] = reason
    out["relevance_counts"] = counts
    out["literature_use_guidance"] = (
        "Keep these PubMed hits visible, but use paper-level relevance_label values. "
        "Only direct_support_candidate papers should be cited as direct support; general_support papers are context; "
        "weak_background papers are background only."
    )
    return out


def fetch_pubmed_context(
    *,
    phenotype: str,
    context: Dict[str, Any],
    triage: Dict[str, Any],
    programs: Dict[str, Any],
    max_papers: int = 5,
) -> Dict[str, Any]:
    query_candidates = _build_fallback_queries(
        phenotype=phenotype,
        context=context,
        programs=programs,
    )

    out: Dict[str, Any] = {
        "query": query_candidates[0] if query_candidates else "",
        "query_candidates": query_candidates,
        "query_used": "",
        "query_strategy": "",
        "papers": [],
        "top_terms_used": _pick_top_program_terms(programs),
        "top_genes_used": _pick_top_genes(programs),
        "source": "PubMed via NCBI E-utilities",
        "status": "not_run",
        "retrieval_quality": "not_assessed",
        "retrieval_quality_reason": "PubMed retrieval has not run yet.",
        "literature_use_guidance": "Use PubMed hits as background unless relevance labels indicate stronger alignment.",
        "relevance_counts": {},
    }

    try:
        pmids: List[str] = []
        used_query = ""
        used_strategy = ""

        for i, query in enumerate(query_candidates):
            try:
                pmids = _esearch(query, retmax=max_papers)
            except Exception as e:
                out.setdefault("search_errors", []).append({
                    "query": query,
                    "error": str(e),
                })
                pmids = []

            if pmids:
                used_query = query
                used_strategy = "strict" if i == 0 else f"fallback_{i}"
                break

        if not pmids:
            out["status"] = "no_hits"
            out["retrieval_quality"] = "none"
            out["retrieval_quality_reason"] = "No PubMed hits were retrieved for the query candidates."
            return out

        # brief polite pause before metadata fetch
        time.sleep(0.4)

        try:
            summaries = _with_retry_esummary(pmids, retries=1, delay_s=1.2)
        except Exception as e:
            out["status"] = "partial_error"
            out["query_used"] = used_query
            out["query_strategy"] = used_strategy
            out["error"] = f"esummary failed: {e}"

            # fallback: still expose PMIDs and PubMed URLs
            out["papers"] = [
                {
                    "pmid": pmid,
                    "title": "",
                    "pubdate": "",
                    "source": "",
                    "authors": [],
                    "doi": "",
                    "abstract": "",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                }
                for pmid in pmids
            ]
            return _annotate_pubmed_relevance(out, phenotype=phenotype, context=context, programs=programs)

        time.sleep(0.4)

        try:
            abstracts = _with_retry_abstracts(pmids, retries=1, delay_s=1.2)
            abstract_status = "ok"
        except Exception as e:
            abstracts = {}
            abstract_status = "failed"
            out["abstract_error"] = str(e)

        papers: List[Dict[str, Any]] = []
        for s in summaries:
            pmid = s["pmid"]
            papers.append({
                **s,
                "abstract": abstracts.get(pmid, ""),
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })

        out["papers"] = papers
        out["query_used"] = used_query
        out["query_strategy"] = used_strategy
        out["abstract_status"] = abstract_status
        out["status"] = "ok" if papers else "partial_error"
        return _annotate_pubmed_relevance(out, phenotype=phenotype, context=context, programs=programs)

    except Exception as e:
        out["status"] = "error"
        out["error"] = str(e)
        return out
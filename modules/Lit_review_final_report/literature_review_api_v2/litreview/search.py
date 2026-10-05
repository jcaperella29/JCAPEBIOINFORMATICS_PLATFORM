from __future__ import annotations

import httpx

from .biology_playbook import DEFAULT_PLAYBOOK
from .models import Paper, StudyFinding


TIMEOUT = 30.0


def build_queries(f: StudyFinding) -> list[str]:
    """Liberal retrieval by design; strictness belongs in evidence classification."""
    core = ([f.target] if f.target else []) + f.genes[:6]
    biology = f.pathways[:4] + f.phenotypes[:4]
    context = [x for x in [f.tissue_or_cell_type, f.disease_or_context] if x]

    # Add nearby biological consequences to retrieval, but do not let their presence
    # determine evidence class.  The classifier separately verifies causal distance.
    expanded = sorted(DEFAULT_PLAYBOOK.terms_reachable_from(core + biology, max_distance=2))[:6]

    queries: list[str] = []

    # Assay-aware searches
    for term in (core[:3] + biology[:2]):
        if term:
            queries.append(f'"{term}" AND "{f.assay_type}"')

    # Assay-free searches are intentional: they catch mechanistic/contextual work
    # that a strict retrieval query would miss.
    for term in core[:3]:
        if term:
            queries.append(f'"{term}"')
    for term in biology[:3]:
        if term:
            queries.append(f'"{term}"')
    for term in expanded[:4]:
        if term:
            queries.append(f'"{term}"')

    # Pair a target with biology/context when available, without requiring assay.
    if core and biology:
        queries.append(f'"{core[0]}" AND "{biology[0]}"')
    if core and context:
        queries.append(f'"{core[0]}" AND "{context[0]}"')
    if biology and context:
        queries.append(f'"{biology[0]}" AND "{context[0]}"')

    return list(dict.fromkeys(q for q in queries if q))


async def search_pubmed(query: str, limit: int = 10) -> list[Paper]:
    async with httpx.AsyncClient(timeout=TIMEOUT) as client:
        es = await client.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={
            "db": "pubmed", "term": query, "retmode": "json", "retmax": limit,
        })
        es.raise_for_status()
        ids = es.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        sm = await client.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params={
            "db": "pubmed", "id": ",".join(ids), "retmode": "json"
        })
        sm.raise_for_status()
        data = sm.json().get("result", {})
        papers = []
        for pmid in ids:
            item = data.get(pmid, {})
            papers.append(Paper(
                source="pubmed", source_id=pmid, pmid=pmid,
                title=item.get("title", ""), journal=item.get("fulljournalname"),
                year=_year(item.get("pubdate")),
                authors=[a.get("name", "") for a in item.get("authors", [])],
                publication_type=item.get("pubtype", []),
                url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", raw=item,
            ))
        return papers


async def search_europe_pmc(query: str, limit: int = 10) -> list[Paper]:
    async with httpx.AsyncClient(timeout=TIMEOUT) as client:
        r = await client.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search", params={
            "query": query, "format": "json", "pageSize": limit, "resultType": "core"
        })
        r.raise_for_status()
        out = []
        for item in r.json().get("resultList", {}).get("result", []):
            out.append(Paper(
                source="europe_pmc",
                source_id=str(item.get("id") or item.get("pmid") or item.get("doi") or item.get("title")),
                pmid=item.get("pmid"), doi=item.get("doi"), title=item.get("title", ""),
                abstract=item.get("abstractText"), journal=item.get("journalTitle"),
                year=_safe_int(item.get("pubYear")),
                authors=[x.strip() for x in (item.get("authorString") or "").split(",") if x.strip()],
                publication_type=item.get("pubTypeList", {}).get("pubType", []) if isinstance(item.get("pubTypeList"), dict) else [],
                url=f"https://europepmc.org/article/MED/{item.get('pmid')}" if item.get("pmid") else None,
                raw=item,
            ))
        return out


async def search_openalex(query: str, limit: int = 10) -> list[Paper]:
    async with httpx.AsyncClient(timeout=TIMEOUT) as client:
        r = await client.get("https://api.openalex.org/works", params={"search": query, "per-page": limit})
        r.raise_for_status()
        out = []
        for item in r.json().get("results", []):
            doi = item.get("doi")
            out.append(Paper(
                source="openalex",
                source_id=item.get("id", ""),
                title=item.get("display_name", ""),
                year=item.get("publication_year"),
                journal=((item.get("primary_location") or {}).get("source") or {}).get("display_name"),
                doi=doi.replace("https://doi.org/", "") if doi else None,
                authors=[a.get("author", {}).get("display_name", "") for a in item.get("authorships", [])],
                publication_type=[item.get("type")] if item.get("type") else [],
                url=item.get("id"), raw=item,
            ))
        return out


def _safe_int(x):
    try:
        return int(x)
    except Exception:
        return None


def _year(x):
    if not x:
        return None
    for token in str(x).split():
        if token[:4].isdigit():
            return int(token[:4])
    return None

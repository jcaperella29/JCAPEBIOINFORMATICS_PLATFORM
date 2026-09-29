from __future__ import annotations

import asyncio
from collections import Counter

from .classifier import classify_evidence
from .models import EvidenceAssessment, ReviewRequest, ReviewResponse, Paper
from .search import build_queries, search_pubmed, search_europe_pmc, search_openalex


def _dedup_key(p: Paper) -> str:
    if p.doi:
        return "doi:" + p.doi.lower().strip()
    if p.pmid:
        return "pmid:" + p.pmid.strip()
    return "title:" + " ".join(p.title.lower().split())


async def run_review(req: ReviewRequest) -> ReviewResponse:
    assessments: list[EvidenceAssessment] = []

    for finding in req.handoff.findings:
        queries = build_queries(finding)
        papers: list[Paper] = []
        for q in queries:
            batches = await asyncio.gather(
                search_pubmed(q, req.max_papers_per_finding),
                search_europe_pmc(q, req.max_papers_per_finding),
                search_openalex(q, req.max_papers_per_finding),
                return_exceptions=True,
            )
            for batch in batches:
                if isinstance(batch, list):
                    papers.extend(batch)

        deduped: dict[str, Paper] = {}
        for p in papers:
            deduped.setdefault(_dedup_key(p), p)

        these = [classify_evidence(finding, p) for p in deduped.values()]
        these.sort(key=lambda x: x.relevance_score, reverse=True)
        assessments.extend(these[: req.max_papers_per_finding])

    counts = Counter(a.level.value for a in assessments)
    return ReviewResponse(
        run_id=req.handoff.run_id,
        assessments=assessments,
        summary_counts=dict(counts),
    )

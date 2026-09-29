from __future__ import annotations

from enum import Enum
from typing import Any, Literal
from pydantic import BaseModel, Field


class EvidenceLevel(str, Enum):
    DIRECT = "direct"
    MECHANISTIC = "mechanistic_support"
    CONTEXTUAL = "contextual_support"
    BACKGROUND = "background"
    EXCLUDE = "exclude"


class StudyFinding(BaseModel):
    finding_id: str
    assay_type: str = Field(description="Assay used by the current experiment, e.g. perturb-seq, scRNA-seq, bulk RNA-seq")
    target: str | None = None
    genes: list[str] = []
    pathways: list[str] = []
    phenotypes: list[str] = []
    organism: str | None = None
    tissue_or_cell_type: str | None = None
    disease_or_context: str | None = None
    summary: str


class PipelineHandoff(BaseModel):
    run_id: str
    analysis_type: str
    findings: list[StudyFinding]
    quantitative_results: dict[str, Any] = {}
    network_results: dict[str, Any] = {}
    llm_triage: dict[str, Any] = {}
    metadata: dict[str, Any] = {}


class Paper(BaseModel):
    source: Literal["pubmed", "europe_pmc", "openalex"]
    source_id: str
    title: str
    abstract: str | None = None
    year: int | None = None
    journal: str | None = None
    doi: str | None = None
    pmid: str | None = None
    url: str | None = None
    authors: list[str] = []
    publication_type: list[str] = []
    raw: dict[str, Any] = {}


class EvidenceAssessment(BaseModel):
    paper: Paper
    finding_id: str
    level: EvidenceLevel
    assay_match: bool = False
    same_target_or_question: bool = False
    directly_impacted_pathway_or_phenotype: bool = False
    playbook_match: bool = False
    playbook_distance: int | None = None
    playbook_path: list[str] = []
    playbook_relations: list[str] = []
    playbook_explanation: str | None = None
    organism_match: bool = False
    context_match: bool = False
    review_article: bool = False
    reasons: list[str] = []
    caution: str | None = None
    relevance_score: float = 0.0


class ReviewRequest(BaseModel):
    handoff: PipelineHandoff
    max_papers_per_finding: int = 20


class ReviewResponse(BaseModel):
    run_id: str
    assessments: list[EvidenceAssessment]
    summary_counts: dict[str, int]

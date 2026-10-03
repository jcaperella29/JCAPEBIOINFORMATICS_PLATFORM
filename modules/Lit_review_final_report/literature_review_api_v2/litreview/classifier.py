from __future__ import annotations

import re
from .biology_playbook import DEFAULT_PLAYBOOK
from .models import EvidenceAssessment, EvidenceLevel, Paper, StudyFinding


ASSAY_ALIASES = {
    "perturb-seq": {"perturb-seq", "perturb seq", "crispr single-cell", "crispr single cell", "single-cell crispr", "single cell crispr", "pooled crispr single-cell"},
    "scrna-seq": {"scrna-seq", "single-cell rna-seq", "single cell rna seq", "single-cell transcriptomics", "single cell transcriptomics"},
    "bulk rna-seq": {"bulk rna-seq", "rna-seq", "rna sequencing", "transcriptome sequencing"},
}

REVIEW_TERMS = {"review", "systematic review", "meta-analysis", "meta analysis"}


def _norm(text: str | None) -> str:
    return re.sub(r"\s+", " ", (text or "").lower()).strip()


def _paper_text(paper: Paper) -> str:
    return _norm(" ".join([paper.title or "", paper.abstract or "", " ".join(paper.publication_type)]))


def _contains_term(text: str, value: str) -> bool:
    """Match a biological term without allowing gene-symbol substring accidents.

    Example: MET must not match "metabolic" and CAT must not match "catalysis".
    Multi-word pathway/phenotype phrases still match normally.
    """
    term = _norm(value)
    if not term:
        return False
    pattern = rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])"
    return re.search(pattern, text) is not None


def _contains_any(text: str, values: list[str] | set[str]) -> bool:
    return any(_contains_term(text, v) for v in values)


def assay_matches(current_assay: str, text: str) -> bool:
    ca = _norm(current_assay)
    aliases = set()
    for key, vals in ASSAY_ALIASES.items():
        if ca == key or ca in vals:
            aliases |= vals
            aliases.add(key)
    if not aliases:
        aliases = {ca}
    return any(a in text for a in aliases if a)


def classify_evidence(finding: StudyFinding, paper: Paper) -> EvidenceAssessment:
    text = _paper_text(paper)
    assay_match = assay_matches(finding.assay_type, text)

    target_terms = [x for x in [finding.target] if x] + finding.genes
    target_match = _contains_any(text, target_terms)
    pathway_match = _contains_any(text, finding.pathways)
    phenotype_match = _contains_any(text, finding.phenotypes)
    explicit_direct_biology = pathway_match or phenotype_match

    # Biological consequence playbook: a short, interpretable relationship chain can
    # establish that the paper studies an immediately impacted pathway/phenotype even
    # if it does not repeat the exact finding term.
    biology_seeds = target_terms + finding.pathways + finding.phenotypes
    direct_path = DEFAULT_PLAYBOOK.find_path_to_text(
        biology_seeds,
        text,
        max_distance=2,
        allowed_relations=DEFAULT_PLAYBOOK.direct_relations,
    )
    mechanistic_path = direct_path or DEFAULT_PLAYBOOK.find_path_to_text(
        biology_seeds,
        text,
        max_distance=3,
        allowed_relations=DEFAULT_PLAYBOOK.mechanistic_relations,
    )
    playbook_direct = direct_path is not None
    playbook_mechanistic = mechanistic_path is not None
    direct_biology = explicit_direct_biology or playbook_direct

    organism_match = bool(finding.organism and _norm(finding.organism) in text)
    context_terms = [x for x in [finding.tissue_or_cell_type, finding.disease_or_context] if x]
    context_match = _contains_any(text, context_terms)

    pubtypes = {_norm(x) for x in paper.publication_type}
    review_article = bool(pubtypes & REVIEW_TERMS) or any(term in _norm(paper.title) for term in REVIEW_TERMS)

    reasons: list[str] = []
    cautions: list[str] = []

    # DIRECT is intentionally strict:
    # 1) same/equivalent assay, AND
    # 2) same target/question OR same named pathway/phenotype OR a <=2-edge
    #    consequence path using only direct causal/mechanistic relations,
    # 3) primary experimental literature.
    if assay_match and (target_match or direct_biology) and not review_article:
        level = EvidenceLevel.DIRECT
        reasons.append("Same or equivalent assay detected.")
        if target_match:
            reasons.append("Paper examines the same target/question or named gene(s).")
        if explicit_direct_biology:
            reasons.append("Paper explicitly examines a pathway or phenotype named by the current finding.")
        if playbook_direct and direct_path:
            reasons.append(f"Biological playbook links the finding to a measured paper endpoint within {direct_path.distance} step(s).")
        if not context_match and context_terms:
            cautions.append("Direct by assay/biology rule, but experimental context differs or is not explicit.")

    # Mechanistic evidence may use a slightly longer (<=3-edge) consequence chain,
    # but it cannot be upgraded to direct without assay equivalence.
    elif not review_article and ((target_match and (explicit_direct_biology or playbook_mechanistic)) or playbook_mechanistic):
        level = EvidenceLevel.MECHANISTIC
        if target_match:
            reasons.append("Same target/gene is present.")
        if explicit_direct_biology:
            reasons.append("An implicated pathway/phenotype is explicitly present.")
        if playbook_mechanistic and mechanistic_path:
            reasons.append(f"Biological playbook provides a mechanistic chain of {mechanistic_path.distance} step(s).")
        if not assay_match:
            cautions.append("Assay equivalence was not established, so this cannot count as direct evidence.")

    elif not review_article and (target_match or explicit_direct_biology):
        level = EvidenceLevel.CONTEXTUAL
        reasons.append("Biologically related evidence was found, but it does not satisfy the direct-evidence assay/causal-distance rule.")

    elif review_article and (target_match or explicit_direct_biology or playbook_mechanistic):
        level = EvidenceLevel.BACKGROUND
        reasons.append("Review-level/background literature; cannot count as direct experimental evidence.")

    else:
        level = EvidenceLevel.EXCLUDE
        reasons.append("Insufficient match to the current finding for evidence use.")

    path = direct_path or mechanistic_path

    # Ranking is separate from evidence class.  Retrieval remains liberal; the score
    # simply orders the papers inside the evidence framework.
    score = 0.0
    score += 4.0 if assay_match else 0.0
    score += 3.0 if target_match else 0.0
    score += 2.5 if explicit_direct_biology else 0.0
    if path:
        score += {1: 2.5, 2: 1.75, 3: 1.0}.get(path.distance, 0.0)
    score += 1.0 if context_match else 0.0
    score += 0.5 if organism_match else 0.0
    score -= 2.0 if review_article else 0.0

    return EvidenceAssessment(
        paper=paper,
        finding_id=finding.finding_id,
        level=level,
        assay_match=assay_match,
        same_target_or_question=target_match,
        directly_impacted_pathway_or_phenotype=direct_biology,
        playbook_match=path is not None,
        playbook_distance=path.distance if path else None,
        playbook_path=list(path.nodes) if path else [],
        playbook_relations=list(path.relations) if path else [],
        playbook_explanation=path.explanation if path else None,
        organism_match=organism_match,
        context_match=context_match,
        review_article=review_article,
        reasons=reasons,
        caution=" ".join(cautions) if cautions else None,
        relevance_score=score,
    )

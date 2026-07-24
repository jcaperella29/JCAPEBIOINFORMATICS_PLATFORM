# biofit.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import re
import math

import numpy as np


_WORD_RE = re.compile(r"[A-Za-z0-9\-\+]+")
_WS_RE = re.compile(r"\s+")


def _norm_text(s: str) -> str:
    s = (s or "").lower()
    s = s.replace("_", " ")
    s = _WS_RE.sub(" ", s).strip()
    return s


def _tokens(s: str) -> List[str]:
    s = _norm_text(s)
    return [t.lower() for t in _WORD_RE.findall(s)]


def _has_any(text: str, patterns: List[str]) -> bool:
    t = _norm_text(text)
    return any(p in t for p in patterns)


def _count_any(text: str, patterns: List[str]) -> int:
    t = _norm_text(text)
    return sum(1 for p in patterns if p in t)


def _clip01(x: float) -> float:
    return float(np.clip(x, 0.0, 1.0))


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


# -------------------------
# Small default ontology
# -------------------------

@dataclass
class Program:
    name: str
    # keywords that identify the program in pathway/term names
    term_keywords: List[str]
    # phenotype keywords that imply this program should matter
    phenotype_keywords: List[str]
    # gene family hints (regex) that indicate this program is actually represented in driver genes
    gene_regex: List[str]
    # common confounder note (optional)
    confounder: Optional[str] = None


DEFAULT_PROGRAMS: List[Program] = [
    Program(
        name="ECM_FIBROSIS",
        term_keywords=["extracellular matrix", "collagen", "ecm", "focal adhesion", "integrin", "tgf", "wound healing", "myofibroblast"],
        phenotype_keywords=["fibrosis", "scar", "ecm", "collagen", "stiffness", "matrix", "myofibroblast", "tgf"],
        gene_regex=[r"^COL\d+", r"^FN1$", r"^POSTN$", r"^SPARC$", r"^TGFB\d", r"^SMAD\d", r"^ITGA", r"^ITGB", r"^MMP\d", r"^TIMP\d"],
        confounder="Often core to fibroblast activation; can also reflect cell-state shifts/composition.",
    ),
    Program(
        name="INNATE_IFN_ANTIVIRAL",
        term_keywords=["interferon", "ifn", "antiviral", "defense response to virus", "pattern recognition", "rig-i", "tlr", "jak-stat"],
        phenotype_keywords=["viral", "ifn", "interferon", "innate immune", "antiviral", "inflammation"],
        gene_regex=[r"^IFIT", r"^ISG\d+", r"^OAS\d", r"^MX\d", r"^STAT\d", r"^IRF\d", r"^DDX\d+", r"^TLR\d", r"^CXCL\d+"],
        confounder="Strong IFN programs can dominate many datasets; check if it's expected vs contamination/stimulation artifact.",
    ),
    Program(
        name="ADAPTIVE_T_CELL",
        term_keywords=["t cell", "t-cell", "cd3", "antigen presentation", "mhc", "cytotoxic", "pd-1", "checkpoint", "exhaustion"],
        phenotype_keywords=["t cell", "exhaustion", "checkpoint", "pd-1", "cytotoxic", "antigen", "mhc"],
        gene_regex=[r"^CD3", r"^TRAC$", r"^TRBC", r"^GZMB$", r"^PRF1$", r"^LAG3$", r"^PDCD1$", r"^CTLA4$", r"^HLA-", r"^B2M$"],
    ),
    Program(
        name="DDR_P53_APOPTOSIS",
        term_keywords=["dna damage", "p53", "apoptosis", "checkpoint", "g2m", "g1s", "repair", "homologous recombination", "nhej"],
        phenotype_keywords=["dna damage", "genotoxic", "apoptosis", "senescence", "p53", "repair"],
        gene_regex=[r"^TP53$", r"^CDKN1A$", r"^GADD45", r"^BAX$", r"^BBC3$", r"^ATM$", r"^ATR$", r"^BRCA\d", r"^RAD\d", r"^CHEK\d"],
        confounder="Cell cycle & DDR can be downstream of stress/proliferation changes—interpret causality carefully.",
    ),
    Program(
        name="UPR_PROTEOSTASIS",
        term_keywords=["unfolded protein", "upr", "er stress", "endoplasmic reticulum", "protein folding", "proteasome", "heat shock"],
        phenotype_keywords=["er stress", "upr", "proteostasis", "misfold", "heat shock", "proteasome"],
        gene_regex=[r"^HSPA", r"^HSPB", r"^DNAJ", r"^XBP1$", r"^ATF4$", r"^DDIT3$", r"^HERPUD", r"^PSM", r"^UBB$"],
        confounder="Common ‘stress signature’—can be real biology or handling/viability artifact.",
    ),
    Program(
        name="MITO_OXPHOS",
        term_keywords=["mitochond", "oxidative phosphorylation", "oxphos", "respiratory chain", "tca", "electron transport"],
        phenotype_keywords=["mitochond", "oxphos", "energy", "respiration", "tca"],
        gene_regex=[r"^MT-", r"^NDUF", r"^COX\d", r"^ATP5", r"^SDH", r"^UQCR", r"^CS$", r"^IDH\d", r"^PDHA"],
        confounder="Mito terms swing with composition and overall expression quality; check for global shifts.",
    ),
]

# programs we often treat as confounders unless phenotype explicitly wants them
DEFAULT_ARTIFACT_BUCKETS = {
    "RIBOSOME_TRANSLATION": {
        "term_keywords": ["ribosome", "translation", "rRNA", "ribosomal", "srp"],
        "gene_regex": [r"^RPL", r"^RPS", r"^EEF", r"^EIF"],
        "note": "Often reflects growth rate, stress, or library composition; treat as downstream unless phenotype supports it.",
    },
    "CELL_CYCLE_PROLIFERATION": {
        "term_keywords": ["cell cycle", "mitotic", "g2m", "g1s", "dna replication", "m phase", "e2f"],
        "gene_regex": [r"^MKI67$", r"^TOP2A$", r"^CCN", r"^CDK", r"^PCNA$", r"^MCM\d", r"^UBE2C$"],
        "note": "Extremely common confounder; real if phenotype is proliferation, growth, tumor, regeneration.",
    },
}


@dataclass
class BioFitConfig:
    programs: List[Program] = None
    artifact_buckets: Dict = None

    # weights (sum doesn't need to be 1; we'll rescale)
    w_term_program: float = 0.45
    w_pheno_program: float = 0.25
    w_gene_support: float = 0.25
    w_system_penalty: float = 0.10

    # penalties
    artifact_penalty: float = 0.25      # subtract fraction if looks like confounder
    tiny_overlap_penalty: float = 0.20  # subtract fraction if overlap is tiny
    min_genes_for_confidence: int = 3

    def __post_init__(self):
        if self.programs is None:
            self.programs = DEFAULT_PROGRAMS
        if self.artifact_buckets is None:
            self.artifact_buckets = DEFAULT_ARTIFACT_BUCKETS


def _gene_family_hits(genes: List[str], regexes: List[str]) -> int:
    if not genes or not regexes:
        return 0
    hits = 0
    for g in genes:
        for rx in regexes:
            if re.search(rx, g, flags=re.IGNORECASE):
                hits += 1
                break
    return hits


def _program_match_strength(term: str, program: Program) -> float:
    # Term keyword match: count then squash
    c = _count_any(term, program.term_keywords)
    return _clip01(c / 2.0)  # 0, 0.5, 1.0...


def _phenotype_support_strength(phenotype: str, program: Program) -> float:
    c = _count_any(phenotype, program.phenotype_keywords)
    return _clip01(c / 2.0)


def _system_penalty(term: str, context: Dict) -> float:
    """
    Lightweight plausibility checks.
    You can expand this later (cell-type marker sets, tissue whitelist, etc.).
    """
    tissue = _norm_text(context.get("tissue", ""))
    celltype = _norm_text(context.get("cell_type", ""))
    assay = _norm_text(context.get("assay", ""))

    t = _norm_text(term)

    penalty = 0.0

    # neuron/synapse in non-neural contexts
    if any(k in t for k in ["synapse", "neuron", "axonal", "dopamine", "glutamate"]):
        if not any(k in (tissue + " " + celltype) for k in ["brain", "neuron", "astro", "microglia"]):
            penalty += 0.6

    # muscle contraction in non-muscle
    if any(k in t for k in ["muscle contraction", "sarcomere", "myofibril"]):
        if not any(k in (tissue + " " + celltype) for k in ["muscle", "cardio", "myocyte"]):
            penalty += 0.5

    # immunoglobulin / BCR in non-immune
    if any(k in t for k in ["b cell", "b-cell", "immunoglobulin", "bcr signaling"]):
        if not any(k in (tissue + " " + celltype) for k in ["b cell", "lymph", "spleen", "pbmc", "immune"]):
            penalty += 0.4

    # ATAC/ChIP: translation/ribosome often less interpretable as “primary”
    if any(k in assay for k in ["atac", "chip"]):
        if any(k in t for k in ["ribosome", "translation"]):
            penalty += 0.2

    return float(np.clip(penalty, 0.0, 1.0))


def _artifact_likelihood(term: str, genes: List[str], cfg: BioFitConfig) -> Tuple[float, List[str]]:
    """
    Returns artifact_likelihood in [0,1] and reasons.
    """
    t = _norm_text(term)
    reasons = []
    likelihood = 0.0

    for name, bucket in cfg.artifact_buckets.items():
        if _has_any(t, bucket["term_keywords"]):
            hits = _gene_family_hits(genes, bucket["gene_regex"])
            # if most genes are from that family, likely confounder bucket
            frac = hits / max(len(genes), 1)
            # be more confident when overlap is mostly that family
            if frac >= 0.6 and len(genes) >= cfg.min_genes_for_confidence:
                likelihood = max(likelihood, 0.85)
                reasons.append(f"{name.lower()}_dominant_genes")
            else:
                likelihood = max(likelihood, 0.55)
                reasons.append(f"{name.lower()}_term_match")

    return float(np.clip(likelihood, 0.0, 1.0)), reasons


def biofit_score(
    term: str,
    genes: List[str],
    overlap_k: int,
    phenotype: str,
    context: Optional[Dict] = None,
    cfg: Optional[BioFitConfig] = None,
) -> Dict:
    """
    Deterministic biological plausibility score.
    Returns dict with score and diagnostic components for UI/debugging.
    """
    cfg = cfg or BioFitConfig()
    context = context or {}

    phenotype_n = _norm_text(phenotype)
    term_n = _norm_text(term)

    # program alignment
    best_program = None
    best_prog = 0.0
    best_pheno = 0.0
    best_gene = 0.0

    for prog in cfg.programs:
        pm = _program_match_strength(term_n, prog)
        ph = _phenotype_support_strength(phenotype_n, prog)
        gh = 0.0
        if genes:
            hits = _gene_family_hits(genes, prog.gene_regex)
            gh = hits / max(len(genes), 1)
            gh = _clip01(gh / 0.35)  # saturate once ~35% of overlap matches family

        # a program "wins" when term matches + (phenotype OR genes support)
        combined = 0.55 * pm + 0.25 * ph + 0.20 * gh
        if combined > best_prog:
            best_prog = combined
            best_program = prog.name
            best_pheno = ph
            best_gene = gh

    # system plausibility penalty
    sys_pen = _system_penalty(term_n, context)

        # artifact likelihood
    art_like, art_reasons = _artifact_likelihood(term_n, genes, cfg)

    # base score from components (rescale weights)
    wsum = (cfg.w_term_program + cfg.w_pheno_program + cfg.w_gene_support + cfg.w_system_penalty)
    if wsum <= 0:
        wsum = 1.0

    term_component = best_prog  # already 0..1
    pheno_component = best_pheno  # 0..1
    gene_component = best_gene  # 0..1

    # system penalty: higher penalty should reduce score
    system_component = 1.0 - sys_pen  # 1 is good, 0 is bad

    base01 = (
        cfg.w_term_program * term_component +
        cfg.w_pheno_program * pheno_component +
        cfg.w_gene_support * gene_component +
        cfg.w_system_penalty * system_component
    ) / wsum

    base01 = _clip01(base01)

    # penalties
    penalties = 0.0
    flags: List[str] = []

    if overlap_k <= 2:
        penalties += cfg.tiny_overlap_penalty
        flags.append("tiny_overlap_biofit")

    if art_like >= 0.8:
        penalties += cfg.artifact_penalty
        flags.append("likely_artifact")
    elif art_like >= 0.55:
        penalties += cfg.artifact_penalty * 0.6
        flags.append("possible_artifact")

    final01 = _clip01(base01 * (1.0 - penalties))

    # score 0..100
    score = float(np.clip(final01 * 100.0, 0.0, 100.0))

    out = {
        "biofit_score": score,
        "best_program": best_program,
        "components": {
            "term_program": float(term_component),
            "phenotype_program": float(pheno_component),
            "gene_support": float(gene_component),
            "system_plausibility": float(system_component),
            "artifact_likelihood": float(art_like),
            "penalties_total": float(penalties),
            "base01": float(base01),
            "final01": float(final01),
        },
        "artifact_reasons": art_reasons,
        "flags": flags,
    }

    return out

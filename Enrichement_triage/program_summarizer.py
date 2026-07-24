# program_summarizer.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import re
from collections import Counter, defaultdict

import numpy as np

# Reuse Program / BioFitConfig ontology from biofit.py
# from biofit import BioFitConfig, Program, DEFAULT_PROGRAMS, DEFAULT_ARTIFACT_BUCKETS
# If you prefer not to import, you can pass programs/artifacts in config.

_WORD_RE = re.compile(r"[A-Za-z0-9\-\+]+")


def _norm(s: str) -> str:
    return (s or "").lower().replace("_", " ").strip()


def _count_term_hits(term: str, keywords: List[str]) -> int:
    t = _norm(term)
    return sum(1 for k in keywords if k in t)


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


@dataclass
class ProgramRule:
    """
    A lighter rule wrapper that mirrors Program but includes weights.
    """
    name: str
    term_keywords: List[str]
    gene_regex: List[str]
    # Some programs are "confounders unless phenotype matches"
    confounder_like: bool = False


@dataclass
class ProgramSummaryConfig:
    # program rules
    programs: List[ProgramRule]

    # optional phenotype text to decide whether to downweight confounders
    phenotype: str = ""

    # label used when no good match
    other_label: str = "OTHER"

    # assignment thresholds
    min_term_hit: int = 1
    min_total_support: float = 0.30  # combined term+gene support threshold

    # scoring weights
    w_term: float = 0.55
    w_gene: float = 0.45

    # when combining member terms into a program score:
    # use top_k mean + max so one strong term can lift the program but not dominate
    top_k_for_mean: int = 5
    w_prog_max: float = 0.60
    w_prog_mean: float = 0.40

    # confounder penalty if phenotype doesn't suggest it
    confounder_penalty: float = 0.20


def default_program_rules() -> List[ProgramRule]:
    # You can expand this list over time; this is a strong starting set.
    return [
        ProgramRule(
            name="GLUCOCORTICOID_GR_RESPONSE",
            term_keywords=[
                "glucocorticoid", "dexamethasone", "corticosteroid",
                "steroid hormone", "glucocorticoid receptor", "nr3c1",
                "response to steroid", "response to corticosteroid",
                "airway smooth muscle response to corticosteroid",
            ],
            gene_regex=[
                r"^FKBP5$", r"^TSC22D3$", r"^DUSP1$", r"^SGK1$",
                r"^PER1$", r"^KLF9$", r"^ZBTB16$", r"^NR3C1$",
                r"^CRISPLD2$", r"^BCL6$",
            ],
        ),
        ProgramRule(
            name="MAPK_STRESS_KINASE_ATTENUATION",
            term_keywords=[
                "mapk", "stress kinase", "p38", "jnk", "map kinase",
                "phosphatase", "kinase attenuation", "stress-kinase",
            ],
            gene_regex=[
                r"^DUSP\d+$", r"^MAPK", r"^FOS$", r"^JUN$", r"^KLF9$",
            ],
        ),
        ProgramRule(
            name="ECM_FIBROSIS",
            term_keywords=["extracellular matrix", "collagen", "ecm", "integrin", "focal adhesion", "tgf", "wound healing", "myofibroblast"],
            gene_regex=[r"^COL\d+", r"^FN1$", r"^POSTN$", r"^SPARC$", r"^TGFB\d", r"^SMAD\d", r"^ITGA", r"^ITGB", r"^MMP\d", r"^TIMP\d"],
        ),
        ProgramRule(
            name="INNATE_IFN_ANTIVIRAL",
            term_keywords=["interferon", "ifn", "antiviral", "defense response to virus", "jak-stat", "tlr", "rig-i", "pattern recognition"],
            gene_regex=[r"^IFIT", r"^ISG\d+", r"^OAS\d", r"^MX\d", r"^STAT\d", r"^IRF\d", r"^DDX\d+", r"^TLR\d", r"^CXCL\d+"],
        ),
        ProgramRule(
            name="ADAPTIVE_T_CELL",
            term_keywords=["t cell", "t-cell", "cd3", "antigen presentation", "mhc", "cytotoxic", "pd-1", "checkpoint", "exhaustion"],
            gene_regex=[r"^CD3", r"^TRAC$", r"^TRBC", r"^GZMB$", r"^PRF1$", r"^LAG3$", r"^PDCD1$", r"^CTLA4$", r"^HLA-", r"^B2M$"],
        ),
        ProgramRule(
            name="DDR_P53_APOPTOSIS",
            term_keywords=["dna damage", "p53", "apoptosis", "checkpoint", "repair", "homologous recombination", "nhej", "senescence"],
            gene_regex=[r"^TP53$", r"^CDKN1A$", r"^GADD45", r"^BAX$", r"^BBC3$", r"^ATM$", r"^ATR$", r"^BRCA\d", r"^RAD\d", r"^CHEK\d"],
        ),
        ProgramRule(
            name="UPR_PROTEOSTASIS",
            term_keywords=["unfolded protein", "upr", "er stress", "endoplasmic reticulum", "protein folding", "proteasome", "heat shock"],
            gene_regex=[r"^HSPA", r"^HSPB", r"^DNAJ", r"^XBP1$", r"^ATF4$", r"^DDIT3$", r"^HERPUD", r"^PSM", r"^UBB$"],
        ),
        ProgramRule(
            name="MITO_OXPHOS",
            term_keywords=["mitochond", "oxidative phosphorylation", "oxphos", "respiratory chain", "tca", "electron transport"],
            gene_regex=[r"^MT-", r"^NDUF", r"^COX\d", r"^ATP5", r"^SDH", r"^UQCR", r"^CS$", r"^IDH\d", r"^PDHA"],
        ),
        ProgramRule(
            name="INFLAMMATION_NFKB",
            term_keywords=["nf-kb", "nfkb", "tnf", "il-1", "inflammatory", "cytokine", "chemokine", "toll-like receptor"],
            gene_regex=[r"^NFKB", r"^RELA$", r"^TNF$", r"^IL1", r"^CXCL", r"^CCL", r"^ICAM1$", r"^SELE$", r"^PTGS2$"],
        ),
        ProgramRule(
            name="ANGIO_HYPOXIA",
            term_keywords=["hypoxia", "hif", "angiogenesis", "vegf", "vascular", "endothelial"],
            gene_regex=[r"^HIF1A$", r"^VEGFA$", r"^KDR$", r"^FLT1$", r"^ANGPT", r"^EGLN", r"^CA9$"],
        ),
        # Confounder-like buckets:
        ProgramRule(
            name="CELL_CYCLE_PROLIFERATION",
            term_keywords=["cell cycle", "mitotic", "g2m", "g1s", "dna replication", "m phase", "e2f"],
            gene_regex=[r"^MKI67$", r"^TOP2A$", r"^CCN", r"^CDK", r"^PCNA$", r"^MCM\d", r"^UBE2C$"],
            confounder_like=True,
        ),
        ProgramRule(
            name="RIBOSOME_TRANSLATION",
            term_keywords=["ribosome", "translation", "ribosomal", "rRNA", "srp"],
            gene_regex=[r"^RPL", r"^RPS", r"^EEF", r"^EIF"],
            confounder_like=True,
        ),
    ]


def assign_program(term: str, genes: List[str], cfg: ProgramSummaryConfig) -> Tuple[str, Dict]:
    """
    Deterministic assignment of a term to a program bucket.
    Returns (program_name, debug_info).
    """
    best = cfg.other_label
    best_support = 0.0
    best_dbg = {"term_hits": 0, "gene_hits": 0, "gene_frac": 0.0, "support": 0.0}

    for prog in cfg.programs:
        term_hits = _count_term_hits(term, prog.term_keywords)
        if term_hits < cfg.min_term_hit:
            # allow gene-only mapping if gene evidence is very strong
            pass

        gene_hits = _gene_family_hits(genes, prog.gene_regex)
        gene_frac = gene_hits / max(len(genes), 1)

        # normalize components
        term_component = min(1.0, term_hits / 2.0)   # 0, 0.5, 1...
        gene_component = min(1.0, gene_frac / 0.35)  # saturate around 35%

        support = cfg.w_term * term_component + cfg.w_gene * gene_component

        # require some evidence
        if support >= cfg.min_total_support and support > best_support:
            best_support = support
            best = prog.name
            best_dbg = {
                "term_hits": term_hits,
                "gene_hits": gene_hits,
                "gene_frac": float(gene_frac),
                "support": float(support),
                "confounder_like": prog.confounder_like,
            }

    return best, best_dbg


def summarize_programs(
    rows: List[Dict],
    phenotype: str = "",
    config: Optional[ProgramSummaryConfig] = None,
) -> Dict:
    """
    rows: output rows from triage layer (each has term, genes_list, scores)
    Returns program-level summary suitable for UI + GPT-5 conditioning.
    """
    if config is None:
        config = ProgramSummaryConfig(programs=default_program_rules(), phenotype=phenotype or "")

    ph = _norm(phenotype or "")
    phenotype_mentions_cycle = any(k in ph for k in ["prolifer", "cell cycle", "growth", "tumor", "regenerat"])
    phenotype_mentions_ribo = any(k in ph for k in ["translation", "ribosom", "proteosynth"])

    # group members
    program_members = defaultdict(list)
    program_debug = defaultdict(list)

    for r in rows:
        prog, dbg = assign_program(r["term"], r.get("genes_list", []), config)
        r["program_label"] = prog
        r["program_assign_debug"] = dbg
        program_members[prog].append(r)
        program_debug[prog].append(dbg)

    # build summaries
    programs_out = []
    for prog, members in program_members.items():
        # skip OTHER if it’s empty-ish? keep it but lower
        scores = [m.get("combined_pre_gpt_score", m.get("triage_score", 0.0)) for m in members]
        scores_sorted = sorted(scores, reverse=True)
        top_k = scores_sorted[: config.top_k_for_mean]
        mean_top = float(np.mean(top_k)) if top_k else 0.0
        max_score = float(scores_sorted[0]) if scores_sorted else 0.0

        prog_score = config.w_prog_max * max_score + config.w_prog_mean * mean_top

        # confounder penalty unless phenotype supports it
        conf_like = False
        for pr in config.programs:
            if pr.name == prog:
                conf_like = pr.confounder_like
                break

        if conf_like:
            if prog == "CELL_CYCLE_PROLIFERATION" and not phenotype_mentions_cycle:
                prog_score *= (1.0 - config.confounder_penalty)
            if prog == "RIBOSOME_TRANSLATION" and not phenotype_mentions_ribo:
                prog_score *= (1.0 - config.confounder_penalty)

        # driver genes: frequency-weighted, breaking ties by appearing in high-score terms
        gene_counter = Counter()
        gene_weight = Counter()

        for m in members:
            s = float(m.get("combined_pre_gpt_score", m.get("triage_score", 0.0)))
            for g in m.get("genes_list", []):
                gene_counter[g] += 1
                gene_weight[g] += s

        top_genes = sorted(gene_counter.keys(), key=lambda g: (gene_counter[g], gene_weight[g]), reverse=True)[:25]

        # representative terms: top by score
        rep_terms = sorted(members, key=lambda m: m.get("combined_pre_gpt_score", m.get("triage_score", 0.0)), reverse=True)[:10]
        rep_terms_out = [{
            "row_id": m["row_id"],
            "term": m["term"],
            "score": float(m.get("combined_pre_gpt_score", m.get("triage_score", 0.0))),
            "biofit_score": float(m.get("biofit_score", 0.0)),
            "flags": m.get("flags", []),
            "overlap_k": int(m.get("overlap_k", 0)),
        } for m in rep_terms]

        programs_out.append({
            "program": prog,
            "program_score": float(np.clip(prog_score, 0.0, 100.0)),
            "member_count": len(members),
            "top_genes": top_genes,
            "representative_terms": rep_terms_out,
        })

    programs_out = sorted(programs_out, key=lambda x: x["program_score"], reverse=True)

    return {
        "programs": programs_out,
        "meta": {
            "n_programs": len(programs_out),
            "phenotype": phenotype or "",
        }
    }

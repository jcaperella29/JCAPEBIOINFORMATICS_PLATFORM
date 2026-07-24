from __future__ import annotations

from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Literal, Optional
import hashlib
import json
import os

Role = Literal["Likely driver", "Likely reactive", "Likely artifact / confounded", "Uncertain"]
EvidenceStrength = Literal["Stronger", "Moderate", "Weak"]
LiteratureStatus = Literal["direct_support", "general_support", "conflicting", "background", "not_assessed"]

SUPPORTED_ASSAYS = {
    "bulk_rnaseq": "bulk RNA-seq",
    "scrnaseq": "scRNA-seq",
    "perturbseq": "Perturb-seq",
    "atacseq": "ATAC-seq",
    "mirnaseq": "miRNA-seq",
    "gwas": "GWAS",
    "dna_methylation": "DNA methylation",
}


def norm_assay(s: str) -> str:
    x = (s or "").strip().lower().replace("-", "").replace("_", "").replace(" ", "")
    mapping = {
        "bulkrnaseq": "bulk_rnaseq",
        "rnaseq": "bulk_rnaseq",
        "scrnaseq": "scrnaseq",
        "singlecellrnaseq": "scrnaseq",
        "perturbseq": "perturbseq",
        "atacseq": "atacseq",
        "mirnaseq": "mirnaseq",
        "gwas": "gwas",
        "dnamethylation": "dna_methylation",
        "methylation": "dna_methylation",
    }
    return mapping.get(x, x)


def stable_json_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass
class EvidenceLink:
    terms: List[str] = field(default_factory=list)
    row_ids: List[str] = field(default_factory=list)
    genes: List[str] = field(default_factory=list)
    pmids: List[str] = field(default_factory=list)
    literature_status: LiteratureStatus = "not_assessed"


@dataclass
class ValidationStep:
    experiment: str
    readout: str
    control: str
    expected_result_if_causal: str = ""
    purpose: str = "test causality or validate the interpretation"


@dataclass
class InterpretationClaim:
    claim: str
    program: str
    role: Role = "Uncertain"
    evidence_strength: EvidenceStrength = "Weak"
    rationale: str = ""
    evidence: EvidenceLink = field(default_factory=EvidenceLink)
    confounders: List[str] = field(default_factory=list)
    validation: List[ValidationStep] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class RunMetadata:
    app_version: str
    prompt_version: str
    playbook_version: str
    model_name: str
    timestamp_utc: str
    input_hash: str
    phenotype_hash: str
    context_hash: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def build_run_metadata(*, df_columns: List[str], n_rows: int, phenotype: str, context: Dict[str, Any], model_name: str) -> Dict[str, Any]:
    input_fingerprint = {"columns": df_columns, "n_rows": n_rows}
    return RunMetadata(
        app_version=os.environ.get("APP_VERSION", "0.2.0-productization"),
        prompt_version=os.environ.get("PROMPT_VERSION", "claims-v1"),
        playbook_version=os.environ.get("PLAYBOOK_VERSION", "2026-05-productization"),
        model_name=model_name,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        input_hash=stable_json_hash(input_fingerprint),
        phenotype_hash=stable_json_hash({"phenotype": phenotype}),
        context_hash=stable_json_hash(context),
    ).to_dict()


def claim_schema_for_prompt() -> Dict[str, Any]:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "headline": {"type": "string"},
            "executive_summary": {"type": "string"},
            "claims": {
                "type": "array",
                "items": {
                    "type": "object",
                    "additionalProperties": False,
                    "properties": {
                        "claim": {"type": "string"},
                        "program": {"type": "string"},
                        "role": {"type": "string", "enum": ["Likely driver", "Likely reactive", "Likely artifact / confounded", "Uncertain"]},
                        "evidence_strength": {"type": "string", "enum": ["Stronger", "Moderate", "Weak"]},
                        "rationale": {"type": "string"},
                        "supporting_terms": {"type": "array", "items": {"type": "string"}},
                        "supporting_row_ids": {"type": "array", "items": {"type": "string"}},
                        "supporting_genes": {"type": "array", "items": {"type": "string"}},
                        "supporting_pmids": {"type": "array", "items": {"type": "string"}},
                        "literature_status": {"type": "string", "enum": ["direct_support", "general_support", "conflicting", "background", "not_assessed"]},
                        "confounders": {"type": "array", "items": {"type": "string"}},
                        "validation": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "additionalProperties": False,
                                "properties": {
                                    "experiment": {"type": "string"},
                                    "readout": {"type": "string"},
                                    "control": {"type": "string"},
                                    "expected_result_if_causal": {"type": "string"},
                                    "purpose": {"type": "string"},
                                },
                                "required": ["experiment", "readout", "control", "expected_result_if_causal", "purpose"],
                            },
                        },
                    },
                    "required": [
                        "claim", "program", "role", "evidence_strength", "rationale",
                        "supporting_terms", "supporting_row_ids", "supporting_genes",
                        "supporting_pmids", "literature_status", "confounders", "validation",
                    ],
                },
            },
            "assay_limitations": {"type": "array", "items": {"type": "string"}},
            "main_confounders": {"type": "array", "items": {"type": "string"}},
            "next_best_experiment": {"type": "string"},
            "markdown_report": {"type": "string"},
        },
        "required": ["headline", "executive_summary", "claims", "assay_limitations", "main_confounders", "next_best_experiment", "markdown_report"],
    }



def _role_rank(role: str) -> int:
    role = (role or "").lower()
    if "driver" in role:
        return 4
    if "reactive" in role:
        return 3
    if "uncertain" in role:
        return 2
    if "artifact" in role or "confounded" in role:
        return 1
    return 0


def _evidence_rank(evidence_strength: str) -> int:
    x = (evidence_strength or "").lower()
    if "strong" in x:
        return 3
    if "moderate" in x:
        return 2
    if "weak" in x:
        return 1
    return 0


def _claim_priority(c: Dict[str, Any]) -> int:
    return 10 * _role_rank(c.get("role", "")) + _evidence_rank(c.get("evidence_strength", ""))


def _is_handling_or_ieg_artifact(c: Dict[str, Any]) -> bool:
    text_bits = [
        c.get("claim", ""), c.get("rationale", ""), c.get("program", ""),
        " ".join(c.get("confounders", []) or []),
        " ".join(c.get("supporting_genes", []) or []),
        " ".join(c.get("supporting_terms", []) or []),
    ]
    text = " ".join(str(x) for x in text_bits).lower()
    return any(k in text for k in ["handling", "immediate-early", "immediate early", "ieg", "fos", "jun", "media-change", "media change", "procedural stress"])


def deduplicate_claims(claims: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Keep one primary claim per exact program label.

    If GPT emits a second, lower-priority claim for the same program with a conflicting
    role, preserve the useful reasoning but relabel it as an explicit alternative
    explanation. This prevents the PDF/UI from looking contradictory while still showing
    caveats such as handling-stress / immediate-early gene artifacts.
    """
    by_program: Dict[str, List[Dict[str, Any]]] = {}
    for c in claims:
        program = (c.get("program") or "Unassigned").strip() or "Unassigned"
        by_program.setdefault(program, []).append(c)

    out: List[Dict[str, Any]] = []
    for program, group in by_program.items():
        if len(group) == 1:
            out.append(group[0])
            continue

        group_sorted = sorted(group, key=_claim_priority, reverse=True)
        primary = group_sorted[0]
        out.append(primary)

        primary_role = (primary.get("role") or "").lower()
        for dup in group_sorted[1:]:
            dup_role = (dup.get("role") or "").lower()
            # If genuinely redundant with same role/evidence, drop it.
            if dup_role == primary_role and (dup.get("evidence_strength") == primary.get("evidence_strength")):
                continue

            alt = dict(dup)
            if _is_handling_or_ieg_artifact(alt):
                alt["program"] = f"Alternative explanation: handling-stress / immediate-early gene artifact ({program})"
            else:
                alt["program"] = f"Alternative explanation: {program}"
            out.append(alt)

    return sorted(out, key=_claim_priority, reverse=True)


def normalize_claims(parsed: Dict[str, Any]) -> List[Dict[str, Any]]:
    claims = []
    for c in parsed.get("claims", []) or []:
        evidence = EvidenceLink(
            terms=list(c.get("supporting_terms", []) or []),
            row_ids=list(c.get("supporting_row_ids", []) or []),
            genes=list(c.get("supporting_genes", []) or []),
            pmids=[str(x) for x in (c.get("supporting_pmids", []) or [])],
            literature_status=c.get("literature_status") or "not_assessed",
        )
        validation = []
        for v in c.get("validation", []) or []:
            validation.append(ValidationStep(
                experiment=v.get("experiment", ""),
                readout=v.get("readout", ""),
                control=v.get("control", ""),
                expected_result_if_causal=v.get("expected_result_if_causal", ""),
                purpose=v.get("purpose", "test causality or validate the interpretation"),
            ))
        claims.append(InterpretationClaim(
            claim=c.get("claim", ""),
            program=c.get("program", ""),
            role=c.get("role") or "Uncertain",
            evidence_strength=c.get("evidence_strength") or "Weak",
            rationale=c.get("rationale", ""),
            evidence=evidence,
            confounders=list(c.get("confounders", []) or []),
            validation=validation,
        ).to_dict())
    return deduplicate_claims(claims)

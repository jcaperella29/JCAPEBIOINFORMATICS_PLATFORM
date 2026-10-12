#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import io
import json
import zipfile
from pathlib import Path


def load_json(path):
    path = Path(path)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def read_csv_rows(path, limit=None):
    path = Path(path)
    if not path.exists():
        return []

    try:
        with path.open(newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
    except Exception:
        return []

    return rows[:limit] if limit else rows


def read_zip_member(zip_path, member):
    zip_path = Path(zip_path)
    if not zip_path.exists():
        return None

    try:
        with zipfile.ZipFile(zip_path) as z:
            if member not in z.namelist():
                return None
            return z.read(member)
    except Exception:
        return None


def json_from_zip(zip_path, member):
    raw = read_zip_member(zip_path, member)
    if raw is None:
        return None

    try:
        return json.loads(raw.decode("utf-8"))
    except Exception:
        return None


def csv_from_zip(zip_path, member, limit=None):
    raw = read_zip_member(zip_path, member)
    if raw is None:
        return []

    try:
        txt = raw.decode("utf-8-sig")
        rows = list(csv.DictReader(io.StringIO(txt)))
    except Exception:
        return []

    return rows[:limit] if limit else rows


def find_target_root(target_dir):
    target_dir = Path(target_dir)

    direct = target_dir / f"target_{target_dir.name}"
    if direct.exists():
        return direct

    candidates = sorted(target_dir.glob("target_*"))
    if candidates:
        return candidates[0]

    return target_dir


def first_present(row, names):
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    return None


def summarize_mixscape(root):
    mix = root / "01_mixscape"

    provenance = load_json(mix / "target_provenance.json") or {}
    run_summary = load_json(mix / "run_summary.json") or {}
    handoff = load_json(mix / "crispr_mixscape_literature_handoff.json") or {}

    de_rows = read_csv_rows(mix / "ko_genes.csv", limit=10)
    enrichment_rows = read_csv_rows(mix / "enrichment_all.csv", limit=10)

    perturbation_summary = handoff.get("perturbation_summary") or []

    selected = None
    selected_target = (
        (handoff.get("context") or {}).get("selected_perturbation")
        or provenance.get("selected_perturbation")
        or provenance.get("target")
    )

    if isinstance(perturbation_summary, list):
        for row in perturbation_summary:
            if not isinstance(row, dict):
                continue
            if str(row.get("perturbation", "")).upper() == str(selected_target or "").upper():
                selected = row
                break

    return {
        "selected_target": selected_target,
        "classification_summary": selected,
        "analysis_summary": handoff.get("analysis_summary") or run_summary,
        "top_de": de_rows,
        "top_enrichment": enrichment_rows,
    }


def summarize_network(root):
    net = root / "02_network"

    consensus = load_json(net / "consensus.json") or {}
    consensus_rows = read_csv_rows(net / "consensus_candidates.csv", limit=10)
    bip_rows = read_csv_rows(net / "bipartite_top_candidates.csv", limit=10)
    proj_rows = read_csv_rows(net / "projection_top_candidates.csv", limit=10)

    return {
        "consensus": consensus,
        "consensus_candidates": consensus_rows,
        "bipartite_candidates": bip_rows,
        "projection_candidates": proj_rows,
    }


def summarize_llm(root):
    z = root / "03_llm_triage" / "llm_bundle.zip"

    analysis = json_from_zip(z, "analysis_result.json") or {}
    claims = csv_from_zip(z, "claims.csv", limit=10)
    triage = csv_from_zip(z, "triage_rows.csv", limit=10)
    programs = csv_from_zip(z, "programs.csv", limit=10)

    return {
        "analysis_result": analysis,
        "claims": claims,
        "triage_rows": triage,
        "programs": programs,
    }


def summarize_literature(root):
    lit = root / "04_literature_review" / "literature_review_response.json"
    response = load_json(lit) or {}

    assessments = response.get("assessments") or []

    by_level = {}
    for a in assessments:
        if not isinstance(a, dict):
            continue
        level = str(a.get("level", "unknown"))
        by_level[level] = by_level.get(level, 0) + 1

    ranked = sorted(
        [a for a in assessments if isinstance(a, dict)],
        key=lambda x: float(x.get("relevance_score") or 0),
        reverse=True,
    )[:10]

    compact = []

    for a in ranked:
        paper = a.get("paper") or {}

        compact.append({
            "finding_id": a.get("finding_id"),
            "level": a.get("level"),
            "relevance_score": a.get("relevance_score"),
            "title": paper.get("title"),
            "year": paper.get("year"),
            "journal": paper.get("journal"),
            "pmid": paper.get("pmid"),
            "doi": paper.get("doi"),
            "same_target_or_question": a.get("same_target_or_question"),
            "assay_match": a.get("assay_match"),
            "organism_match": a.get("organism_match"),
            "context_match": a.get("context_match"),
            "reasons": a.get("reasons") or [],
        })

    return {
        "run_id": response.get("run_id"),
        "summary_counts": response.get("summary_counts") or {},
        "evidence_levels": by_level,
        "top_assessments": compact,
        "n_assessments": len(assessments),
    }


def best_label(row):
    if not isinstance(row, dict):
        return None

    return first_present(
        row,
        [
            "term_name",
            "name",
            "description",
            "term",
            "pathway",
            "label",
            "candidate",
            "node",
            "gene",
            "symbol",
        ],
    )


def best_gene(row):
    if not isinstance(row, dict):
        return None

    return first_present(
        row,
        ["gene", "symbol", "feature", "Gene", "gene_symbol"],
    )


def markdown_target(target, data):
    mix = data["mixscape"]
    network = data["network"]
    llm = data["triage"]
    literature = data["literature"]

    lines = [f"## {target}", ""]

    cls = mix.get("classification_summary")
    if isinstance(cls, dict):
        ko = first_present(cls, ["ko_fraction", "KO_fraction", "fraction_ko"])
        np = first_present(cls, ["np_fraction", "NP_fraction", "fraction_np"])

        bits = []
        if ko is not None:
            bits.append(f"KO fraction: {ko}")
        if np is not None:
            bits.append(f"NP fraction: {np}")

        if bits:
            lines.append("**Mixscape:** " + "; ".join(bits))
            lines.append("")

    de = [
        best_gene(r)
        for r in mix.get("top_de", [])
        if best_gene(r)
    ][:5]

    if de:
        lines.append("**Top DE genes:** " + ", ".join(de))
        lines.append("")

    enrich = [
        best_label(r)
        for r in mix.get("top_enrichment", [])
        if best_label(r)
    ][:5]

    if enrich:
        lines.append("**Top enriched programs:** " + "; ".join(enrich))
        lines.append("")

    candidates = [
        best_label(r)
        for r in network.get("consensus_candidates", [])
        if best_label(r)
    ][:5]

    if candidates:
        lines.append("**Top network candidates:** " + ", ".join(candidates))
        lines.append("")

    claims = llm.get("claims") or []

    if claims:
        summaries = []

        for row in claims[:3]:
            text = first_present(
                row,
                ["claim", "statement", "summary", "text"],
            )
            if text:
                summaries.append(str(text))

        if summaries:
            lines.append("**Triage highlights:**")
            for text in summaries:
                lines.append(f"- {text}")
            lines.append("")

    lit = literature

    lines.append(
        f"**Literature evidence:** {lit.get('n_assessments', 0)} assessments"
    )

    if lit.get("evidence_levels"):
        levels = ", ".join(
            f"{k}: {v}"
            for k, v in sorted(lit["evidence_levels"].items())
        )
        lines.append(f"- Evidence levels: {levels}")

    for paper in lit.get("top_assessments", [])[:3]:
        title = paper.get("title")
        level = paper.get("level")
        score = paper.get("relevance_score")

        if title:
            lines.append(
                f"- {title} [{level}; relevance={score}]"
            )

    lines.append("")
    return lines


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--targets-dir", required=True)
    ap.add_argument("--context", required=True)
    ap.add_argument("--cell-state-handoff")
    ap.add_argument("--outdir", required=True)

    args = ap.parse_args()

    targets_dir = Path(args.targets_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    context = load_json(args.context) or {}
    cell_state = (
        load_json(args.cell_state_handoff)
        if args.cell_state_handoff
        else None
    ) or {}

    targets = {}

    for target_dir in sorted(targets_dir.iterdir()):
        if not target_dir.is_dir():
            continue

        target = target_dir.name
        root = find_target_root(target_dir)

        targets[target] = {
            "mixscape": summarize_mixscape(root),
            "network": summarize_network(root),
            "triage": summarize_llm(root),
            "literature": summarize_literature(root),
        }

    final = {
        "schema_version": "jcap-perturbseq-final-report-v1",
        "experiment": context,
        "cell_state": cell_state,
        "target_count": len(targets),
        "targets": targets,
    }

    json_path = outdir / "perturbseq_final_report.json"
    json_path.write_text(
        json.dumps(final, indent=2),
        encoding="utf-8",
    )

    md = [
        "# JCAP Perturb-seq Integrated Final Report",
        "",
        f"Targets summarized: {len(targets)}",
        "",
        "## Experiment Context",
        "",
        "```json",
        json.dumps(context, indent=2),
        "```",
        "",
        "## Cell-State Context",
        "",
    ]

    if cell_state:
        md.extend([
            "```json",
            json.dumps(cell_state, indent=2),
            "```",
            "",
        ])
    else:
        md.extend([
            "No cell-state handoff was available.",
            "",
        ])

    md.extend([
        "## Cross-target Summary",
        "",
        "| Target | Literature assessments |",
        "|---|---:|",
    ])

    for target, data in targets.items():
        md.append(
            f"| {target} | {data['literature'].get('n_assessments', 0)} |"
        )

    md.append("")

    for target, data in targets.items():
        md.extend(markdown_target(target, data))

    md.extend([
        "## Interpretation Notes",
        "",
        "- Mixscape classifications estimate perturbation-responsive versus non-perturbed cells.",
        "- Differential expression and enrichment summarize downstream transcriptional consequences and do not by themselves establish direct molecular interaction.",
        "- Network and LLM outputs are prioritization aids rather than independent experimental validation.",
        "- Literature evidence is ranked contextual support and should be interpreted with assay, organism, and biological-context matching in mind.",
        "",
    ])

    md_path = outdir / "perturbseq_final_report.md"
    md_path.write_text("\n".join(md), encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()

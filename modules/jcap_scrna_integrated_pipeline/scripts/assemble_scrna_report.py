#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def jload(p):
    p = Path(p)
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def csvload(p):
    p = Path(p)
    if not p.exists():
        return []
    with p.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def find(root, name):
    hits = list(Path(root).rglob(name))
    return hits[0] if hits else None


def first(row, keys):
    for k in keys:
        v = row.get(k)
        if v not in (None, ""):
            return str(v).strip()
    return None


def number(v):
    try:
        return float(v)
    except Exception:
        return None



def build_finding_labels(obj):
    labels = {}

    def walk(x):
        if isinstance(x, dict):
            fid = x.get("finding_id") or x.get("id")

            if isinstance(fid, str) and fid.startswith("finding_"):
                label = None

                for key in (
                    "claim",
                    "finding",
                    "finding_text",
                    "question",
                    "summary",
                    "description",
                    "text",
                    "title",
                    "hypothesis",
                ):
                    val = x.get(key)

                    if isinstance(val, str) and val.strip():
                        label = val.strip()
                        break

                if label:
                    labels[fid] = label

            for v in x.values():
                walk(v)

        elif isinstance(x, list):
            for v in x:
                walk(v)

    walk(obj)
    return labels


def literature_markdown(request_obj, response_obj):
    assessments = response_obj.get("assessments") or []

    if not isinstance(assessments, list) or not assessments:
        return ["_No literature assessments were returned._"]

    labels = build_finding_labels(request_obj)

    grouped = {}

    for a in assessments:
        if not isinstance(a, dict):
            continue

        fid = str(a.get("finding_id") or "unassigned")
        grouped.setdefault(fid, []).append(a)

    out = []

    def finding_sort_key(fid):
        try:
            return int(fid.rsplit("_", 1)[-1])
        except Exception:
            return 10**9

    for fid in sorted(grouped, key=finding_sort_key):
        rows = grouped[fid]

        label = labels.get(fid)

        if label:
            out.append(f"### {label}")
            out.append("")
            out.append(f"*Literature target: `{fid}`*")
        else:
            pretty = fid.replace("_", " ").title()
            out.append(f"### {pretty}")

        out.append("")

        counts = {}
        for row in rows:
            level = str(row.get("level") or "unclassified")
            counts[level] = counts.get(level, 0) + 1

        count_text = ", ".join(
            f"**{k}**: {v}"
            for k, v in sorted(counts.items())
        )

        out.append(
            f"**Evidence screen:** {len(rows)} paper"
            + ("s" if len(rows) != 1 else "")
            + (f" ({count_text})" if count_text else "")
        )
        out.append("")

        # Strongest / most relevant papers first.
        def paper_sort(a):
            level = str(a.get("level") or "").lower()

            level_rank = {
                "direct": 0,
                "support": 1,
                "supporting": 1,
                "indirect": 2,
                "context": 3,
                "exclude": 9,
            }.get(level, 5)

            try:
                rel = float(a.get("relevance_score"))
            except Exception:
                rel = 0.0

            return (level_rank, -rel)

        for a in sorted(rows, key=paper_sort):
            paper = a.get("paper") or {}

            title = paper.get("title") or "Untitled paper"
            year = paper.get("year")
            journal = paper.get("journal")
            pmid = paper.get("pmid")
            level = a.get("level") or "unclassified"
            score = a.get("relevance_score")

            bits = [f"**{title}**"]

            meta = []
            if year:
                meta.append(str(year))
            if journal:
                meta.append(str(journal))
            if pmid:
                meta.append(f"PMID {pmid}")

            if meta:
                bits.append(" — " + "; ".join(meta))

            out.append("- " + "".join(bits))

            detail = f"  - Evidence level: **{level}**"

            if score is not None:
                try:
                    detail += f"; relevance score: **{float(score):.2f}**"
                except Exception:
                    detail += f"; relevance score: **{score}**"

            out.append(detail)

            reasons = a.get("reasons") or []

            if isinstance(reasons, list):
                for reason in reasons:
                    if reason:
                        out.append(f"  - {reason}")

            caution = a.get("caution")

            if caution:
                out.append(f"  - Caution: {caution}")

        out.append("")

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scrna-dir", required=True)
    ap.add_argument("--network-dir", required=True)
    ap.add_argument("--llm-dir", required=True)
    ap.add_argument("--literature-dir", required=True)
    ap.add_argument("--context-json", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ctx = jload(a.context_json)

    run_summary_path = find(a.scrna_dir, "run_summary.json")
    run_summary = jload(run_summary_path) if run_summary_path else {}

    condition_only_path = find(a.scrna_dir, "condition_only_de_table.csv")
    condition_de_path = find(a.scrna_dir, "condition_de_table.csv")

    de_rows = (
        csvload(condition_only_path)
        if condition_only_path
        else csvload(condition_de_path) if condition_de_path else []
    )

    network_path = (
        find(a.network_dir, "consensus_candidates.csv")
        or find(a.network_dir, "projection_top_candidates.csv")
        or find(a.network_dir, "bipartite_top_candidates.csv")
    )
    network_rows = csvload(network_path) if network_path else []

    claims_path = find(a.llm_dir, "claims.csv")
    claims_rows = csvload(claims_path) if claims_path else []

    programs_path = find(a.llm_dir, "programs.csv")
    programs_rows = csvload(programs_path) if programs_path else []

    llm_report_path = find(a.llm_dir, "report.md")
    llm_report = (
        llm_report_path.read_text(encoding="utf-8")
        if llm_report_path else ""
    )

    lit_request_path = find(a.literature_dir, "literature_review_request.json")
    lit_request = jload(lit_request_path) if lit_request_path else {}

    lit_path = find(a.literature_dir, "literature_review_response.json")
    lit = jload(lit_path) if lit_path else {}

    cs = ctx.get("cell_state_analysis", {})
    state_rows = cs.get("group_state_summary") or []

    groups = sorted({
        str(r.get("group"))
        for r in state_rows
        if r.get("group") is not None
    })

    states = []
    for r in state_rows:
        s = r.get("state")
        if s and s not in states:
            states.append(s)

    by_state = {}
    for r in state_rows:
        by_state[(str(r.get("group")), str(r.get("state")))] = r

    dataset = ctx.get("dataset", "scRNA-seq dataset")
    organism = ctx.get("organism", "")
    assay = ctx.get("assay", "single-cell RNA-seq")
    study = ctx.get("study_context", "")
    target_ct = ctx.get("target_cell_type", "")
    comparison = ctx.get("comparison") or {}
    ca = comparison.get("condition_a", "")
    cb = comparison.get("condition_b", "")

    n_cells = (
        run_summary.get("n_cells")
        or cs.get("n_cells")
    )
    n_genes = run_summary.get("n_genes")

    claims = []
    for row in claims_rows:
        text = first(row, [
            "claim", "claim_text", "summary", "finding",
            "description", "interpretation", "rationale"
        ])
        if text and text not in claims:
            claims.append(text)

    programs = []
    for row in programs_rows:
        name = first(row, [
            "program", "pathway", "term", "name", "label"
        ])
        summary = first(row, [
            "summary", "description", "interpretation", "rationale"
        ])
        if name or summary:
            programs.append((name, summary))

    def de_gene(row):
        g = first(row, [
            "gene", "Gene", "symbol", "Symbol",
            "feature", "Feature"
        ])
        if g:
            return g

        for k, v in row.items():
            if k in (None, "", "X") and v:
                return str(v)

        return None

    def de_fc(row):
        for k in ["avg_log2FC", "logFC", "log2FoldChange"]:
            x = number(row.get(k))
            if x is not None:
                return x
        return None

    def de_p(row):
        for k in ["p_val_adj", "adj.P.Val", "padj", "FDR"]:
            x = number(row.get(k))
            if x is not None:
                return x
        return None

    de_top = []
    for row in sorted(
        de_rows,
        key=lambda r: de_p(r) if de_p(r) is not None else 1e99
    ):
        g = de_gene(row)
        if g:
            de_top.append((g, de_fc(row), de_p(row)))
        if len(de_top) >= 12:
            break

    network_top = []
    for row in network_rows[:10]:
        name = first(row, [
            "gene", "Gene", "candidate", "Candidate",
            "node", "Node", "name", "Name",
            "symbol", "Symbol"
        ])

        if not name:
            name = next(
                (str(v) for v in row.values() if v not in (None, "")),
                "unknown"
            )

        score = first(row, [
            "consensus_score", "score", "Score",
            "rank_score", "diffusion_score", "projection_score"
        ])

        network_top.append((name, score))

    shifts = []

    if len(groups) == 2:
        for state in states:
            x = by_state.get((groups[0], state), {})
            y = by_state.get((groups[1], state), {})

            fx = x.get("fraction_within_group")
            fy = y.get("fraction_within_group")

            if fx is not None and fy is not None:
                shifts.append(
                    (state, 100 * (float(fy) - float(fx)))
                )

        shifts.sort(key=lambda x: abs(x[1]), reverse=True)

    md = []

    md.append("# Integrated scRNA-seq Analysis")
    md.append("")
    md.append("## Executive Summary")
    md.append("")

    summary_bits = []

    summary_bits.append(
        f"This analysis integrates differential expression, enrichment, "
        f"cell-state composition, network prioritization, LLM biological "
        f"triage, and literature review for **{dataset}**."
    )

    if n_cells:
        summary_bits.append(
            f"The dataset contains **{int(n_cells):,} cells**"
            + (f" and **{int(n_genes):,} genes**." if n_genes else ".")
        )

    if ca and cb:
        summary_bits.append(
            f"The experimental comparison is **{ca} vs {cb}**."
        )

    if target_ct:
        summary_bits.append(
            f"The configured biological focus is **{target_ct}**."
        )

    if shifts:
        top = ", ".join(
            f"{name} ({delta:+.2f} percentage points)"
            for name, delta in shifts[:3]
        )
        summary_bits.append(
            "The largest cell-state composition differences were "
            + top + "."
        )

    if claims:
        summary_bits.append(
            f"LLM triage produced **{len(claims)} structured biological claims**."
        )

    assessments = lit.get("assessments") or []

    if isinstance(assessments, list):
        summary_bits.append(
            f"The literature stage produced **{len(assessments)} evidence assessments**."
        )

    md.extend(summary_bits)

    md.append("")
    md.append("### PI-facing Biological Synthesis")
    md.append("")

    if claims:
        for claim in claims[:8]:
            md.append(f"- {claim}")
    elif llm_report:
        md.append(llm_report[:5000])
    else:
        md.append("_No LLM synthesis was available._")

    md.append("")
    md.append("## Experimental Context")
    md.append("")
    md.append(f"- **Dataset:** {dataset}")
    md.append(f"- **Organism:** {organism}")
    md.append(f"- **Assay:** {assay}")

    if study:
        md.append(f"- **Study context:** {study}")

    if ca or cb:
        md.append(f"- **Comparison:** {ca} vs {cb}")

    if target_ct:
        md.append(f"- **Target cell type:** {target_ct}")

    md.append("")
    md.append("## Cell-State / Lineage Context")
    md.append("")
    md.append(
        "Condition-aware composition is shown below so transcriptional "
        "changes can be interpreted alongside possible cell-state shifts."
    )
    md.append("")

    if groups and states:
        hdr = ["State"]
        for g in groups:
            hdr += [f"{g} n", f"{g} %"]

        if len(groups) == 2:
            hdr += ["Δ percentage points"]

        md.append("| " + " | ".join(hdr) + " |")
        md.append("|" + "|".join(
            ["---"] + ["---:" for _ in hdr[1:]]
        ) + "|")

        for state in states:
            row = [state]

            for g in groups:
                d = by_state.get((g, state), {})
                n = d.get("n_cells", "")
                f = d.get("fraction_within_group")

                row.append(str(n))
                row.append(
                    f"{100*float(f):.1f}%"
                    if f is not None else ""
                )

            if len(groups) == 2:
                x = by_state.get((groups[0], state), {})
                y = by_state.get((groups[1], state), {})

                fx = x.get("fraction_within_group")
                fy = y.get("fraction_within_group")

                row.append(
                    f"{100*(float(fy)-float(fx)):+.2f}"
                    if fx is not None and fy is not None else ""
                )

            md.append("| " + " | ".join(row) + " |")
    else:
        md.append("_Condition-specific cell-state table unavailable._")

    md.append("")
    md.append("## Differential Expression Across Conditions")
    md.append("")

    if de_top:
        md.append("| Gene | log2FC | adjusted p |")
        md.append("|---|---:|---:|")

        for gene, fc, p in de_top:
            md.append(
                f"| {gene} | "
                f"{'' if fc is None else f'{fc:.3f}'} | "
                f"{'' if p is None else f'{p:.3g}'} |"
            )
    else:
        md.append("_No DE table available._")

    md.append("")
    md.append("## Cross-Condition Biological Themes")
    md.append("")

    if programs:
        for name, summary in programs[:12]:
            if name and summary:
                md.append(f"- **{name}** — {summary}")
            elif name:
                md.append(f"- **{name}**")
            elif summary:
                md.append(f"- {summary}")
    elif claims:
        for claim in claims[:8]:
            md.append(f"- {claim}")

    md.append("")
    md.append("## Condition-Level Findings")
    md.append("")
    md.append(f"### {ca or 'Condition A'} vs {cb or 'Condition B'}")
    md.append("")
    md.append("#### Biological interpretation")
    md.append("")

    for claim in claims[:8]:
        md.append(f"- {claim}")

    md.append("")
    md.append("#### Key supporting terms/genes")
    md.append("")

    if de_top:
        md.append(
            "- **Leading DE genes:** "
            + ", ".join(x[0] for x in de_top[:10])
        )

    if programs:
        names = [x[0] for x in programs if x[0]]
        if names:
            md.append(
                "- **Programs/pathways:** "
                + ", ".join(names[:10])
            )

    md.append("")
    md.append("#### Top confounders to monitor")
    md.append("")
    md.append(
        "- Cell-state composition differences can contribute to apparent "
        "pathway changes."
    )
    md.append(
        "- Cell-level DE is exploratory for biological inference; "
        "replicate-aware pseudobulk is preferable for formal inference."
    )
    md.append(
        "- Enrichment and network ranking prioritize hypotheses rather "
        "than establish direct causality."
    )

    md.append("")
    md.append("## Network Prioritization")
    md.append("")

    if network_top:
        md.append("| Candidate | Score |")
        md.append("|---|---:|")

        for name, score in network_top:
            md.append(f"| {name} | {score or ''} |")
    else:
        md.append("_No network candidate table available._")

    md.append("")
    md.append("## Literature Evidence")
    md.append("")
    md.extend(literature_markdown(lit_request, lit))
    md.append("")

    md.append("## Supporting Material")
    md.append("")
    md.append("- `01_scrna/` — scRNA analysis, DE, enrichment, classifier, PCA/UMAP")
    md.append("- `02_network/` — network and candidate prioritization")
    md.append("- `03_llm_triage/` — LLM claims, programs, and report")
    md.append("- `04_literature_review/` — literature request and evidence response")
    md.append("- `cell_state/` — condition-aware state context")
    md.append("")
    md.append(
        "Numerical source tables remain the authoritative record for "
        "individual statistics."
    )

    report = "\n".join(md).rstrip() + "\n"

    (out / "scrna_final_report.md").write_text(
        report,
        encoding="utf-8"
    )

    payload = {
        "dataset": dataset,
        "context": ctx,
        "run_summary": run_summary,
        "cell_state_rows": state_rows,
        "leading_de": [
            {"gene": g, "log2fc": fc, "adjusted_p": p}
            for g, fc, p in de_top
        ],
        "network_candidates": network_top,
        "llm_claims": claims,
        "llm_programs": programs,
        "literature_assessments": assessments,
    }

    (out / "scrna_final_report.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8"
    )

    print(out / "scrna_final_report.md")
    print(out / "scrna_final_report.json")


if __name__ == "__main__":
    main()

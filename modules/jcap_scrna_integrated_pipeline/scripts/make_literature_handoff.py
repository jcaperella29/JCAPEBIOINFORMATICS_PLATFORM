#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def first_value(row, names):
    for name in names:
        v = row.get(name)
        if v is not None and str(v).strip():
            return str(v).strip()
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--llm-dir", required=True)
    ap.add_argument("--context", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    llm = Path(args.llm_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    context = json.loads(Path(args.context).read_text())

    targets = []

    # Prefer explicit LLM claims.
    claims = next(llm.rglob("claims.csv"), None)
    if claims:
        with claims.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                summary = first_value(
                    row,
                    [
                        "claim",
                        "claim_text",
                        "summary",
                        "finding",
                        "description",
                        "interpretation",
                        "rationale",
                    ],
                )
                gene = first_value(
                    row,
                    ["gene", "symbol", "target_gene"],
                )
                pathway = first_value(
                    row,
                    ["pathway", "program", "term", "name"],
                )

                if summary or gene or pathway:
                    item = {}
                    if summary:
                        item["summary"] = summary
                    if gene:
                        item["gene"] = gene
                    if pathway:
                        item["pathway"] = pathway
                    targets.append(item)

    # Also capture named LLM programs if we still have room.
    programs = next(llm.rglob("programs.csv"), None)
    if programs and len(targets) < 20:
        with programs.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                pathway = first_value(
                    row,
                    ["program", "pathway", "term", "name", "label"],
                )
                summary = first_value(
                    row,
                    ["summary", "description", "interpretation", "rationale"],
                )

                if pathway or summary:
                    item = {}
                    if pathway:
                        item["pathway"] = pathway
                    if summary:
                        item["summary"] = summary
                    targets.append(item)

                if len(targets) >= 20:
                    break

    # Remove exact duplicates while preserving order.
    dedup = []
    seen = set()

    for item in targets:
        key = json.dumps(item, sort_keys=True)
        if key not in seen:
            seen.add(key)
            dedup.append(item)

    targets = dedup[:20]

    if not targets:
        raise SystemExit(
            "No literature targets could be extracted from claims.csv or programs.csv"
        )

    handoff = {
        "source_api": "enrichment_triage",
        "context": context,
        "literature_targets": targets,
    }

    dest = out / "enrichment_triage_literature_handoff.json"
    dest.write_text(json.dumps(handoff, indent=2) + "\n")

    print(f"Wrote {dest}")
    print(f"Literature targets: {len(targets)}")


if __name__ == "__main__":
    main()

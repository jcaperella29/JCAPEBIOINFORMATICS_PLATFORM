#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import shutil
import subprocess
from pathlib import Path
from datetime import datetime, timezone

def find_enrichment(root: Path) -> Path:
    files = list(root.rglob("*.csv"))
    preferred = ["enrichment_all.csv"]
    for name in preferred:
        hits = [p for p in files if p.name == name]
        if hits:
            return sorted(hits)[0]
    hits = sorted(p for p in files if "enrichment_all" in p.name.lower() and "network_edges" not in p.name.lower())
    if hits:
        return hits[0]
    hits = sorted(p for p in files if "network_edges" in p.name.lower() and "all" in p.name.lower())
    if hits:
        return hits[0]
    raise SystemExit("No RNA-seq enrichment CSV was found for LLM triage.")

def as_bool_text(v: str) -> str:
    return "true" if str(v).lower() in {"1", "true", "yes", "y"} else "false"

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--api-url", required=True)
    p.add_argument("--rnaseq-dir", required=True)
    p.add_argument("--context-json", required=True)
    p.add_argument("--mode", default="full")
    p.add_argument("--make-pdf", default="true")
    p.add_argument("--outdir", required=True)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    enrichment = find_enrichment(Path(a.rnaseq_dir))
    context = json.loads(Path(a.context_json).read_text(encoding="utf-8"))

    phenotype = str(context.pop("phenotype", context.pop("experiment_context", ""))).strip()
    if not phenotype:
        raise SystemExit("context_json must include a non-empty 'phenotype' or 'experiment_context' value.")

    standard = {
        "organism": str(context.pop("organism", "")),
        "assay": str(context.pop("assay", "")),
        "tissue": str(context.pop("tissue", "")),
        "cell_type": str(context.pop("cell_type", "")),
        "perturbation": str(context.pop("perturbation", "")),
        "timepoint": str(context.pop("timepoint", "")),
    }
    extra_context = json.dumps(context)

    input_format = "long_edges" if "network_edges" in enrichment.name.lower() else "auto"
    bundle = out / "llm_triage_bundle.zip"
    cmd = [
        "curl", "--fail-with-body", "--silent", "--show-error",
        "-X", "POST", f"{a.api_url.rstrip('/')}/analyze-bundle",
        "-F", f"file=@{enrichment}",
        "-F", f"phenotype={phenotype}",
        "-F", f"organism={standard['organism']}",
        "-F", f"assay={standard['assay']}",
        "-F", f"tissue={standard['tissue']}",
        "-F", f"cell_type={standard['cell_type']}",
        "-F", f"perturbation={standard['perturbation']}",
        "-F", f"timepoint={standard['timepoint']}",
        "-F", f"extra_context_json={extra_context}",
        "-F", f"input_format={input_format}",
        "-F", f"mode={a.mode}",
        "-F", f"make_pdf={as_bool_text(a.make_pdf)}",
        "-o", str(bundle),
    ]
    subprocess.run(cmd, check=True)

    extracted = out / "results"
    extracted.mkdir(exist_ok=True)
    shutil.unpack_archive(str(bundle), str(extracted))
    (out / "selected_input.txt").write_text(str(enrichment), encoding="utf-8")
    (out / "submitted_context.json").write_text(json.dumps({
        "phenotype": phenotype, **standard, **context
    }, indent=2), encoding="utf-8")

    manifest = {
        "stage": "llm_triage",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "api_url": a.api_url,
        "selected_rnaseq_input": str(enrichment),
        "mode": a.mode,
        "make_pdf": as_bool_text(a.make_pdf) == "true",
    }
    (out / "pipeline_stage_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

def request_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--api-url", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--phenotype", required=True)
    p.add_argument("--phenotype-column", required=True)
    p.add_argument("--species", default="hsapiens")
    p.add_argument("--contrast", default="")
    p.add_argument("--reference-group", default="")
    p.add_argument("--fdr", default="0.05")
    p.add_argument("--logfc-cutoff", default="1.0")
    p.add_argument("--enrich-backend", default="auto")
    p.add_argument("--gprof-sources", default="GO:BP,GO:MF,GO:CC,KEGG,REAC")
    p.add_argument("--enrichr-db", default="GO_Biological_Process_2023")
    p.add_argument("--classifier-validation", default="cv")
    p.add_argument("--classifier-model", default="auto")
    p.add_argument("--poll-seconds", type=int, default=5)
    p.add_argument("--timeout-minutes", type=int, default=240)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    raw_zip = out / "rnaseq_api_bundle.zip"

    cmd = [
        "curl", "--fail-with-body", "--silent", "--show-error",
        "-X", "POST", f"{a.api_url.rstrip('/')}/rnaseq/jobs",
        "-F", f"counts=@{a.counts}",
        "-F", f"phenotype=@{a.phenotype}",
        "-F", f"phenotype_column={a.phenotype_column}",
        "-F", f"species={a.species}",
        "-F", "cli_script=rnaseq_cli.R",
        "-F", f"fdr={a.fdr}",
        "-F", f"logfc_cutoff={a.logfc_cutoff}",
        "-F", f"enrich_backend={a.enrich_backend}",
        "-F", f"gprof_sources={a.gprof_sources}",
        "-F", f"enrichr_db={a.enrichr_db}",
        "-F", f"classifier_validation={a.classifier_validation}",
        "-F", f"classifier_model={a.classifier_model}",
        "-F", "zip_outputs=true",
    ]
    if a.contrast:
        cmd += ["-F", f"contrast={a.contrast}"]
    if a.reference_group:
        cmd += ["-F", f"reference_group={a.reference_group}"]

    submitted = subprocess.run(cmd, check=True, text=True, capture_output=True)
    job = json.loads(submitted.stdout)
    job_id = job["job_id"]
    (out / "submission.json").write_text(json.dumps(job, indent=2), encoding="utf-8")

    deadline = time.time() + a.timeout_minutes * 60
    status = job
    while time.time() < deadline:
        status = request_json(f"{a.api_url.rstrip('/')}/rnaseq/jobs/{job_id}")
        (out / "status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
        if status.get("status") == "completed":
            break
        if status.get("status") == "failed":
            raise SystemExit(f"RNA-seq job failed: {json.dumps(status, indent=2)}")
        time.sleep(a.poll_seconds)
    else:
        raise SystemExit(f"RNA-seq job timed out after {a.timeout_minutes} minutes: {job_id}")

    subprocess.run([
        "curl", "--fail-with-body", "--silent", "--show-error", "-L",
        f"{a.api_url.rstrip('/')}/rnaseq/jobs/{job_id}/download",
        "-o", str(raw_zip)
    ], check=True)

    extracted = out / "results"
    extracted.mkdir(exist_ok=True)
    shutil.unpack_archive(str(raw_zip), str(extracted))

    manifest = {
        "stage": "rnaseq",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "job_id": job_id,
        "api_url": a.api_url,
        "counts": Path(a.counts).name,
        "phenotype": Path(a.phenotype).name,
        "phenotype_column": a.phenotype_column,
        "species": a.species,
        "contrast": a.contrast or None,
        "reference_group": a.reference_group or None,
        "bundle": raw_zip.name,
    }
    (out / "pipeline_stage_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()

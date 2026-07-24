from __future__ import annotations

import json
import shutil
import subprocess
import uuid
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse
from pydantic import BaseModel

ROOT = Path(__file__).resolve().parents[1]
JOBS = ROOT / "jobs"
R_SCRIPT = ROOT / "R" / "crispr_mixscape_cli.R"

app = FastAPI(title="JCAP CRISPR Mixscape API", version="0.1.0")


class JobSubmitResponse(BaseModel):
    job_id: str
    status: str


def job_dir(job_id: str) -> Path:
    return JOBS / job_id


def read_status(job_id: str) -> dict:
    p = job_dir(job_id) / "status.json"
    if not p.exists():
        return {"status": "queued", "message": "Job submitted."}
    return json.loads(p.read_text())


@app.get("/health")
def health():
    return {"status": "ok", "module": "crispr_mixscape"}


@app.post("/jobs", response_model=JobSubmitResponse)
async def submit_job(
    counts: UploadFile = File(...),
    metadata: UploadFile = File(...),
    species: str = Form("hsapiens"),
    neighbors: int = Form(20),
    min_de_genes: int = Form(3),
    iter_num: int = Form(20),
    ko_label: Optional[str] = Form(None),
):
    job_id = uuid.uuid4().hex
    d = job_dir(job_id)
    d.mkdir(parents=True, exist_ok=False)

    counts_path = d / "counts.csv"
    metadata_path = d / "metadata.csv"
    outdir = d / "output"
    outdir.mkdir(exist_ok=True)

    with counts_path.open("wb") as f:
        shutil.copyfileobj(counts.file, f)
    with metadata_path.open("wb") as f:
        shutil.copyfileobj(metadata.file, f)

    cmd = [
        "Rscript",
        str(R_SCRIPT),
        "--counts", str(counts_path),
        "--metadata", str(metadata_path),
        "--outdir", str(outdir),
        "--species", species,
        "--neighbors", str(neighbors),
        "--min_de_genes", str(min_de_genes),
        "--iter_num", str(iter_num),
    ]
    if ko_label:
        cmd.extend(["--ko_label", ko_label])

    log = (d / "run.log").open("wb")
    subprocess.Popen(cmd, cwd=str(ROOT), stdout=log, stderr=subprocess.STDOUT)

    (d / "status.json").write_text(json.dumps({"status": "running", "message": "Job started."}, indent=2))
    return JobSubmitResponse(job_id=job_id, status="running")


@app.get("/jobs/{job_id}")
def get_job(job_id: str):
    d = job_dir(job_id)
    if not d.exists():
        raise HTTPException(status_code=404, detail="Job not found.")

    status = read_status(job_id)
    manifest_path = d / "output" / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None

    return {
        "job_id": job_id,
        "status": status,
        "manifest": manifest,
        "bundle_url": f"/jobs/{job_id}/bundle" if manifest_path.exists() else None,
    }


@app.get("/jobs/{job_id}/files/{path:path}")
def get_file(job_id: str, path: str):
    f = job_dir(job_id) / "output" / path
    if not f.exists() or not f.is_file():
        raise HTTPException(status_code=404, detail="File not found.")
    return FileResponse(f)


@app.get("/jobs/{job_id}/bundle")
def get_bundle(job_id: str):
    d = job_dir(job_id)
    out = d / "output"
    if not out.exists():
        raise HTTPException(status_code=404, detail="Output not found.")

    zip_path = d / "crispr_mixscape_bundle.zip"
    if zip_path.exists():
        zip_path.unlink()

    shutil.make_archive(str(zip_path.with_suffix("")), "zip", out)
    return FileResponse(zip_path, filename="crispr_mixscape_bundle.zip")

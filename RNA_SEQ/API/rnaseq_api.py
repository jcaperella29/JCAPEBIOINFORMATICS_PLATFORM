from __future__ import annotations


import json
import os
import shutil
import subprocess
import uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from fastapi import BackgroundTasks, FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse
from pydantic import BaseModel


APP_ROOT = Path(os.getenv("APP_ROOT", "/work")).resolve()
CLI_DIR = Path(
    os.getenv(
        "CLI_DIR",
        str(APP_ROOT / "CLI" / "RNA_seq_engine_modular_v7_1"),
    )
).resolve()
JOB_ROOT = Path(os.getenv("JOB_ROOT", str(APP_ROOT / "JOBS"))).resolve()

JOB_ROOT.mkdir(parents=True, exist_ok=True)

app = FastAPI(
    title="JCAP RNA-seq API",
    version="0.3.0",
    description="FastAPI wrapper around the modular JCAP RNA-seq engine v7.1.",
)


class JobStatus(BaseModel):
    job_id: str
    status: str
    created_at: str
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    exit_code: Optional[int] = None
    error: Optional[str] = None
    cli_script: str
    output_dir: str
    zip_path: Optional[str] = None


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_filename(name: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in {".", "_", "-"} else "_" for ch in name)
    cleaned = cleaned.strip("._")
    return cleaned or "uploaded_file"


def job_dir(job_id: str) -> Path:
    return JOB_ROOT / job_id


def status_path(job_id: str) -> Path:
    return job_dir(job_id) / "status.json"


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def read_json_if_exists(path: Path) -> Optional[dict]:
    if not path.exists() or not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


def read_status(job_id: str) -> dict:
    path = status_path(job_id)
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    return json.loads(path.read_text())


def write_status(job_id: str, payload: dict) -> None:
    write_json(status_path(job_id), payload)


def resolve_cli_script(script_name: str) -> Path:
    allowed = {"rnaseq_cli.R"}
    if script_name not in allowed:
        raise HTTPException(status_code=400, detail=f"cli_script must be one of: {sorted(allowed)}")

    cli_path = (CLI_DIR / script_name).resolve()

    if not str(cli_path).startswith(str(CLI_DIR)):
        raise HTTPException(status_code=400, detail="Invalid CLI script path")

    if not cli_path.exists():
        raise HTTPException(status_code=500, detail=f"CLI script not found: {cli_path}")

    return cli_path


def validate_optional_choice(value: str, allowed: set[str], field_name: str) -> None:
    if value not in allowed:
        raise HTTPException(status_code=400, detail=f"{field_name} must be one of: {sorted(allowed)}")


def run_rnaseq_job(
    job_id: str,
    cli_script: str,
    counts_path: Path,
    phenotype_path: Path,
    phenotype_column: str,
    species: str,
    contrast: Optional[str],
    reference_group: Optional[str],
    fdr: float,
    logfc_cutoff: float,
    enrich_backend: str,
    gprof_sources: str,
    enrichr_db: str,
    effect_size: float,
    power_test_type: str,
    curve_n_min: int,
    curve_n_max: int,
    classifier_validation: str,
    classifier_model: str,
    zip_outputs: bool,
    drop_library_outliers: bool,
    skip_enrichment: bool,
    skip_annotation: bool,
    skip_classifier: bool,
    skip_plots: bool,
) -> None:
    jd = job_dir(job_id)
    outdir = jd / "results"
    stdout_path = jd / "stdout.log"
    stderr_path = jd / "stderr.log"

    status = read_status(job_id)
    status.update({"status": "running", "started_at": now_iso()})
    write_status(job_id, status)

    cli_path = resolve_cli_script(cli_script)

    cmd = [
        "Rscript",
        str(cli_path),
        "--counts",
        str(counts_path),
        "--phenotype",
        str(phenotype_path),
        "--phenotype-column",
        phenotype_column,
        "--species",
        species,
        "--outdir",
        str(outdir),
        "--fdr",
        str(fdr),
        "--logfc-cutoff",
        str(logfc_cutoff),
        "--effect-size",
        str(effect_size),
        "--power-test-type",
        power_test_type,
        "--curve-n-min",
        str(curve_n_min),
        "--curve-n-max",
        str(curve_n_max),
        "--enrich-backend",
        enrich_backend,
        "--gprof-sources",
        gprof_sources,
        "--enrichr-db",
        enrichr_db,
        "--classifier-validation",
        classifier_validation,
        "--classifier-model",
        classifier_model,
    ]

    if contrast:
        cmd.extend(["--contrast", contrast])
    if reference_group:
        cmd.extend(["--reference-group", reference_group])
    if zip_outputs:
        cmd.append("--zip")
    if drop_library_outliers:
        cmd.append("--drop-library-outliers")
    if skip_enrichment:
        cmd.append("--skip-enrichment")
    if skip_annotation:
        cmd.append("--skip-annotation")
    if skip_classifier:
        cmd.append("--skip-classifier")
    if skip_plots:
        cmd.append("--skip-plots")

    write_json(jd / "command.json", {"cmd": cmd})

    try:
        with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
            proc = subprocess.run(
                cmd,
                cwd=str(cli_path.parent),
                stdout=stdout,
                stderr=stderr,
                text=True,
                check=False,
            )

        cli_run_status = read_json_if_exists(outdir / "run_status.json")
        manifest = read_json_if_exists(outdir / "manifest.json")

        zip_path: Optional[str] = None
        if cli_run_status and cli_run_status.get("zip_path"):
            maybe_zip = Path(str(cli_run_status["zip_path"]))
            if not maybe_zip.is_absolute():
                # v6 writes a path like /work/JOBS/<job>/results.zip when outdir is absolute,
                # but keep this safe for relative paths too.
                maybe_zip = (APP_ROOT / maybe_zip).resolve()
            if maybe_zip.exists():
                zip_path = str(maybe_zip)
        elif (jd / "results.zip").exists():
            zip_path = str(jd / "results.zip")

        status = read_status(job_id)
        status.update(
            {
                "finished_at": now_iso(),
                "exit_code": proc.returncode,
                "status": "completed" if proc.returncode == 0 else "failed",
                "error": None if proc.returncode == 0 else f"Rscript exited with code {proc.returncode}",
                "zip_path": zip_path,
                "cli_run_status": cli_run_status,
                "manifest_available": manifest is not None,
            }
        )
        write_status(job_id, status)

    except Exception as exc:
        status = read_status(job_id)
        status.update(
            {
                "status": "failed",
                "finished_at": now_iso(),
                "exit_code": None,
                "error": str(exc),
            }
        )
        write_status(job_id, status)


@app.get("/health")
def health() -> dict:
    return {
        "ok": True,
        "app_root": str(APP_ROOT),
        "cli_dir": str(CLI_DIR),
        "job_root": str(JOB_ROOT),
        "api_version": "0.3.0",
        "expected_cli": "modular rnaseq_cli.R v7.1",
        "available_cli_scripts": [p.name for p in CLI_DIR.glob("*.R")] if CLI_DIR.exists() else [],
    }


@app.post("/rnaseq/jobs", response_model=JobStatus)
async def create_rnaseq_job(
    background_tasks: BackgroundTasks,
    counts: UploadFile = File(...),
    phenotype: UploadFile = File(...),
    phenotype_column: str = Form(...),
    species: str = Form("hsapiens"),
    cli_script: str = Form("rnaseq_cli.R"),
    contrast: Optional[str] = Form(None),
    reference_group: Optional[str] = Form(None),
    fdr: float = Form(0.05),
    logfc_cutoff: float = Form(1.0),
    enrich_backend: str = Form("auto"),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    effect_size: float = Form(0.8),
    power_test_type: str = Form("ttest"),
    curve_n_min: int = Form(2),
    curve_n_max: int = Form(30),
    classifier_validation: str = Form("cv"),
    classifier_model: str = Form("auto"),
    zip_outputs: bool = Form(True),
    drop_library_outliers: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_annotation: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_plots: bool = Form(False),
) -> JobStatus:
    resolve_cli_script(cli_script)

    validate_optional_choice(species, {"hsapiens", "mmusculus", "drerio", "dmelanogaster"}, "species")
    validate_optional_choice(power_test_type, {"ttest", "anova"}, "power_test_type")
    validate_optional_choice(enrich_backend, {"auto", "gprof", "enrichr"}, "enrich_backend")
    validate_optional_choice(classifier_validation, {"cv", "train_test", "none"}, "classifier_validation")
    validate_optional_choice(classifier_model, {"auto", "rf", "logistic"}, "classifier_model")

    if not 0 < fdr <= 1:
        raise HTTPException(status_code=400, detail="fdr must be > 0 and <= 1")
    if logfc_cutoff < 0:
        raise HTTPException(status_code=400, detail="logfc_cutoff must be >= 0")
    if curve_n_min < 2 or curve_n_max < curve_n_min:
        raise HTTPException(status_code=400, detail="curve_n_min must be >= 2 and curve_n_max must be >= curve_n_min")
    job_id = uuid.uuid4().hex
    jd = job_dir(job_id)
    input_dir = jd / "inputs"
    result_dir = jd / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    counts_path = input_dir / safe_filename(counts.filename or "counts.csv")
    phenotype_path = input_dir / safe_filename(phenotype.filename or "phenotype.csv")

    with counts_path.open("wb") as f:
        shutil.copyfileobj(counts.file, f)

    with phenotype_path.open("wb") as f:
        shutil.copyfileobj(phenotype.file, f)

    status = {
        "job_id": job_id,
        "status": "queued",
        "created_at": now_iso(),
        "started_at": None,
        "finished_at": None,
        "exit_code": None,
        "error": None,
        "cli_script": cli_script,
        "output_dir": str(result_dir),
        "zip_path": None,
        "parameters": {
            "phenotype_column": phenotype_column,
            "species": species,
            "contrast": contrast,
            "reference_group": reference_group,
            "fdr": fdr,
            "logfc_cutoff": logfc_cutoff,
            "enrich_backend": enrich_backend,
            "gprof_sources": gprof_sources,
            "enrichr_db": enrichr_db,
            "classifier_validation": classifier_validation,
            "classifier_model": classifier_model,
            "zip_outputs": zip_outputs,
        },
    }
    write_status(job_id, status)

    background_tasks.add_task(
        run_rnaseq_job,
        job_id,
        cli_script,
        counts_path,
        phenotype_path,
        phenotype_column,
        species,
        contrast,
        reference_group,
        fdr,
        logfc_cutoff,
        enrich_backend,
        gprof_sources,
        enrichr_db,
        effect_size,
        power_test_type,
        curve_n_min,
        curve_n_max,
        classifier_validation,
        classifier_model,
        zip_outputs,
        drop_library_outliers,
        skip_enrichment,
        skip_annotation,
        skip_classifier,
        skip_plots,
    )

    return JobStatus(**status)


@app.get("/rnaseq/jobs/{job_id}", response_model=JobStatus)
def get_rnaseq_job(job_id: str) -> JobStatus:
    return JobStatus(**read_status(job_id))


@app.get("/rnaseq/jobs/{job_id}/run-status")
def get_rnaseq_cli_run_status(job_id: str) -> dict:
    read_status(job_id)
    path = job_dir(job_id) / "results" / "run_status.json"
    payload = read_json_if_exists(path)
    if payload is None:
        raise HTTPException(status_code=404, detail="CLI run_status.json not found yet")
    return payload


@app.get("/rnaseq/jobs/{job_id}/summary")
def get_rnaseq_summary(job_id: str) -> dict:
    read_status(job_id)
    path = job_dir(job_id) / "results" / "run_summary.json"
    payload = read_json_if_exists(path)
    if payload is None:
        raise HTTPException(status_code=404, detail="run_summary.json not found yet")
    return payload


@app.get("/rnaseq/jobs/{job_id}/manifest")
def get_rnaseq_manifest(job_id: str) -> dict:
    read_status(job_id)
    path = job_dir(job_id) / "results" / "manifest.json"
    payload = read_json_if_exists(path)
    if payload is None:
        raise HTTPException(status_code=404, detail="manifest.json not found yet")
    return payload


@app.get("/rnaseq/jobs/{job_id}/files")
def list_rnaseq_job_files(job_id: str) -> dict:
    read_status(job_id)
    jd = job_dir(job_id)
    files = []
    for path in jd.rglob("*"):
        if path.is_file():
            files.append(str(path.relative_to(jd)))
    return {"job_id": job_id, "files": sorted(files)}


@app.get("/rnaseq/jobs/{job_id}/files/{file_path:path}")
def get_rnaseq_job_file(job_id: str, file_path: str) -> FileResponse:
    read_status(job_id)
    jd = job_dir(job_id).resolve()
    target = (jd / file_path).resolve()

    if not str(target).startswith(str(jd)):
        raise HTTPException(status_code=400, detail="Invalid file path")
    if not target.exists() or not target.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {file_path}")

    return FileResponse(target)


@app.get("/rnaseq/jobs/{job_id}/download")
def download_rnaseq_job_zip(job_id: str) -> FileResponse:
    status = read_status(job_id)
    jd = job_dir(job_id)

    # Prefer the CLI's v6 --zip bundle when present.
    zip_from_status = status.get("zip_path")
    if zip_from_status:
        zp = Path(zip_from_status)
        if zp.exists() and zp.is_file():
            return FileResponse(zp, filename=f"{job_id}_rnaseq_results.zip", media_type="application/zip")

    cli_zip = jd / "results.zip"
    if cli_zip.exists() and cli_zip.is_file():
        return FileResponse(cli_zip, filename=f"{job_id}_rnaseq_results.zip", media_type="application/zip")

    # Fallback: build an API-level bundle of the whole job directory.
    zip_path = jd / f"{job_id}_results.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in jd.rglob("*"):
            if path.is_file() and path != zip_path:
                zf.write(path, arcname=str(path.relative_to(jd)))

    return FileResponse(zip_path, filename=f"{job_id}_results.zip", media_type="application/zip")


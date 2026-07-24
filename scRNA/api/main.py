from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import uuid
import zipfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Literal

from fastapi import BackgroundTasks, FastAPI, File, Form, HTTPException, UploadFile, status
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

APP_VERSION = "0.4.0"
# SCRNA_API_V4_LOCAL_PATH_ENDPOINTS_2026_07_03
# SCRNA_API_V3_ROUTE_BUTTONS_FIXED_2026_07_03


class JobState(str, Enum):
    queued = "queued"
    running = "running"
    complete = "complete"
    failed = "failed"
    deleted = "deleted"


class EnrichBackend(str, Enum):
    gprof = "gprof"
    enrichr = "enrichr"


class Species(str, Enum):
    hsapiens = "hsapiens"
    mmusculus = "mmusculus"
    drerio = "drerio"
    dmelanogaster = "dmelanogaster"


class IdType(str, Enum):
    Symbol = "Symbol"
    Ensembl = "Ensembl"
    Entrez = "Entrez"


@dataclass(frozen=True)
class Settings:
    base_dir: Path = Path(__file__).resolve().parents[1]
    max_upload_mb: int = int(os.environ.get("SCRNA_MAX_UPLOAD_MB", "500"))
    cli_timeout_seconds: int = int(os.environ.get("SCRNA_CLI_TIMEOUT_SECONDS", "7200"))
    jobs_dir: Path = Path(os.environ.get("SCRNA_JOBS_DIR", Path(__file__).resolve().parents[1] / "jobs"))
    cli_script: Path = Path(os.environ.get("SCRNA_CLI_SCRIPT", Path(__file__).resolve().parents[1] / "cli" / "run_scrna.R"))


settings = Settings()
settings.jobs_dir.mkdir(parents=True, exist_ok=True)

app = FastAPI(
    title="JCAP scRNA-seq API",
    version=APP_VERSION,
    description="Production-oriented FastAPI wrapper for the scRNA-seq Seurat CLI module.",
)


class HealthResponse(BaseModel):
    status: Literal["ok", "degraded"]
    version: str
    cli_script: str
    cli_exists: bool
    jobs_dir: str
    jobs_dir_writable: bool


class SubmitResponse(BaseModel):
    job_id: str
    status: JobState
    status_url: str
    download_url: str | None = None
    tables_url: str | None = None


class JobStatusResponse(BaseModel):
    job_id: str
    status: JobState
    created_at: float | None = None
    started_at: float | None = None
    finished_at: float | None = None
    returncode: int | None = None
    error: str | None = None
    summary: dict[str, Any] = Field(default_factory=dict)
    run_status: dict[str, Any] = Field(default_factory=dict)
    has_bundle: bool = False
    download_url: str | None = None
    tables_url: str | None = None


class ValidationResponse(BaseModel):
    job_id: str
    counts_file: str
    meta_file: str
    counts_columns_preview: list[str]
    metadata_columns_preview: list[str]
    warnings: list[str] = Field(default_factory=list)
    notes: list[str]


SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
ALLOWED_TABLE_SUFFIXES = {".csv"}
ALLOWED_UPLOAD_SUFFIXES = {".csv", ".tsv", ".txt", ".rds", ".h5", ".h5ad", ".mtx", ".gz"}


def _now() -> float:
    return time.time()


def _job_id() -> str:
    return str(uuid.uuid4())


def _parse_job_id(job_id: str) -> str:
    try:
        return str(uuid.UUID(job_id))
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="Invalid job_id") from exc


def _job_dir(job_id: str) -> Path:
    return settings.jobs_dir / _parse_job_id(job_id)


def _safe_name(name: str, fallback: str = "upload.dat") -> str:
    raw = Path(name or fallback).name.strip() or fallback
    cleaned = SAFE_FILENAME_RE.sub("_", raw)
    return cleaned[:180]


def _ensure_allowed_upload_name(filename: str) -> None:
    suffixes = [s.lower() for s in Path(filename).suffixes]
    if not suffixes or not any(s in ALLOWED_UPLOAD_SUFFIXES for s in suffixes):
        raise HTTPException(status_code=415, detail=f"Unsupported upload type for {filename}")


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8") as tmp:
        json.dump(payload, tmp, indent=2, sort_keys=True)
        tmp.write("\n")
        tmp_path = Path(tmp.name)
    tmp_path.replace(path)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _write_state(job_dir: Path, **updates: Any) -> dict[str, Any]:
    state_path = job_dir / "job_state.json"
    state = _read_json(state_path)
    state.update(updates)
    _atomic_write_json(state_path, state)
    return state


def _save_upload(upload: UploadFile, dest: Path) -> int:
    filename = _safe_name(upload.filename or dest.name)
    _ensure_allowed_upload_name(filename)
    max_bytes = settings.max_upload_mb * 1024 * 1024
    size = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with dest.open("wb") as fh:
            while chunk := upload.file.read(1024 * 1024):
                size += len(chunk)
                if size > max_bytes:
                    raise HTTPException(status_code=413, detail=f"File exceeds {settings.max_upload_mb} MB limit")
                fh.write(chunk)
    finally:
        upload.file.close()
    if size == 0:
        raise HTTPException(status_code=400, detail=f"Empty upload: {filename}")
    return size


def _csv_header(path: Path, delimiter: str = ",") -> list[str]:
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as fh:
            sample = fh.readline()
            if not sample:
                return []
            return [col.strip() for col in next(csv.reader([sample], delimiter=delimiter))]
    except Exception:
        return []


def _zip_dir(src_dir: Path, zip_path: Path) -> None:
    """Create the bundle atomically.

    Write to a temporary filename first so the status endpoint never reports
    has_bundle=True while the ZIP is still being created.
    """
    tmp_zip = zip_path.with_suffix(zip_path.suffix + ".tmp")
    tmp_zip.unlink(missing_ok=True)

    try:
        with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            for path in src_dir.rglob("*"):
                if path.is_file() and path not in {zip_path, tmp_zip}:
                    zf.write(path, path.relative_to(src_dir))

        # Atomic replacement on the same filesystem.
        tmp_zip.replace(zip_path)
    except Exception:
        tmp_zip.unlink(missing_ok=True)
        raise


def _build_cmd(
    counts_path: Path,
    meta_path: Path,
    out_dir: Path,
    species: Species,
    id_type: IdType,
    enrich_backend: EnrichBackend,
    gprof_sources: str,
    enrichr_db: str,
    top_n_features: int,
    max_pcs: int,
    skip_pca_umap: bool,
    skip_enrichment: bool,
    skip_classifier: bool,
    skip_power: bool,
    condition_col: str | None = None,
    celltype_col: str | None = None,
    annotation_col: str | None = None,
    sample_col: str | None = None,
    condition_a: str | None = None,
    condition_b: str | None = None,
    target_celltype: str | None = None,
    de_mode: str | None = None,
    de_scope: str | None = None,
    classifier_validation: str | None = None,
    classifier_model: str | None = None,
    classifier_selector: str | None = None,
    classifier_feature_source: str | None = None,
) -> list[str]:
    cmd = [
        "Rscript", str(settings.cli_script),
        "--counts", str(counts_path),
        "--metadata", str(meta_path),
        "--outdir", str(out_dir),
        "--species", species.value,
        "--id-type", id_type.value,
        "--enrich-backend", enrich_backend.value,
        "--gprof-sources", gprof_sources,
        "--enrichr-db", enrichr_db,
        "--top-n-features", str(top_n_features),
        "--max-pcs", str(max_pcs),
    ]

    optional_flags = {
        "--condition-col": condition_col,
        "--celltype-col": celltype_col,
        "--annotation-col": annotation_col,
        "--sample-col": sample_col,
        "--condition-a": condition_a,
        "--condition-b": condition_b,
        "--target-celltype": target_celltype,
        "--de-mode": de_mode,
        "--de-scope": de_scope,
        "--classifier-validation": classifier_validation,
        "--classifier-model": classifier_model,
        "--classifier-selector": classifier_selector,
        "--classifier-feature-source": classifier_feature_source,
    }
    for flag, val in optional_flags.items():
        if val not in (None, ""):
            cmd.extend([flag, str(val)])

    if skip_pca_umap:
        cmd.append("--skip-pca-umap")
    if skip_enrichment:
        cmd.append("--skip-enrichment")
    if skip_classifier:
        cmd.append("--skip-classifier")
    if skip_power:
        cmd.append("--skip-power")
    return cmd


def _run_job(job_id: str, cmd: list[str]) -> None:
    job_dir = _job_dir(job_id)
    out_dir = job_dir / "results"
    run_log = job_dir / "run.log"
    _write_state(job_dir, status=JobState.running.value, started_at=_now())

    try:
        result = subprocess.run(
            cmd,
            cwd=str(job_dir),
            capture_output=True,
            text=True,
            timeout=settings.cli_timeout_seconds,
            check=False,
        )
        run_log.write_text(
            "COMMAND:\n" + " ".join(cmd) +
            "\n\nSTDOUT:\n" + result.stdout +
            "\n\nSTDERR:\n" + result.stderr,
            encoding="utf-8",
        )
        if result.returncode != 0:
            _write_state(
                job_dir,
                status=JobState.failed.value,
                finished_at=_now(),
                returncode=result.returncode,
                error=result.stderr[-4000:] or "CLI failed without stderr",
            )
            return

        bundle_path = job_dir / f"scrna_results_{job_id}.zip"

        # The R workflow is finished, but packaging may still take time for a
        # large Seurat object. Keep the top-level job running until the ZIP is
        # fully written and atomically renamed into place.
        _write_state(
            job_dir,
            status=JobState.running.value,
            summary=_read_json(out_dir / "run_summary.json"),
            run_status=_read_json(out_dir / "run_status.json"),
        )
        _zip_dir(out_dir, bundle_path)

        _write_state(
            job_dir,
            status=JobState.complete.value,
            finished_at=_now(),
            returncode=0,
            summary=_read_json(out_dir / "run_summary.json"),
            run_status=_read_json(out_dir / "run_status.json"),
            download_url=f"/jobs/{job_id}/download",
            tables_url=f"/jobs/{job_id}/tables",
        )
    except subprocess.TimeoutExpired:
        _write_state(
            job_dir,
            status=JobState.failed.value,
            finished_at=_now(),
            error=f"CLI timed out after {settings.cli_timeout_seconds} seconds",
        )
    except Exception as exc:
        _write_state(job_dir, status=JobState.failed.value, finished_at=_now(), error=str(exc))


@app.get("/", response_model=HealthResponse)
def root() -> HealthResponse:
    return health()


@app.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    cli_exists = settings.cli_script.exists()
    jobs_dir_writable = os.access(settings.jobs_dir, os.W_OK)
    overall = "ok" if cli_exists and jobs_dir_writable else "degraded"
    return HealthResponse(
        status=overall,
        version=APP_VERSION,
        cli_script=str(settings.cli_script),
        cli_exists=cli_exists,
        jobs_dir=str(settings.jobs_dir),
        jobs_dir_writable=jobs_dir_writable,
    )


@app.post("/validate", response_model=ValidationResponse)
def validate_inputs(
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
) -> ValidationResponse:
    job_id = _job_id()
    job_dir = _job_dir(job_id)
    input_dir = job_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    counts_name = "counts_" + _safe_name(counts_file.filename or "counts.csv")
    meta_name = "metadata_" + _safe_name(meta_file.filename or "metadata.csv")
    counts_path = input_dir / counts_name
    meta_path = input_dir / meta_name
    _save_upload(counts_file, counts_path)
    _save_upload(meta_file, meta_path)

    counts_header = _csv_header(counts_path)
    meta_header = _csv_header(meta_path)
    warnings: list[str] = []
    if not counts_header:
        warnings.append("Could not read a CSV-style header from counts_file.")
    if not meta_header:
        warnings.append("Could not read a CSV-style header from meta_file.")
    for required in ("stim", "cell_type"):
        if meta_header and required not in meta_header:
            warnings.append(f"Metadata header does not contain required column: {required}")

    _write_state(job_dir, status=JobState.complete.value, created_at=_now(), validation_only=True)
    return ValidationResponse(
        job_id=job_id,
        counts_file=counts_path.name,
        meta_file=meta_path.name,
        counts_columns_preview=counts_header[:10],
        metadata_columns_preview=meta_header[:20],
        warnings=warnings,
        notes=[
            "Counts CSV should have genes as rows and cells as columns.",
            "Metadata CSV should have cells as rows; stim and cell_type are required for the full pipeline.",
        ],
    )


@app.post("/jobs", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_job(
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    if not settings.cli_script.exists():
        raise HTTPException(status_code=503, detail="scRNA CLI script is not available")

    job_id = _job_id()
    job_dir = _job_dir(job_id)
    input_dir = job_dir / "inputs"
    out_dir = job_dir / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = input_dir / ("counts_" + _safe_name(counts_file.filename or "counts.csv"))
    meta_path = input_dir / ("metadata_" + _safe_name(meta_file.filename or "metadata.csv"))
    counts_size = _save_upload(counts_file, counts_path)
    meta_size = _save_upload(meta_file, meta_path)

    cmd = _build_cmd(
        counts_path=counts_path,
        meta_path=meta_path,
        out_dir=out_dir,
        species=species,
        id_type=id_type,
        enrich_backend=enrich_backend,
        gprof_sources=gprof_sources,
        enrichr_db=enrichr_db,
        top_n_features=top_n_features,
        max_pcs=max_pcs,
        skip_pca_umap=skip_pca_umap,
        skip_enrichment=skip_enrichment,
        skip_classifier=skip_classifier,
        skip_power=skip_power,
    )
    _write_state(
        job_dir,
        job_id=job_id,
        status=JobState.queued.value,
        created_at=_now(),
        counts_file=counts_path.name,
        meta_file=meta_path.name,
        counts_size_bytes=counts_size,
        meta_size_bytes=meta_size,
        params={
            "species": species.value,
            "id_type": id_type.value,
            "enrich_backend": enrich_backend.value,
            "gprof_sources": gprof_sources,
            "enrichr_db": enrichr_db,
            "top_n_features": top_n_features,
            "max_pcs": max_pcs,
            "skip_pca_umap": skip_pca_umap,
            "skip_enrichment": skip_enrichment,
            "skip_classifier": skip_classifier,
            "skip_power": skip_power,
        },
    )
    background_tasks.add_task(_run_job, job_id, cmd)
    return SubmitResponse(job_id=job_id, status=JobState.queued, status_url=f"/jobs/{job_id}")


@app.post("/run", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def backward_compatible_run(
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return submit_scrna_job(
        background_tasks=background_tasks,
        counts_file=counts_file,
        meta_file=meta_file,
        species=species,
        id_type=id_type,
        enrich_backend=enrich_backend,
        gprof_sources=gprof_sources,
        enrichr_db=enrichr_db,
        top_n_features=top_n_features,
        max_pcs=max_pcs,
        skip_pca_umap=skip_pca_umap,
        skip_enrichment=skip_enrichment,
        skip_classifier=skip_classifier,
        skip_power=skip_power,
    )




SCRNA_ROUTE_DEFAULTS = {
    "cell": {"de_mode": "cell", "de_scope": "top_celltype", "classifier_validation": "none"},
    "pseudobulk": {"de_mode": "pseudobulk", "de_scope": "top_celltype", "classifier_validation": "none"},
    "global_pseudobulk": {"de_mode": "pseudobulk", "de_scope": "global", "classifier_validation": "none"},
    "target_celltype_pseudobulk": {"de_mode": "pseudobulk", "de_scope": "target_celltype", "classifier_validation": "none"},
    "classifier": {"de_mode": "cell", "de_scope": "top_celltype", "classifier_validation": "sample_cv"},
}



def _resolve_container_input_path(path_text: str, label: str) -> Path:
    """Resolve a server/container-side input path for large local datasets.

    The Dash UI may send either an absolute container path such as
    /work/sample_data/file.csv, or a project-relative path such as
    modules/scRNA/sample_data/file.csv. Inside the scRNA API container the
    project module is mounted at /work, so modules/scRNA/... maps to /work/...
    """
    if not path_text or not str(path_text).strip():
        raise HTTPException(status_code=400, detail=f"{label} is required")
    raw = str(path_text).strip().strip('"').strip("'")
    candidates: list[Path] = []
    p = Path(raw)
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.append(settings.base_dir / raw)
        candidates.append(Path("/work") / raw)
        marker = "modules/scRNA/"
        norm = raw.replace("\\", "/")
        if marker in norm:
            candidates.append(Path("/work") / norm.split(marker, 1)[1])
    # also map absolute Windows/WSL project paths containing modules/scRNA
    norm_abs = raw.replace("\\", "/")
    marker = "/modules/scRNA/"
    if marker in norm_abs:
        candidates.append(Path("/work") / norm_abs.split(marker, 1)[1])

    seen = set()
    unique = []
    for c in candidates:
        s = str(c)
        if s not in seen:
            unique.append(c)
            seen.add(s)
    for c in unique:
        if c.exists() and c.is_file():
            return c
    raise HTTPException(status_code=400, detail=f"Could not resolve {label}: {raw}. Tried: " + "; ".join(map(str, unique[:6])))


def _submit_route_job_from_paths(
    background_tasks: BackgroundTasks,
    route: str,
    counts_path: Path,
    meta_path: Path,
    species: Species,
    id_type: IdType,
    condition_col: str,
    celltype_col: str,
    annotation_col: str,
    sample_col: str,
    condition_a: str,
    condition_b: str,
    target_celltype: str,
    enrich_backend: EnrichBackend,
    gprof_sources: str,
    enrichr_db: str,
    top_n_features: int,
    max_pcs: int,
    classifier_validation: str,
    classifier_model: str,
    classifier_selector: str,
    classifier_feature_source: str,
    skip_pca_umap: bool,
    skip_enrichment: bool,
    skip_classifier: bool,
    skip_power: bool,
) -> SubmitResponse:
    route = route.strip()
    if route not in SCRNA_ROUTE_DEFAULTS:
        raise HTTPException(status_code=404, detail=f"Unknown scRNA route: {route}")
    if not settings.cli_script.exists():
        raise HTTPException(status_code=503, detail="scRNA CLI script is not available")

    defaults = SCRNA_ROUTE_DEFAULTS[route]
    effective_classifier_validation = classifier_validation or defaults["classifier_validation"]
    effective_skip_classifier = skip_classifier or (effective_classifier_validation == "none" and route != "classifier")

    job_id = _job_id()
    job_dir = _job_dir(job_id)
    out_dir = job_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = _build_cmd(
        counts_path=counts_path,
        meta_path=meta_path,
        out_dir=out_dir,
        species=species,
        id_type=id_type,
        enrich_backend=enrich_backend,
        gprof_sources=gprof_sources,
        enrichr_db=enrichr_db,
        top_n_features=top_n_features,
        max_pcs=max_pcs,
        skip_pca_umap=skip_pca_umap,
        skip_enrichment=skip_enrichment,
        skip_classifier=effective_skip_classifier,
        skip_power=skip_power,
        condition_col=condition_col or None,
        celltype_col=celltype_col or None,
        annotation_col=annotation_col or None,
        sample_col=sample_col or None,
        condition_a=condition_a or None,
        condition_b=condition_b or None,
        target_celltype=target_celltype or None,
        de_mode=defaults["de_mode"],
        de_scope=defaults["de_scope"],
        classifier_validation=effective_classifier_validation,
        classifier_model=classifier_model or None,
        classifier_selector=classifier_selector or None,
        classifier_feature_source=classifier_feature_source or None,
    )

    _write_state(
        job_dir,
        job_id=job_id,
        route=route,
        status=JobState.queued.value,
        created_at=_now(),
        input_mode="server_path",
        counts_path=str(counts_path),
        meta_path=str(meta_path),
        counts_size_bytes=counts_path.stat().st_size,
        meta_size_bytes=meta_path.stat().st_size,
        command=cmd,
        params={
            "route": route,
            "species": species.value,
            "id_type": id_type.value,
            "condition_col": condition_col,
            "celltype_col": celltype_col,
            "annotation_col": annotation_col,
            "sample_col": sample_col,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "target_celltype": target_celltype,
            "de_mode": defaults["de_mode"],
            "de_scope": defaults["de_scope"],
            "classifier_validation": effective_classifier_validation,
        },
    )
    background_tasks.add_task(_run_job, job_id, cmd)
    return SubmitResponse(job_id=job_id, status=JobState.queued, status_url=f"/jobs/{job_id}", download_url=None, tables_url=None)

def _submit_route_job(
    background_tasks: BackgroundTasks,
    route: str,
    counts_file: UploadFile,
    meta_file: UploadFile,
    species: Species,
    id_type: IdType,
    condition_col: str,
    celltype_col: str,
    annotation_col: str,
    sample_col: str,
    condition_a: str,
    condition_b: str,
    target_celltype: str,
    enrich_backend: EnrichBackend,
    gprof_sources: str,
    enrichr_db: str,
    top_n_features: int,
    max_pcs: int,
    classifier_validation: str,
    classifier_model: str,
    classifier_selector: str,
    classifier_feature_source: str,
    skip_pca_umap: bool,
    skip_enrichment: bool,
    skip_classifier: bool,
    skip_power: bool,
) -> SubmitResponse:
    route = route.strip()
    if route not in SCRNA_ROUTE_DEFAULTS:
        raise HTTPException(status_code=404, detail=f"Unknown scRNA route: {route}")
    if not settings.cli_script.exists():
        raise HTTPException(status_code=503, detail="scRNA CLI script is not available")

    defaults = SCRNA_ROUTE_DEFAULTS[route]
    effective_classifier_validation = classifier_validation or defaults["classifier_validation"]
    effective_skip_classifier = skip_classifier or (effective_classifier_validation == "none" and route != "classifier")

    job_id = _job_id()
    job_dir = _job_dir(job_id)
    input_dir = job_dir / "inputs"
    out_dir = job_dir / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = input_dir / ("counts_" + _safe_name(counts_file.filename or "counts.csv"))
    meta_path = input_dir / ("metadata_" + _safe_name(meta_file.filename or "metadata.csv"))
    counts_size = _save_upload(counts_file, counts_path)
    meta_size = _save_upload(meta_file, meta_path)

    cmd = _build_cmd(
        counts_path=counts_path,
        meta_path=meta_path,
        out_dir=out_dir,
        species=species,
        id_type=id_type,
        enrich_backend=enrich_backend,
        gprof_sources=gprof_sources,
        enrichr_db=enrichr_db,
        top_n_features=top_n_features,
        max_pcs=max_pcs,
        skip_pca_umap=skip_pca_umap,
        skip_enrichment=skip_enrichment,
        skip_classifier=effective_skip_classifier,
        skip_power=skip_power,
        condition_col=condition_col or None,
        celltype_col=celltype_col or None,
        annotation_col=annotation_col or None,
        sample_col=sample_col or None,
        condition_a=condition_a or None,
        condition_b=condition_b or None,
        target_celltype=target_celltype or None,
        de_mode=defaults["de_mode"],
        de_scope=defaults["de_scope"],
        classifier_validation=effective_classifier_validation,
        classifier_model=classifier_model or None,
        classifier_selector=classifier_selector or None,
        classifier_feature_source=classifier_feature_source or None,
    )

    _write_state(
        job_dir,
        job_id=job_id,
        route=route,
        status=JobState.queued.value,
        created_at=_now(),
        counts_file=counts_path.name,
        meta_file=meta_path.name,
        counts_size_bytes=counts_size,
        meta_size_bytes=meta_size,
        command=cmd,
        params={
            "route": route,
            "species": species.value,
            "id_type": id_type.value,
            "condition_col": condition_col,
            "celltype_col": celltype_col,
            "annotation_col": annotation_col,
            "sample_col": sample_col,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "target_celltype": target_celltype,
            "de_mode": defaults["de_mode"],
            "de_scope": defaults["de_scope"],
            "classifier_validation": effective_classifier_validation,
        },
    )
    background_tasks.add_task(_run_job, job_id, cmd)
    return SubmitResponse(job_id=job_id, status=JobState.queued, status_url=f"/jobs/{job_id}", download_url=None, tables_url=None)


def _route_endpoint(
    background_tasks: BackgroundTasks,
    route: str,
    counts_file: UploadFile,
    meta_file: UploadFile,
    species: Species,
    id_type: IdType,
    condition_col: str,
    celltype_col: str,
    annotation_col: str,
    sample_col: str,
    condition_a: str,
    condition_b: str,
    target_celltype: str,
    enrich_backend: EnrichBackend,
    gprof_sources: str,
    enrichr_db: str,
    top_n_features: int,
    max_pcs: int,
    classifier_validation: str,
    classifier_model: str,
    classifier_selector: str,
    classifier_feature_source: str,
    skip_pca_umap: bool,
    skip_enrichment: bool,
    skip_classifier: bool,
    skip_power: bool,
) -> SubmitResponse:
    return _submit_route_job(
        background_tasks, route, counts_file, meta_file, species, id_type,
        condition_col, celltype_col, annotation_col, sample_col,
        condition_a, condition_b, target_celltype,
        enrich_backend, gprof_sources, enrichr_db,
        top_n_features, max_pcs,
        classifier_validation, classifier_model, classifier_selector, classifier_feature_source,
        skip_pca_umap, skip_enrichment, skip_classifier, skip_power,
    )


@app.post("/scrna/{route}", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_route_prefixed(
    route: str,
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    condition_col: str = Form("stim"),
    celltype_col: str = Form("cell_type"),
    annotation_col: str = Form("seurat_annotations"),
    sample_col: str = Form("orig.ident"),
    condition_a: str = Form(""),
    condition_b: str = Form(""),
    target_celltype: str = Form(""),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    classifier_validation: str = Form(""),
    classifier_model: str = Form("auto"),
    classifier_selector: str = Form("auto"),
    classifier_feature_source: str = Form("condition_only"),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return _route_endpoint(background_tasks, route, counts_file, meta_file, species, id_type, condition_col, celltype_col, annotation_col, sample_col, condition_a, condition_b, target_celltype, enrich_backend, gprof_sources, enrichr_db, top_n_features, max_pcs, classifier_validation, classifier_model, classifier_selector, classifier_feature_source, skip_pca_umap, skip_enrichment, skip_classifier, skip_power)


@app.post("/scrna/{route}/run", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_route_prefixed_run(
    route: str,
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    condition_col: str = Form("stim"),
    celltype_col: str = Form("cell_type"),
    annotation_col: str = Form("seurat_annotations"),
    sample_col: str = Form("orig.ident"),
    condition_a: str = Form(""),
    condition_b: str = Form(""),
    target_celltype: str = Form(""),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    classifier_validation: str = Form(""),
    classifier_model: str = Form("auto"),
    classifier_selector: str = Form("auto"),
    classifier_feature_source: str = Form("condition_only"),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return _route_endpoint(background_tasks, route, counts_file, meta_file, species, id_type, condition_col, celltype_col, annotation_col, sample_col, condition_a, condition_b, target_celltype, enrich_backend, gprof_sources, enrichr_db, top_n_features, max_pcs, classifier_validation, classifier_model, classifier_selector, classifier_feature_source, skip_pca_umap, skip_enrichment, skip_classifier, skip_power)


@app.post("/scrna/routes/{route}", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_route_routes(
    route: str,
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    condition_col: str = Form("stim"),
    celltype_col: str = Form("cell_type"),
    annotation_col: str = Form("seurat_annotations"),
    sample_col: str = Form("orig.ident"),
    condition_a: str = Form(""),
    condition_b: str = Form(""),
    target_celltype: str = Form(""),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    classifier_validation: str = Form(""),
    classifier_model: str = Form("auto"),
    classifier_selector: str = Form("auto"),
    classifier_feature_source: str = Form("condition_only"),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return _route_endpoint(background_tasks, route, counts_file, meta_file, species, id_type, condition_col, celltype_col, annotation_col, sample_col, condition_a, condition_b, target_celltype, enrich_backend, gprof_sources, enrichr_db, top_n_features, max_pcs, classifier_validation, classifier_model, classifier_selector, classifier_feature_source, skip_pca_umap, skip_enrichment, skip_classifier, skip_power)



def _route_path_endpoint(
    background_tasks: BackgroundTasks,
    route: str,
    counts_path: str,
    meta_path: str,
    species: Species,
    id_type: IdType,
    condition_col: str,
    celltype_col: str,
    annotation_col: str,
    sample_col: str,
    condition_a: str,
    condition_b: str,
    target_celltype: str,
    enrich_backend: EnrichBackend,
    gprof_sources: str,
    enrichr_db: str,
    top_n_features: int,
    max_pcs: int,
    classifier_validation: str,
    classifier_model: str,
    classifier_selector: str,
    classifier_feature_source: str,
    skip_pca_umap: bool,
    skip_enrichment: bool,
    skip_classifier: bool,
    skip_power: bool,
) -> SubmitResponse:
    cp = _resolve_container_input_path(counts_path, "counts_path")
    mp = _resolve_container_input_path(meta_path, "meta_path")
    return _submit_route_job_from_paths(
        background_tasks, route, cp, mp, species, id_type,
        condition_col, celltype_col, annotation_col, sample_col,
        condition_a, condition_b, target_celltype,
        enrich_backend, gprof_sources, enrichr_db,
        top_n_features, max_pcs,
        classifier_validation, classifier_model, classifier_selector, classifier_feature_source,
        skip_pca_umap, skip_enrichment, skip_classifier, skip_power,
    )


@app.post("/scrna/{route}/local", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_route_local(
    route: str,
    background_tasks: BackgroundTasks,
    counts_path: str = Form(...),
    meta_path: str = Form(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    condition_col: str = Form("stim"),
    celltype_col: str = Form("cell_type"),
    annotation_col: str = Form("seurat_annotations"),
    sample_col: str = Form("orig.ident"),
    condition_a: str = Form(""),
    condition_b: str = Form(""),
    target_celltype: str = Form(""),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    classifier_validation: str = Form(""),
    classifier_model: str = Form("auto"),
    classifier_selector: str = Form("auto"),
    classifier_feature_source: str = Form("condition_only"),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return _route_path_endpoint(background_tasks, route, counts_path, meta_path, species, id_type, condition_col, celltype_col, annotation_col, sample_col, condition_a, condition_b, target_celltype, enrich_backend, gprof_sources, enrichr_db, top_n_features, max_pcs, classifier_validation, classifier_model, classifier_selector, classifier_feature_source, skip_pca_umap, skip_enrichment, skip_classifier, skip_power)


@app.post("/scrna/routes/{route}/local", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_route_routes_local(
    route: str,
    background_tasks: BackgroundTasks,
    counts_path: str = Form(...),
    meta_path: str = Form(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    condition_col: str = Form("stim"),
    celltype_col: str = Form("cell_type"),
    annotation_col: str = Form("seurat_annotations"),
    sample_col: str = Form("orig.ident"),
    condition_a: str = Form(""),
    condition_b: str = Form(""),
    target_celltype: str = Form(""),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    classifier_validation: str = Form(""),
    classifier_model: str = Form("auto"),
    classifier_selector: str = Form("auto"),
    classifier_feature_source: str = Form("condition_only"),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return _route_path_endpoint(background_tasks, route, counts_path, meta_path, species, id_type, condition_col, celltype_col, annotation_col, sample_col, condition_a, condition_b, target_celltype, enrich_backend, gprof_sources, enrichr_db, top_n_features, max_pcs, classifier_validation, classifier_model, classifier_selector, classifier_feature_source, skip_pca_umap, skip_enrichment, skip_classifier, skip_power)

@app.get("/scrna/jobs/{job_id}", response_model=JobStatusResponse)
def scrna_job_status_alias(job_id: str) -> JobStatusResponse:
    return job_status(job_id)


@app.get("/scrna/jobs/{job_id}/status", response_model=JobStatusResponse)
def scrna_job_status_alias2(job_id: str) -> JobStatusResponse:
    return job_status(job_id)


@app.get("/scrna/jobs/{job_id}/download")
def scrna_download_alias(job_id: str) -> FileResponse:
    return download_bundle(job_id)

@app.get("/jobs/{job_id}", response_model=JobStatusResponse)
def job_status(job_id: str) -> JobStatusResponse:
    job_dir = _job_dir(job_id)
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    state = _read_json(job_dir / "job_state.json")
    out_dir = job_dir / "results"
    has_bundle = any(job_dir.glob("scrna_results_*.zip"))
    return JobStatusResponse(
        job_id=_parse_job_id(job_id),
        status=JobState(state.get("status", JobState.failed.value)),
        created_at=state.get("created_at"),
        started_at=state.get("started_at"),
        finished_at=state.get("finished_at"),
        returncode=state.get("returncode"),
        error=state.get("error"),
        summary=state.get("summary") or _read_json(out_dir / "run_summary.json"),
        run_status=state.get("run_status") or _read_json(out_dir / "run_status.json"),
        has_bundle=has_bundle,
        download_url=f"/jobs/{job_id}/download" if has_bundle else None,
        tables_url=f"/jobs/{job_id}/tables" if (out_dir / "tables").exists() else None,
    )


@app.get("/jobs/{job_id}/tables")
def list_tables(job_id: str) -> dict[str, Any]:
    tables_dir = _job_dir(job_id) / "results" / "tables"
    if not tables_dir.exists():
        raise HTTPException(status_code=404, detail="Tables not found")
    return {"job_id": _parse_job_id(job_id), "tables": sorted(p.name for p in tables_dir.glob("*.csv"))}


@app.get("/jobs/{job_id}/tables/{filename}")
def download_table(job_id: str, filename: str) -> FileResponse:
    clean = _safe_name(filename)
    path = _job_dir(job_id) / "results" / "tables" / clean
    if path.suffix.lower() not in ALLOWED_TABLE_SUFFIXES or not path.exists() or not path.is_file():
        raise HTTPException(status_code=404, detail="Table not found")
    return FileResponse(path, media_type="text/csv", filename=path.name)


@app.get("/jobs/{job_id}/download")
def download_bundle(job_id: str) -> FileResponse:
    job_dir = _job_dir(job_id)
    bundles = sorted(job_dir.glob("scrna_results_*.zip"))
    if not bundles:
        raise HTTPException(status_code=404, detail="Bundle not found")
    return FileResponse(bundles[0], media_type="application/zip", filename=bundles[0].name)


@app.delete("/jobs/{job_id}")
def delete_job(job_id: str) -> dict[str, Any]:
    job_dir = _job_dir(job_id)
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    shutil.rmtree(job_dir)
    return {"job_id": _parse_job_id(job_id), "deleted": True}

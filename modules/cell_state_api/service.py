from pathlib import Path
import json
import pandas as pd

from models import (
    CellStateRequest,
    CellStateResponse,
    StateCount,
    GroupStateCount,
)


def _overall_summary(df, state_col):
    counts = df[state_col].fillna("Unknown").astype(str).value_counts()
    total = len(df)

    return [
        StateCount(
            state=state,
            n_cells=int(n),
            fraction=float(n / total) if total else 0.0,
        )
        for state, n in counts.items()
    ]


def _group_summary(df, group_col, state_col):
    out = []

    for group, gdf in df.groupby(group_col, dropna=False):
        counts = gdf[state_col].fillna("Unknown").astype(str).value_counts()
        total = len(gdf)

        for state, n in counts.items():
            out.append(
                GroupStateCount(
                    group=str(group),
                    state=str(state),
                    n_cells=int(n),
                    fraction_within_group=float(n / total) if total else 0.0,
                )
            )

    return out


def analyze_cell_states(req: CellStateRequest) -> CellStateResponse:
    metadata_path = Path(req.metadata_path)

    if not metadata_path.exists():
        return CellStateResponse(
            assay=req.assay,
            status="error",
            annotation_mode=req.annotation_mode,
            warnings=[f"Metadata file does not exist: {metadata_path}"],
        )

    df = pd.read_csv(metadata_path)

    if req.cell_id_column not in df.columns:
        return CellStateResponse(
            assay=req.assay,
            status="error",
            annotation_mode=req.annotation_mode,
            warnings=[
                f"Missing cell ID column: {req.cell_id_column}"
            ],
        )

    # ---------- PROVIDED ANNOTATION ----------
    if req.annotation_mode == "provided":

        if not req.state_column or req.state_column not in df.columns:
            return CellStateResponse(
                assay=req.assay,
                status="needs_annotation",
                annotation_mode=req.annotation_mode,
                n_cells=len(df),
                warnings=[
                    "No valid state_column was provided.",
                    "Pipeline should continue without state-aware analysis "
                    "or supply a reference/marker annotation method."
                ],
            )

        state_col = req.state_column

    # ---------- FUTURE MODES ----------
    elif req.annotation_mode in {"reference_mapping", "marker_rules"}:

        return CellStateResponse(
            assay=req.assay,
            status="needs_annotation",
            annotation_mode=req.annotation_mode,
            n_cells=len(df),
            warnings=[
                f"{req.annotation_mode} is declared but not yet implemented."
            ],
        )

    else:
        return CellStateResponse(
            assay=req.assay,
            status="error",
            annotation_mode=req.annotation_mode,
            warnings=["Unsupported annotation mode."],
        )

    overall = _overall_summary(df, state_col)

    # Assay-aware grouping
    if req.assay == "scrna":
        group_col = req.condition_column

        if not group_col or group_col not in df.columns:
            group_summary = []
            warnings = [
                "No valid condition column supplied; "
                "overall state frequencies were computed only."
            ]
        else:
            group_summary = _group_summary(df, group_col, state_col)
            warnings = []

    else:  # perturbseq
        group_col = req.perturbation_column

        if not group_col or group_col not in df.columns:
            group_summary = []
            warnings = [
                "No valid perturbation column supplied; "
                "overall state frequencies were computed only."
            ]
        else:
            group_summary = _group_summary(df, group_col, state_col)
            warnings = []

    outdir = metadata_path.parent / "cell_state_output"
    outdir.mkdir(parents=True, exist_ok=True)

    annotated_path = outdir / "metadata_with_cell_state.csv"
    df.to_csv(annotated_path, index=False)

    handoff = {
        "analysis_type": "cell_state",
        "assay": req.assay,
        "annotation_mode": req.annotation_mode,
        "state_column": state_col,
        "group_column": group_col,
        "n_cells": len(df),
        "states": sorted(
            df[state_col].dropna().astype(str).unique().tolist()
        ),
        "overall_state_summary": [x.model_dump() for x in overall],
        "group_state_summary": [x.model_dump() for x in group_summary],
        "warnings": warnings,
    }

    handoff_path = outdir / "cell_state_handoff.json"
    handoff_path.write_text(json.dumps(handoff, indent=2))

    return CellStateResponse(
        assay=req.assay,
        status="ok",
        annotation_mode=req.annotation_mode,
        annotation_source=f"metadata:{state_col}",
        n_cells=len(df),
        states=handoff["states"],
        overall_state_summary=overall,
        group_state_summary=group_summary,
        annotated_metadata_path=str(annotated_path),
        handoff_path=str(handoff_path),
        warnings=warnings,
    )
from typing import Any, Dict, List, Literal, Optional
from pydantic import BaseModel, Field


class CellStateRequest(BaseModel):
    assay: Literal["scrna", "perturbseq"]

    metadata_path: str

    # If trusted annotations already exist
    cell_id_column: str = "cell_id"
    state_column: Optional[str] = None
    lineage_column: Optional[str] = None

    # scRNA-specific
    condition_column: Optional[str] = None

    # Perturb-seq-specific
    perturbation_column: Optional[str] = "target"
    mixscape_class_column: Optional[str] = "mixscape_class"
    control_label: Optional[str] = "CTRL"

    # Annotation source
    annotation_mode: Literal[
        "provided",
        "reference_mapping",
        "marker_rules"
    ] = "provided"

    reference_path: Optional[str] = None
    marker_rules_path: Optional[str] = None


class StateCount(BaseModel):
    state: str
    n_cells: int
    fraction: float


class GroupStateCount(BaseModel):
    group: str
    state: str
    n_cells: int
    fraction_within_group: float


class CellStateResponse(BaseModel):
    assay: str
    status: Literal["ok", "needs_annotation", "error"]

    annotation_mode: str
    annotation_source: Optional[str] = None

    n_cells: int = 0
    states: List[str] = Field(default_factory=list)

    overall_state_summary: List[StateCount] = Field(default_factory=list)
    group_state_summary: List[GroupStateCount] = Field(default_factory=list)

    annotated_metadata_path: Optional[str] = None
    handoff_path: Optional[str] = None

    warnings: List[str] = Field(default_factory=list)
    metadata: Dict[str, Any] = Field(default_factory=dict)

from __future__ import annotations
from typing import Any, Literal
from pydantic import BaseModel, Field

RankingMode = Literal["balanced", "evidence", "connectivity"]
ProjectionMethod = Literal["jaccard", "shared_count", "weighted_shared"]
LayoutMode = Literal["bipartite", "force"]

class BuildGraphOptions(BaseModel):
    preset: str = Field("custom_long", description="custom_long, enrichr, gprofiler, clusterprofiler, or gsea_msigdb")
    item_col: str | None = None
    group_col: str | None = None
    weight_col: str | None = None
    apply_preset: bool = True

class FilterOptions(BaseModel):
    search: str = ""
    min_degree: int = 0
    min_weight: float = 0.0
    max_groups: int = 50
    largest_component_only: bool = False
    layout_mode: LayoutMode = "bipartite"
    show_labels: bool = False
    thickness_by_weight: bool = False
    return_figure: bool = True

class GraphStoreRequest(BaseModel):
    graph: dict[str, Any]

class FilterGraphRequest(GraphStoreRequest):
    options: FilterOptions = Field(default_factory=FilterOptions)
    highlight_nodes: dict[str, Any] | None = None

class DiffusionRequest(GraphStoreRequest):
    seed_node: str | None = None
    alpha: float = 0.85
    ranking_mode: RankingMode = "balanced"
    top_n: int = 50
    candidate_top_n: int = 30
    candidate_node_type: Literal["group", "item"] = "group"

class ProjectionBuildRequest(GraphStoreRequest):
    method: ProjectionMethod = "jaccard"
    return_figure: bool = True
    show_labels: bool = True

class ProjectionDiffusionRequest(BaseModel):
    projection_graph: dict[str, Any]
    seed_node: str | None = None
    alpha: float = 0.85
    ranking_mode: RankingMode = "balanced"
    top_n: int = 50
    candidate_top_n: int = 30

class ConsensusRequest(BaseModel):
    bipartite_results: dict[str, Any]
    projection_results: dict[str, Any]
    top_n: int = 30

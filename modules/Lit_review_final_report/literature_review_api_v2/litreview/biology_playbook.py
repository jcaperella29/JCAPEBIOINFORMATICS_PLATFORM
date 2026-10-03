from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Any
import re

import yaml


@dataclass(frozen=True)
class BiologyEdge:
    source: str
    relation: str
    target: str
    direction: str | None = None
    context: str | None = None
    caveat: str | None = None
    aliases: tuple[str, ...] = ()
    evidence_note: str | None = None


@dataclass(frozen=True)
class BiologyPath:
    nodes: tuple[str, ...]
    relations: tuple[str, ...]
    distance: int
    explanation: str


def _norm(text: str | None) -> str:
    return re.sub(r"[^a-z0-9+\- ]+", " ", (text or "").lower()).strip()


def _pretty(edge: BiologyEdge) -> str:
    bits = [edge.source, edge.relation.lower().replace("_", " "), edge.target]
    if edge.direction:
        bits.append(f"({edge.direction})")
    if edge.context:
        bits.append(f"[context: {edge.context}]")
    if edge.caveat:
        bits.append(f"[caveat: {edge.caveat}]")
    return " ".join(bits)


class PlaybookConfigError(ValueError):
    pass


class BiologyPlaybook:
    def __init__(
        self,
        edges: Iterable[BiologyEdge],
        *,
        direct_relations: set[str],
        mechanistic_relations: set[str],
        context_only_relations: set[str],
        source_path: str | None = None,
        version: str | None = None,
    ):
        self.edges = tuple(edges)
        self.direct_relations = set(direct_relations)
        self.mechanistic_relations = set(mechanistic_relations)
        self.context_only_relations = set(context_only_relations)
        self.source_path = source_path
        self.version = version

        self._adj: dict[str, list[BiologyEdge]] = {}
        for edge in self.edges:
            terms = [edge.source, *edge.aliases]
            for term in terms:
                nterm = _norm(term)
                if nterm:
                    self._adj.setdefault(nterm, []).append(edge)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "BiologyPlaybook":
        path = Path(path)
        with path.open("r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh) or {}

        policies = raw.get("relation_policy") or {}
        direct = set(policies.get("direct_relations") or [])
        mechanistic_extra = set(policies.get("mechanistic_only_relations") or [])
        context_only = set(policies.get("context_only_relations") or [])
        mechanistic = direct | mechanistic_extra

        if not direct:
            raise PlaybookConfigError("relation_policy.direct_relations must not be empty")
        if not context_only:
            raise PlaybookConfigError("relation_policy.context_only_relations must not be empty")

        edges_raw = raw.get("edges") or []
        edges: list[BiologyEdge] = []
        valid_relations = direct | mechanistic_extra | context_only
        for i, item in enumerate(edges_raw, start=1):
            if not isinstance(item, dict):
                raise PlaybookConfigError(f"Edge #{i} must be a mapping")
            for required in ("source", "relation", "target"):
                if not item.get(required):
                    raise PlaybookConfigError(f"Edge #{i} is missing required field: {required}")
            relation = str(item["relation"]).strip().upper()
            if relation not in valid_relations:
                raise PlaybookConfigError(
                    f"Edge #{i} uses unknown relation {relation!r}. Add it to relation_policy first."
                )
            aliases = tuple(str(x) for x in (item.get("aliases") or []))
            edges.append(
                BiologyEdge(
                    source=str(item["source"]),
                    relation=relation,
                    target=str(item["target"]),
                    direction=item.get("direction"),
                    context=item.get("context"),
                    caveat=item.get("caveat"),
                    aliases=aliases,
                    evidence_note=item.get("evidence_note"),
                )
            )

        if not edges:
            raise PlaybookConfigError("Playbook must contain at least one edge")

        return cls(
            edges,
            direct_relations=direct,
            mechanistic_relations=mechanistic,
            context_only_relations=context_only,
            source_path=str(path),
            version=str(raw.get("version")) if raw.get("version") is not None else None,
        )

    def terms_reachable_from(self, seeds: Iterable[str], max_distance: int = 2) -> set[str]:
        out: set[str] = set()
        queue = deque()
        seen: set[tuple[str, int]] = set()
        for s in seeds:
            ns = _norm(s)
            if ns:
                queue.append((ns, 0))
                seen.add((ns, 0))
        while queue:
            node, d = queue.popleft()
            if d >= max_distance:
                continue
            for edge in self._adj.get(node, []):
                target = _norm(edge.target)
                if target:
                    out.add(edge.target)
                    state = (target, d + 1)
                    if state not in seen:
                        seen.add(state)
                        queue.append(state)
        return out

    def find_path_to_text(
        self,
        seeds: Iterable[str],
        paper_text: str,
        *,
        max_distance: int,
        allowed_relations: set[str],
    ) -> BiologyPath | None:
        """Find the shortest allowed mechanistic path from a finding term to a term named in a paper."""
        text = _norm(paper_text)
        queue = deque()
        seen: set[str] = set()

        for seed in seeds:
            ns = _norm(seed)
            if not ns:
                continue
            queue.append((ns, (seed,), (), ()))
            seen.add(ns)

        while queue:
            node, nodes, rels, pretty_edges = queue.popleft()
            distance = len(rels)
            if distance >= max_distance:
                continue

            for edge in self._adj.get(node, []):
                if edge.relation not in allowed_relations:
                    continue
                target_norm = _norm(edge.target)
                new_nodes = nodes + (edge.target,)
                new_rels = rels + (edge.relation,)
                new_pretty = pretty_edges + (_pretty(edge),)
                new_distance = len(new_rels)

                if len(target_norm) >= 4 and target_norm in text:
                    return BiologyPath(
                        nodes=new_nodes,
                        relations=new_rels,
                        distance=new_distance,
                        explanation=" -> ".join(new_pretty),
                    )

                if new_distance < max_distance and target_norm not in seen:
                    seen.add(target_norm)
                    queue.append((target_norm, new_nodes, new_rels, new_pretty))
        return None

    def metadata(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "source_path": self.source_path,
            "edge_count": len(self.edges),
            "direct_relations": sorted(self.direct_relations),
            "mechanistic_relations": sorted(self.mechanistic_relations),
            "context_only_relations": sorted(self.context_only_relations),
        }


DEFAULT_PLAYBOOK_PATH = Path(__file__).resolve().parent.parent / "config" / "biology_playbook.yaml"
DEFAULT_PLAYBOOK = BiologyPlaybook.from_yaml(DEFAULT_PLAYBOOK_PATH)
DIRECT_RELATIONS = DEFAULT_PLAYBOOK.direct_relations
MECHANISTIC_RELATIONS = DEFAULT_PLAYBOOK.mechanistic_relations
CONTEXT_ONLY_RELATIONS = DEFAULT_PLAYBOOK.context_only_relations

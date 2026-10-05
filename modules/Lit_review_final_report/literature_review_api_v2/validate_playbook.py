"""Validate the editable biology playbook without starting the API."""
from litreview.biology_playbook import BiologyPlaybook, DEFAULT_PLAYBOOK_PATH

pb = BiologyPlaybook.from_yaml(DEFAULT_PLAYBOOK_PATH)
print(f"OK: {pb.source_path}")
print(f"version={pb.version} edges={len(pb.edges)}")
print("direct relations:", ", ".join(sorted(pb.direct_relations)))
print("mechanistic relations:", ", ".join(sorted(pb.mechanistic_relations)))
print("context-only relations:", ", ".join(sorted(pb.context_only_relations)))

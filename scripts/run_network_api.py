#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import subprocess
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

def post_json(url: str, payload: dict) -> dict:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=600) as r:
        return json.loads(r.read().decode("utf-8"))

def find_input(root: Path) -> Path:
    preferred = [
        "enrichment_all_network_edges.csv",
        "enrichment_network_edges.csv",
    ]
    files = list(root.rglob("*.csv"))
    for name in preferred:
        hits = [p for p in files if p.name == name]
        if hits:
            return sorted(hits)[0]
    hits = sorted(p for p in files if "network_edges" in p.name.lower() and "all" in p.name.lower())
    if hits:
        return hits[0]
    hits = sorted(p for p in files if "network_edges" in p.name.lower())
    if hits:
        return hits[0]
    raise SystemExit("No RNA-seq enrichment network-edge CSV was found.")

def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")

def write_text(path: Path, value: str) -> None:
    path.write_text(value or "", encoding="utf-8")

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--api-url", required=True)
    p.add_argument("--rnaseq-dir", required=True)
    p.add_argument("--ranking-mode", default="balanced")
    p.add_argument("--projection-method", default="jaccard")
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--candidate-top-n", type=int, default=30)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    edge_csv = find_input(Path(a.rnaseq_dir))
    (out / "selected_input.txt").write_text(str(edge_csv), encoding="utf-8")

    options = json.dumps({
        "apply_preset": False,
        "item_col": "gene",
        "group_col": "term",
        "weight_col": "adjusted_pvalue",
    })
    build_cmd = [
        "curl", "--fail-with-body", "--silent", "--show-error",
        "-X", "POST", f"{a.api_url.rstrip('/')}/network/build",
        "-F", f"file=@{edge_csv}",
        "-F", f"options_json={options}",
    ]
    built = json.loads(subprocess.run(build_cmd, check=True, text=True, capture_output=True).stdout)
    write_json(out / "network_build.json", built)
    write_text(out / "main_nodes.csv", built.get("exports", {}).get("nodes_csv", ""))
    write_text(out / "main_edges.csv", built.get("exports", {}).get("edges_csv", ""))

    bip = post_json(f"{a.api_url.rstrip('/')}/diffusion/bipartite", {
        "graph": built["graph"],
        "seed_node": None,
        "alpha": 0.85,
        "ranking_mode": a.ranking_mode,
        "top_n": a.top_n,
        "candidate_top_n": a.candidate_top_n,
        "candidate_node_type": "group",
    })
    write_json(out / "bipartite_diffusion.json", bip)
    write_text(out / "bipartite_diffusion_results.csv", bip.get("csv", ""))
    write_text(out / "bipartite_top_candidates.csv", bip.get("candidate_csv", ""))

    projection = post_json(f"{a.api_url.rstrip('/')}/projection/build", {
        "graph": built["graph"],
        "method": a.projection_method,
        "return_figure": True,
        "show_labels": True,
    })
    write_json(out / "projection_build.json", projection)
    write_text(out / "projection_nodes.csv", projection.get("exports", {}).get("nodes_csv", ""))
    write_text(out / "projection_edges.csv", projection.get("exports", {}).get("edges_csv", ""))

    proj_diff = post_json(f"{a.api_url.rstrip('/')}/diffusion/projection", {
        "projection_graph": projection["projection_graph"],
        "seed_node": None,
        "alpha": 0.85,
        "ranking_mode": a.ranking_mode,
        "top_n": a.top_n,
        "candidate_top_n": a.candidate_top_n,
    })
    write_json(out / "projection_diffusion.json", proj_diff)
    write_text(out / "projection_diffusion_results.csv", proj_diff.get("csv", ""))
    write_text(out / "projection_top_candidates.csv", proj_diff.get("candidate_csv", ""))

    consensus = post_json(f"{a.api_url.rstrip('/')}/consensus", {
        "bipartite_results": bip,
        "projection_results": proj_diff,
        "top_n": a.candidate_top_n,
    })
    write_json(out / "consensus.json", consensus)
    write_text(out / "consensus_candidates.csv", consensus.get("csv", ""))

    bundle_payload = {
        "graph": built["graph"],
        "projection_graph": projection["projection_graph"],
        "bipartite_results": bip,
        "projection_results": proj_diff,
        "consensus_results": consensus,
        "mapped_columns": built.get("mapped_columns", {}),
        "settings": {
            "ranking_mode": a.ranking_mode,
            "projection_method": a.projection_method,
            "top_n": a.top_n,
            "candidate_top_n": a.candidate_top_n,
        },
    }
    payload_path = out / "report_bundle_payload.json"
    write_json(payload_path, bundle_payload)
    subprocess.run([
        "curl", "--fail-with-body", "--silent", "--show-error",
        "-X", "POST", f"{a.api_url.rstrip('/')}/export/report-bundle",
        "-H", "Content-Type: application/json",
        "--data-binary", f"@{payload_path}",
        "-o", str(out / "network_report_bundle.zip"),
    ], check=True)

    manifest = {
        "stage": "network",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "api_url": a.api_url,
        "selected_rnaseq_input": str(edge_csv),
        "ranking_mode": a.ranking_mode,
        "projection_method": a.projection_method,
    }
    write_json(out / "pipeline_stage_manifest.json", manifest)

if __name__ == "__main__":
    main()

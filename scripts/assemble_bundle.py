#!/usr/bin/env python3
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
import zipfile
from pathlib import Path
from datetime import datetime, timezone

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def copy_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rnaseq-dir", required=True)
    p.add_argument("--network-dir", required=True)
    p.add_argument("--llm-dir", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--phenotype", required=True)
    p.add_argument("--context-json", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--zip-path", required=True)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    copy_tree(Path(a.rnaseq_dir), out / "01_rnaseq")
    copy_tree(Path(a.network_dir), out / "02_network")
    copy_tree(Path(a.llm_dir), out / "03_llm_triage")

    inputs = out / "00_inputs"
    inputs.mkdir(exist_ok=True)
    input_records = []
    for source, label in [
        (Path(a.counts), "counts"),
        (Path(a.phenotype), "phenotype"),
        (Path(a.context_json), "context"),
    ]:
        dest = inputs / source.name
        shutil.copy2(source, dest)
        input_records.append({
            "label": label,
            "file": str(dest.relative_to(out)),
            "sha256": sha256(dest),
            "bytes": dest.stat().st_size,
        })

    manifest = {
        "product": "JCAP integrated RNA-seq, network, and LLM triage pipeline",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "layout": {
            "00_inputs": "Original inputs and experiment context",
            "01_rnaseq": "RNA-seq API outputs",
            "02_network": "Network diffusion, projection, and consensus outputs",
            "03_llm_triage": "Deterministic or full PubMed/LLM interpretation outputs",
        },
        "inputs": input_records,
        "caveats": [
            "Network diffusion and consensus scores are prioritization scores, not statistical p-values.",
            "The LLM interpretation should be treated as an evidence-organizing layer, not a causal conclusion.",
            "Original adjusted p-values and RNA-seq model outputs remain the statistical source of truth.",
        ],
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    zip_path = Path(a.zip_path)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for f in sorted(out.rglob("*")):
            if f.is_file():
                z.write(f, arcname=f.relative_to(out.parent))

if __name__ == "__main__":
    main()

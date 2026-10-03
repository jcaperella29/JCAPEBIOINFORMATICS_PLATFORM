#!/usr/bin/env python3

import json
import shutil
import sys
from pathlib import Path


def copy_any(src, dst):
    src = Path(src)
    dst = Path(dst)

    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main():
    if len(sys.argv) != 6:
        raise SystemExit(
            "Usage: bundle.py SCRNA_DIR NETWORK_DIR LLM_DIR METADATA CONTEXT"
        )

    scrna_dir, network_dir, llm_dir, metadata, context = map(Path, sys.argv[1:])

    out = Path("jcap_scrna_integrated_results")
    if out.exists():
        shutil.rmtree(out)

    out.mkdir(parents=True)

    copy_any(scrna_dir, out / "01_scrna")
    copy_any(network_dir, out / "02_network")
    copy_any(llm_dir, out / "03_llm_triage")

    inputs = out / "inputs"
    inputs.mkdir()

    shutil.copy2(metadata, inputs / "metadata_stateaware.csv")
    shutil.copy2(context, inputs / "context_stateaware.json")

    manifest = {
        "pipeline": "jcap_scrna_integrated_pipeline",
        "scrna_results": "01_scrna",
        "network_results": "02_network",
        "llm_triage": "03_llm_triage",
        "metadata": "inputs/metadata_stateaware.csv",
        "context": "inputs/context_stateaware.json"
    }

    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

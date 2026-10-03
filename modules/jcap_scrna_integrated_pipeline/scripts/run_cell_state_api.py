#!/usr/bin/env python3

import argparse
import json
import shutil
from pathlib import Path

import requests


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--metadata", required=True)
    ap.add_argument("--api-url", required=True)
    ap.add_argument("--assay", required=True)

    ap.add_argument("--cell-id-column", default="cell_id")
    ap.add_argument("--state-column", default="")
    ap.add_argument("--lineage-column", default="")
    ap.add_argument("--condition-column", default="")
    ap.add_argument("--perturbation-column", default="")
    ap.add_argument("--mixscape-class-column", default="")
    ap.add_argument("--control-label", default="CTRL")

    ap.add_argument(
        "--annotation-mode",
        default="provided",
        choices=["provided", "reference_mapping", "marker_rules"]
    )

    ap.add_argument("--out-metadata", default="metadata_stateaware.csv")
    ap.add_argument("--out-handoff", default="cell_state_handoff.json")

    args = ap.parse_args()

    payload = {
        "assay": args.assay,
        "metadata_path": str(Path(args.metadata).resolve()),
        "cell_id_column": args.cell_id_column,
        "state_column": args.state_column or None,
        "lineage_column": args.lineage_column or None,
        "condition_column": args.condition_column or None,
        "perturbation_column": args.perturbation_column or None,
        "mixscape_class_column": args.mixscape_class_column or None,
        "control_label": args.control_label,
        "annotation_mode": args.annotation_mode,
    }

    url = args.api_url.rstrip("/") + "/analyze"

    r = requests.post(
        url,
        json=payload,
        timeout=300
    )

    r.raise_for_status()

    result = r.json()

    Path(args.out_handoff).write_text(
        json.dumps(result, indent=2) + "\n"
    )

    status = result.get("status")

    if status == "ok":
        annotated = result.get("annotated_metadata_path")

        if not annotated:
            raise RuntimeError(
                "Cell-state API returned status=ok but no annotated_metadata_path"
            )

        shutil.copyfile(
            annotated,
            args.out_metadata
        )

    elif status == "needs_annotation":
        shutil.copyfile(
            args.metadata,
            args.out_metadata
        )

    else:
        raise RuntimeError(
            f"Cell-state API failed with status={status}: {result}"
        )

    print(
        f"Cell-state analysis status: {status}"
    )


if __name__ == "__main__":
    main()

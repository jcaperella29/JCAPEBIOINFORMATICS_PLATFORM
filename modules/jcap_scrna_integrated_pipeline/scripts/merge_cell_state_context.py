#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--context", required=True)
    ap.add_argument("--cell-state-handoff", required=True)
    ap.add_argument("--out", default="context_stateaware.json")
    args = ap.parse_args()

    context = json.loads(Path(args.context).read_text())
    cell_state = json.loads(Path(args.cell_state_handoff).read_text())

    context["cell_state_analysis"] = cell_state

    context["cell_state_reporting_instruction"] = (
        "When summarizing perturbation biology, use the quantitative "
        "group_state_summary from cell_state_analysis. If discussing lineage "
        "or cell-state composition, report the relevant state names together "
        "with their cell counts and percentages where available. Do not merely "
        "say that composition shifted. Distinguish this perturbation's state "
        "distribution from cross-target comparisons, and do not call a pathway "
        "a direct target-gene mechanism when the same signal could be explained "
        "by state-composition differences."
    )

    Path(args.out).write_text(
        json.dumps(context, indent=2) + "\n"
    )

    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

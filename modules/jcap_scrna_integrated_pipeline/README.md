# JCAP Kang scRNA Integrated Pipeline + Literature Review

This patch extends the working Kang scRNA integrated pipeline:

1. scRNA pseudobulk analysis
2. Network analysis
3. Full LLM triage
4. Literature Review + Final Report API
5. Final integrated bundle

## API ports

Start the existing scRNA, Network, and LLM APIs as before, plus the Literature Review API on port 8006.

The supplied `params.json` uses the same WSL/host gateway pattern as the working Kang parameters:

- scRNA: `http://172.22.32.1:8003`
- Network: `http://172.22.32.1:8001`
- LLM: `http://172.22.32.1:8004`
- Literature: `http://172.22.32.1:8006`

The Kang analysis settings remain unchanged: `stim`, `orig.ident`, STIM vs CTRL, target `CD14+ Monocytes`, full LLM mode.

## Install

Copy these files into the existing `jcap_scrna_integrated_pipeline` directory:

- `main.nf`
- `params.json`
- `nextflow.config`
- `scripts/run_literature_api.py`

Keep the existing `scripts/check.py`, `scripts/scrna.py`, `scripts/network.py`, `scripts/llm.py`, and `scripts/bundle.py`.

## Start literature service

Using the established JCAP literature module:

```bash
cd "/mnt/c/Users/jcape/Downloads/omics_dashboard/omics_pipeline/modules/Lit_review_final_report"
set -a
source ../Enrichement_triage/.env
set +a
python -m uvicorn app:app --host 0.0.0.0 --port 8006
```

## Run

```bash
nextflow run main.nf -params-file params.json -resume
```

The new stage writes `04_literature_review/`, including `literature_review_request.json` and the extracted Literature Review report artifacts. The final bundle also contains `04_literature_review/`.

## First validation

Before trusting the prose, inspect the new request:

```bash
python -m json.tool results/jcap_scrna_integrated_results/04_literature_review/literature_review_request.json | less
```

Confirm that the handoffs include scRNA, Network, and LLM/Triage data and that scRNA context retains the Kang cell type/condition/replicate information.

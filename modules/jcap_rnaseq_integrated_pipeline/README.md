# JCAP Integrated RNA-seq Pipeline

This Nextflow workflow orchestrates three already-running HTTP services:

1. RNA-seq API
2. Enrichment Network API
3. Enrichment LLM Triage API

It returns one ZIP with labeled folders:

- `00_inputs/`
- `01_rnaseq/`
- `02_network/`
- `03_llm_triage/`
- `manifest.json`

## Why this first version uses HTTP

Your three services already work and have stable API contracts. This version therefore gives you an immediately testable end-to-end workflow without rewriting the analysis code. A later version can replace each HTTP call with direct container execution for cloud batch systems.

## Requirements

- Java 17 or newer
- Nextflow
- Python 3
- `curl`
- All three APIs running and reachable

Default URLs:

- RNA-seq: `http://127.0.0.1:8000`
- Network: `http://127.0.0.1:8001`
- LLM triage: `http://127.0.0.1:8002`

The URLs are configurable in `params.example.json`.

## Inputs

### Counts

CSV or whitespace-delimited raw count matrix accepted by your RNA-seq CLI. Genes are rows and samples are columns.

### Phenotype

CSV accepted by the RNA-seq CLI. Sample IDs must be row names / first column and match the count-matrix columns.

### Context JSON

Use `data/context.example.json` as the template. `phenotype` is required for the LLM triage request. Other fields are forwarded as assay and experiment context.

## Run

Copy and edit the examples:

```bash
cp params.example.json params.json
cp data/context.example.json data/context.json
```

Then:

```bash
nextflow run main.nf -params-file params.json -resume
```

The final deliverables are published to the configured `outdir`, normally:

```text
results/
├── jcap_integrated_results/
├── jcap_integrated_results.zip
├── nextflow_report.html
├── nextflow_timeline.html
├── nextflow_trace.tsv
└── nextflow_dag.html
```

## Test API health first

```bash
curl http://127.0.0.1:8000/health
curl http://127.0.0.1:8001/health
curl http://127.0.0.1:8002/health
```

## LLM modes

Use:

```json
"llm_mode": "pretriage"
```

for deterministic validation, ranking, and program grouping without OpenAI or PubMed.

Use:

```json
"llm_mode": "full"
```

for the full PubMed/OpenAI interpretation. The LLM API container must receive `OPENAI_API_KEY` and the relevant NCBI configuration.

## Resume behavior

Nextflow caches completed processes. After fixing a failed API or changing only a downstream parameter, rerun with:

```bash
nextflow run main.nf -params-file params.json -resume
```

## Current boundary

The workflow calls long-running API containers over HTTP. It does not start or stop those containers. The next hardening step is a direct-container profile for AWS Batch, Google Batch, Slurm, or Kubernetes.

## Included airway demonstration context

The bundle now includes:

- `data/context.json`
- `data/airway_context.json`
- `params.airway.example.json`

The context is matched to the Bioconductor airway experiment: primary human airway smooth muscle
cells, 1 micromolar dexamethasone for 18 hours, with treated and untreated samples from each of
four cell lines.

The LLM API mapping is:

- `phenotype` -> required phenotype form field
- `organism`
- `assay`
- `tissue`
- `cell_type`
- `perturbation`
- `timepoint`
- every remaining JSON key -> `extra_context_json`

Before running, confirm that the phenotype CSV actually uses:

- phenotype column: `dex`
- treated value: `trt`
- reference value: `untrt`

If your exported phenotype table uses different names, edit `params.json` and the context labels
to match the file rather than renaming biological groups silently.

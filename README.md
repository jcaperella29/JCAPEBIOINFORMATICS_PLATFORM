# JCAP Bioinformatics Platform

An integrated, pipeline-first bioinformatics platform for **bulk RNA-seq**, **single-cell RNA-seq**, and **Perturb-seq** analysis, followed by enrichment-network prioritization, evidence-aware biological triage, literature review, and structured reporting.

The platform is built around a simple idea:

> **Numbers. Narrative. Next Steps.**

First generate and validate the quantitative evidence. Then turn that evidence into a biological narrative. Finally, use the combined statistical, network, and literature context to prioritize what to investigate next.

Raw omics data are not handed directly to an LLM. Statistical analysis, enrichment, network prioritization, deterministic biological rules, and literature retrieval provide the evidence base. The LLM layer is used later to synthesize that evidence into cautious, auditable interpretation and follow-up suggestions.

---

## Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Supported workflows](#supported-workflows)
- [Core components](#core-components)
- [Service ports](#service-ports)
- [Requirements](#requirements)
- [Running the integrated pipelines](#running-the-integrated-pipelines)
- [Dashboard](#dashboard)
- [Network rendering and performance](#network-rendering-and-performance)
- [Evidence-aware interpretation](#evidence-aware-interpretation)
- [Literature evidence](#literature-evidence)
- [Outputs](#outputs)
- [Bundled test data](#bundled-test-data)
- [Configuration](#configuration)
- [Development notes](#development-notes)
- [Scientific-use note](#scientific-use-note)

---

## Overview

JCAP Bioinformatics Platform connects several analysis layers into reproducible Nextflow workflows:

```text
                         JCAP Bioinformatics Platform

   ┌───────────────┐    ┌────────────────┐    ┌─────────────────┐
   │ Bulk RNA-seq  │    │   scRNA-seq    │    │   Perturb-seq   │
   └───────┬───────┘    └───────┬────────┘    └────────┬────────┘
           │                     │                     │
           ▼                     ▼                     ▼
      Omics analysis       Omics + cell state    Mixscape + cell state
           │                     │                     │
           └───────────────┬─────┴─────────────────────┘
                           ▼
                  Enrichment Network
                           │
                           ▼
                Evidence-aware Triage
                           │
                           ▼
                   Literature Review
                           │
                           ▼
                Structured Final Report
                           │
                           ▼
                     Dash Dashboard
```

The three primary integrated workflows share downstream infrastructure while preserving assay-specific analysis.

---

## Architecture

The platform is composed of independent analysis/API modules orchestrated by Nextflow and surfaced through a Dash interface.

### High-level flow

```text
Input data
   │
   ├── Bulk RNA-seq ──► differential expression + enrichment
   │
   ├── scRNA-seq ─────► Seurat analysis + DE + enrichment + cell-state context
   │
   └── Perturb-seq ───► Mixscape + target-specific DE/enrichment + cell-state context
   │
   ▼
Gene/pathway network construction
   │
   ▼
Bipartite diffusion + pathway projection + candidate prioritization
   │
   ▼
Deterministic biological triage / BioFit / program summarization
   │
   ▼
Context retrieval + evidence-aware LLM synthesis
   │
   ▼
Independent literature review / evidence classification
   │
   ▼
Final JSON + Markdown reports
   │
   ▼
Dash visualization and result browsing
```

---

## Supported workflows

### 1. Bulk RNA-seq

The integrated bulk RNA-seq workflow combines:

1. RNA-seq quantitative analysis
2. differential expression
3. enrichment analysis
4. enrichment-network construction
5. bipartite and pathway-projection ranking
6. evidence-aware LLM triage
7. literature review
8. final structured report

The bundled configuration uses the airway example dataset and supports parameters including FDR, log-fold-change cutoff, enrichment backend, network ranking mode, and projection method.

### 2. Single-cell RNA-seq

The scRNA-seq workflow supports:

- counts + metadata input
- Seurat-based analysis
- PCA / UMAP
- cell-level or pseudobulk differential expression
- configurable condition and cell-type comparisons
- condition-only DE analysis
- enrichment
- classifier / feature-prioritization outputs
- power-analysis outputs
- cell-state context
- enrichment-network analysis
- evidence-aware triage
- literature review
- final reporting

The bundled Kang configuration uses pseudobulk differential expression in a target cell type.

### 3. Perturb-seq

The Perturb-seq workflow is target-aware and is built around Mixscape.

At a high level it performs:

```text
API checks
  │
  ▼
Whole-dataset Mixscape classification
  │
  ├──► Cell-state analysis
  │       │
  │       ▼
  │    State-aware context
  │
  ▼
Select requested perturbations that are actually available
  │
  ▼
For each supported target:
  ├── Target-specific Mixscape / DE / enrichment
  ├── Network analysis
  ├── LLM triage
  └── Literature review
  │
  ▼
Assemble target bundles
  │
  ▼
Collect targets
  │
  ▼
Multi-target final report
```

The bundled Giladi subset is configured for:

```text
CEBPA
KLF4
FLT3
FCGR3
CSF1R
IRF8
```

with mouse data, `CTRL` controls, and provided `lineage_state` annotations.

---

## Core components

### RNA-seq module

Provides the assay-specific bulk RNA-seq analysis used by the integrated RNA workflow.

### scRNA-seq module

The single-cell module is an R/Seurat analysis pipeline with configurable statistical and reporting behavior.

Supported options include:

- cell-level or pseudobulk DE
- global, top-cell-type, or selected-cell-type DE scopes
- Wilcoxon and other Seurat-supported tests
- limma-voom pseudobulk analysis
- g:Profiler / Enrichr enrichment
- PCA and UMAP
- feature selection / classification
- sample-aware validation
- power analysis
- reusable tables, HTML plots, and R objects

### Cell-State API

Provides a common cell-state / lineage interface for scRNA-seq and Perturb-seq.

Current `provided` annotation mode can:

- validate a state column
- calculate overall state frequencies
- summarize state frequencies by condition for scRNA-seq
- summarize state frequencies by perturbation for Perturb-seq
- emit annotated metadata
- emit a reusable JSON handoff

`reference_mapping` and `marker_rules` are declared extension modes but are not currently implemented.

### Enrichment Network API

Converts enrichment results into graph-based prioritization.

Supported input presets include:

- custom long-format gene ↔ pathway tables
- Enrichr
- g:Profiler / gprofiler2
- clusterProfiler
- GSEA / MSigDB

The API builds a weighted bipartite gene/pathway network. When adjusted p-values are supplied, they are converted to:

```text
weight = -log10(adjusted p-value)
```

Network analysis includes:

- filtering by degree, edge weight, search, and component
- graph statistics
- weighted PageRank / personalized PageRank diffusion
- separate pathway and gene rankings
- pathway-only projection graphs
- projection diffusion
- candidate prioritization
- consensus ranking
- CSV / ZIP export

Three diffusion ranking modes are available:

| Mode | Behavior |
|---|---|
| `evidence` | emphasizes enrichment evidence |
| `connectivity` | emphasizes graph topology |
| `balanced` | combines evidence and connectivity while damping extreme single-edge effects |

Candidate scores are **prioritization scores**, not statistical p-values.

### Evidence-Aware Enrichment Triage

The triage layer is deliberately more than a raw LLM call.

```text
enrichment input
   │
   ▼
schema validation / normalization
   │
   ▼
deterministic statistical + BioFit triage
   │
   ▼
biological program summarization
   │
   ├──► PubMed context retrieval
   │
   └──► biological playbook retrieval
   │
   ▼
evidence-aware LLM reasoning
   │
   ▼
structured interpretation claims
```

Input handling supports enrichment-style tables including Enrichr-like, clusterProfiler/GSEA-like, and network-derived tables.

The deterministic layer uses pathway terms, gene-family support, phenotype context, overlap information, statistical evidence, and explicit biological program rules before model synthesis.

Structured claims can include:

- likely driver
- likely reactive
- likely artifact / confounded
- uncertain

along with:

- evidence strength
- supporting terms and genes
- linked PMIDs
- rationale
- possible confounders
- validation experiments
- controls and readouts

The reasoning layer is instructed not to claim causality from enrichment alone.

### Biological playbooks

The triage system can retrieve assay- and biology-specific guidance rather than injecting an entire playbook into every prompt.

Current assay-aware paths include support for areas such as:

- bulk RNA-seq
- scRNA-seq
- Perturb-seq
- ATAC-seq
- miRNA-seq
- GWAS
- DNA methylation

OpenAI vector-store indexing utilities are also included for optional file-search/RAG workflows.

### Literature Review API

The literature layer is separate from the triage model.

It retrieves candidate papers from:

- PubMed
- Europe PMC
- OpenAlex

and classifies evidence into:

```text
direct
mechanistic_support
contextual_support
background
exclude
```

Classification considers assay match, target/gene match, pathway or phenotype match, organism/context information, review status, and short paths through the configurable biology playbook.

Retrieval is intentionally permissive; evidence classification is stricter.

This separation helps avoid two failure modes:

1. retrieval that is so strict that useful literature disappears
2. retrieval that is so permissive that loosely related papers look like direct evidence

---

## Service ports

The integrated workflows expect the following local services by default:

| Service | Default port |
|---|---:|
| Bulk RNA-seq API | `8000` |
| Enrichment Network API | `8001` |
| scRNA-seq API | `8003` |
| Evidence-aware triage API | `8004` |
| Literature Review API | `8006` |
| Cell-State API | `8007` |
| Dash dashboard | `8050` |

The integrated parameter files use `127.0.0.1` URLs by default.

Before running a workflow, make sure its required services are running and their `/health` endpoints respond successfully.

---

## Requirements

The platform currently spans Python, R, and Nextflow.

At minimum, expect to need:

- Git
- Python
- R
- Java
- Nextflow
- `curl`
- assay-specific R packages
- module-specific Python dependencies

The scRNA module uses packages including:

```text
Seurat
jsonlite
ggplot2
plotly
htmlwidgets
gprofiler2
randomForest
caret
pROC
pwr
```

The triage service uses packages including:

```text
FastAPI
Uvicorn
pandas
numpy
scipy
statsmodels
scikit-learn
openai
requests
python-dotenv
reportlab
```

Install the dependency/environment files supplied with each module where present.

### OpenAI configuration

Full LLM triage and optional vector-store/RAG functionality require an OpenAI API key supplied through the environment.

For example:

```bash
export OPENAI_API_KEY="..."
```

Do **not** commit API keys or local `.env` files to the repository.

---

## Running the integrated pipelines

Run commands from the corresponding integrated-pipeline directory.

### Bulk RNA-seq

```bash
cd modules/jcap_rnaseq_integrated_pipeline

nextflow run main.nf \
  -params-file params.json
```

The supplied parameters point to the bundled airway example data.

To override inputs:

```bash
nextflow run main.nf \
  -params-file params.json \
  --counts /path/to/counts.csv \
  --phenotype /path/to/phenotype.csv
```

### scRNA-seq

```bash
cd modules/jcap_scrna_integrated_pipeline_kang/jcap_scrna_integrated_pipeline

nextflow run main.nf \
  -params-file params.json
```

The supplied configuration points to the bundled Kang subset.

Typical custom execution:

```bash
nextflow run main.nf \
  -params-file params.json \
  --counts /path/to/counts.csv \
  --metadata /path/to/metadata.csv
```

### Perturb-seq

```bash
cd modules/jcap_perturbseq_integrated_pipeline/jcap_perturbseq_integrated_pipeline

nextflow run main_reconstructed.nf \
  -params-file params.giladi.json
```

The supplied parameters point to the bundled Giladi Mixscape subset.

The workflow expects Matrix Market counts plus cell, gene, metadata, state-metadata, context, and target information as configured in `params.giladi.json`.

---

## Dashboard

The Dash application provides one workspace for:

- launching integrated pipelines
- selecting built-in test data or uploaded data
- following Nextflow job status
- viewing staged results
- browsing assay-specific outputs

Run the dashboard from its module directory, for example:

```bash
python3 dashboard_v3_scrna_demo_wired.py
```

Then open:

```text
http://127.0.0.1:8050
```

The dashboard can use environment variables to override pipeline and data locations, which is preferable to editing source paths when moving the repository between machines.

### Output tabs

#### Bulk RNA-seq

```text
Omics
Network
Triage
Literature
Final Report
```

#### scRNA-seq

```text
Omics
Cell State
Network
Triage
Literature
Final Report
```

#### Perturb-seq

```text
Mixscape / Omics
Cell State
Network
Triage
Literature
Final Report
```

Plot-heavy sections use browsable HTML/Plotly outputs rather than forcing every figure onto the page simultaneously.

---

## Network rendering and performance

Network **analysis** and network **figure rendering** are intentionally separate.

By default:

```json
"render_network": false
```

The pipelines still build graphs, calculate diffusion rankings, generate projection graphs, export network tables, and produce candidate rankings. They do **not** render the network figure unless explicitly requested.

This is intentional because the dashboard can reconstruct and render the graph from the saved network data. Avoiding duplicate figure generation keeps the pipeline itself lighter and prevents unnecessary rendering work.

### Dashboard performance note

Large enrichment networks can contain many nodes and edges. When a network is large or structurally complex, rendering it interactively in the dashboard UI may take substantial time even though the underlying network analysis has already completed successfully.

For especially complex enrichment results, the **CLI / Nextflow pipeline outputs may be the better primary analysis path**. The exported node tables, edge tables, diffusion rankings, pathway projections, candidate rankings, and other network results can be inspected without waiting for a large interactive graph to render.

In other words:

> **A slow dashboard graph does not mean the network analysis is slow or incomplete. Rendering and analysis are separate steps.**

Use the dashboard when interactive visualization is useful, but prefer CLI/pipeline outputs for very large networks or when rapid access to the quantitative network results matters more than visualization.

To explicitly request API-side network figures:

### RNA-seq

```bash
nextflow run main.nf \
  -params-file params.json \
  --render_network true
```

### scRNA-seq

```bash
nextflow run main.nf \
  -params-file params.json \
  --render_network true
```

### Perturb-seq

```bash
nextflow run main_reconstructed.nf \
  -params-file params.giladi.json \
  --render_network true
```

The dashboard can reconstruct network figures from serialized graph data when appropriate.

For very large graphs, interactive browser rendering may still be slow. In those cases, prefer the pipeline/CLI outputs and exported node/edge tables, or filter the network before visualization.

Useful network-reduction strategies include:

- increasing minimum node degree
- increasing minimum edge weight
- reducing the number of displayed pathways
- searching for a pathway family
- avoiding large force-directed layouts

---

## Evidence-aware interpretation

The platform deliberately separates **measurement**, **prioritization**, and **interpretation**.

### Quantitative layer

Assay-specific statistical analysis produces DE, enrichment, and related quantitative outputs.

### Network layer

The network module asks questions such as:

- Which pathways are broadly supported by the gene/pathway structure?
- Which nodes remain important when statistical evidence and connectivity are considered together?
- Which pathways overlap through shared genes?
- Which candidates appear consistently across bipartite and projection analyses?

These rankings are intended to guide attention. They are not replacement p-values.

### Deterministic triage layer

BioFit and program summarization provide explicit, inspectable rules for recognizing biologically coherent programs, possible confounders, and gene-family support.

### LLM synthesis layer

The LLM receives already-structured quantitative, network, literature, and playbook context.

It is used to produce cautious biological interpretations and follow-up ideas rather than to replace upstream statistical analysis.

---

## Literature evidence

The platform uses two complementary literature mechanisms.

### Triage-time literature context

The triage pipeline retrieves PubMed context to help the reasoner interpret biological programs while retaining relevance/caution labels.

### Independent literature review

The Literature Review API performs its own retrieval and evidence classification across PubMed, Europe PMC, and OpenAlex.

This gives final reports a literature layer that can distinguish:

- same-assay / same-question evidence
- mechanistic support
- contextual support
- review/background material
- unrelated results that should be excluded

---

## Outputs

Exact output contents vary by assay, but the integrated workflows are organized around consistent stages.

### Bulk RNA-seq

Typical result families:

```text
omics
network
triage
literature
final report
```

### scRNA-seq

The scRNA analysis may emit tables such as:

```text
condition_de_table.csv
celltype_de_table.csv
condition_only_de_table.csv
pca_coordinates.csv
umap_coordinates.csv
enrichment_all.csv
enrichment_up.csv
enrichment_down.csv
classifier metrics / predictions / feature tables
power-analysis tables
```

and interactive HTML plots for PCA, UMAP, volcano plots, and enrichment.

A serialized Seurat object and run/manifest metadata are also produced for reproducibility.

### Perturb-seq

The integrated Perturb-seq result tree is target-oriented.

Conceptually:

```text
results/
└── jcap_perturbseq_integrated_results/
    ├── cell_state/
    ├── context_stateaware.json
    ├── target_selection.tsv
    ├── skipped_targets.tsv
    ├── targets/
    │   ├── CEBPA/
    │   │   └── target_CEBPA/
    │   │       ├── 01_mixscape/
    │   │       ├── 02_network/
    │   │       ├── 03_llm_triage/
    │   │       └── 04_literature_review/
    │   ├── KLF4/
    │   └── ...
    ├── perturbseq_final_report.json
    └── perturbseq_final_report.md
```

Only perturbations supported by the classified dataset proceed into target-specific analysis. Requested targets that cannot be analyzed are recorded separately.

### Final reports

Integrated workflows produce machine-readable and human-readable final report artifacts, centered on:

```text
JSON
Markdown
```

The dashboard exposes these outputs directly.

---

## Bundled test data

The repository includes compact datasets for exercising the integrated workflows.

### Bulk RNA-seq

Airway example data:

```text
counts_data.csv
phenotype_data.csv
```

### scRNA-seq

Kang subset:

```text
kang18_8vs8_counts_subset.csv
kang18_8vs8_metadata_subset.csv
```

The supplied integrated configuration uses a stimulated-vs-control comparison in `CD14+ Monocytes` with pseudobulk `limma_voom`.

### Perturb-seq

Giladi Mixscape subset:

```text
giladi_rna_counts.mtx
giladi_cells.tsv
giladi_genes.tsv
giladi_metadata.csv
giladi_metadata_with_lineage.csv
```

The bundled target list is:

```text
CEBPA,KLF4,FLT3,FCGR3,CSF1R,IRF8
```

---

## Configuration

Each integrated pipeline has a JSON parameter file.

### RNA-seq

```text
params.json
```

Example settings include:

```json
{
  "rnaseq_url": "http://127.0.0.1:8000",
  "network_url": "http://127.0.0.1:8001",
  "llm_url": "http://127.0.0.1:8004",
  "fdr": 0.05,
  "logfc_cutoff": 1.0,
  "network_ranking_mode": "balanced",
  "network_projection_method": "jaccard",
  "render_network": false
}
```

### scRNA-seq

```text
params.json
```

Example analysis settings include:

```json
{
  "de_mode": "pseudobulk",
  "de_scope": "target_celltype",
  "pseudobulk_method": "limma_voom",
  "network_url": "http://127.0.0.1:8001",
  "literature_url": "http://127.0.0.1:8006",
  "render_network": false
}
```

### Perturb-seq

```text
params.giladi.json
```

Example cell-state / target settings include:

```json
{
  "species": "mmusculus",
  "nt_label": "CTRL",
  "cell_state_url": "http://127.0.0.1:8007",
  "state_column": "lineage_state",
  "perturbation_column": "target",
  "cell_id_column": "cell_id",
  "annotation_mode": "provided",
  "render_network": false
}
```

---

## Development notes

### Pipeline-first design

The Nextflow workflows are the reproducible execution layer. The dashboard is a user interface over those workflows rather than a separate analytical implementation.

### API health checks

Integrated workflows perform downstream API checks before expensive analysis. If a required service is unavailable, fix the service before rerunning rather than allowing the workflow to fail far downstream.

### Reproducibility

Where available, outputs include:

- parameters
- run-status metadata
- manifests
- input/context handoffs
- intermediate CSV/JSON artifacts
- serialized analysis objects
- final JSON reports

### Large network graphs

Very large networks are valid analytical outputs even when they are inconvenient to render interactively. Preserve the graph tables and rankings rather than treating rendering failure as analysis failure.

### Virtual environments

Do not commit local virtual environments.

A repository `.gitignore` should exclude common local artifacts such as:

```gitignore
venv/
.venv/
VENV~~1/
__pycache__/
*.pyc
.env
.Rhistory
.RData
.Rproj.user/
work/
```

Generated result directories may also be excluded when they are large or machine-specific, depending on the intended repository policy.

---

## Scientific-use note

JCAP Bioinformatics Platform is research software under active development.

Network diffusion scores, candidate scores, literature classifications, and LLM-generated interpretations are **prioritization and interpretation aids**. They should be evaluated alongside the underlying statistical results, experimental design, domain expertise, and independent validation.

The platform is not intended to replace experimental validation or serve as a clinical diagnostic system.


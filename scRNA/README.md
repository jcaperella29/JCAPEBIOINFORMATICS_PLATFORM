# JCAP scRNA-seq API/CLI Module

This converts the multi-species scRNA-seq Shiny workflow into a reusable backend module: one R CLI is the analysis engine, and FastAPI is the job wrapper that accepts uploads and returns a result bundle.

The original Shiny app had these capabilities:

- Human, mouse, fly, and zebrafish species selection
- Seurat object creation from counts + metadata CSVs
- PCA and UMAP
- Differential expression by condition and by cell type
- Condition-only DE filtering
- g:Profiler or Enrichr enrichment
- Term-gene edge exports
- Random Forest feature selection/classification
- Power analysis
- Table and plot downloads

This API/CLI version keeps those concepts but makes the module callable from a dashboard, another service, an agent, Docker, local shell, or HPC.

## Folder layout

```text
scrna_api_cli/
  api/main.py              # FastAPI wrapper
  cli/run_scrna.R          # Seurat analysis CLI; source of truth
  sample_data/             # Tiny smoke-test CSVs
  examples/run_cli.sh      # local CLI smoke test
  examples/run_api.sh      # API smoke test with curl
  Dockerfile.cli           # CLI container
  Dockerfile.api           # API container
  requirements.txt         # Python API deps
```

## Expected inputs

### Counts CSV

Genes are rows, cells are columns. First column is used as row names.

```csv
gene,cell_01,cell_02,cell_03
CXCL10,0,2,5
ISG15,1,1,8
```

### Metadata CSV

Cells are rows. The first column is used as row names. For the full pipeline, include:

- `stim`: condition/group label
- `cell_type`: cell type label

```csv
cell_id,stim,cell_type
cell_01,ctrl,T_cell
cell_02,stim,T_cell
```

## CLI usage

```bash
Rscript cli/run_scrna.R \
  --counts sample_data/demo_counts.csv \
  --metadata sample_data/demo_metadata.csv \
  --species hsapiens \
  --enrich-backend gprof \
  --gprof-sources GO:BP,GO:MF,GO:CC,KEGG,REAC \
  --outdir results/scrna_demo
```

Quick smoke test without web enrichment or classifier:

```bash
bash examples/run_cli.sh
```

Useful flags:

```bash
--skip-pca-umap
--skip-enrichment
--skip-classifier
--skip-power
--zip
```

## API usage

Start the API:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn api.main:app --reload --port 8000
```

Then run a job:

```bash
curl -X POST http://127.0.0.1:8000/run \
  -F "counts_file=@sample_data/demo_counts.csv" \
  -F "meta_file=@sample_data/demo_metadata.csv" \
  -F "species=hsapiens" \
  -F "enrich_backend=gprof" \
  -F "gprof_sources=GO:BP,GO:MF,GO:CC,KEGG,REAC"
```

The response includes:

```json
{
  "job_id": "...",
  "status": "complete",
  "download_url": "/jobs/{job_id}/download",
  "tables_url": "/jobs/{job_id}/tables"
}
```

## API endpoints

| Method | Route | Purpose |
|---|---|---|
| GET | `/health` | service/version check |
| POST | `/validate` | upload sanity check without running Seurat |
| POST | `/run` | run the full CLI pipeline and create a bundle |
| GET | `/jobs/{job_id}` | read job status and summary |
| GET | `/jobs/{job_id}/tables` | list CSV outputs |
| GET | `/jobs/{job_id}/tables/{filename}` | download one CSV table |
| GET | `/jobs/{job_id}/download` | download the full ZIP bundle |
| DELETE | `/jobs/{job_id}` | remove a job directory |

## Output bundle

Each run writes:

```text
results/
  run_summary.json
  run_status.json
  tables/
    pca_coordinates.csv
    umap_coordinates.csv
    condition_de_table.csv
    celltype_de_table.csv
    condition_only_de_table.csv
    condition_only_volcano_table.csv
    enrichment_all.csv
    enrichment_up.csv
    enrichment_down.csv
    enrichment_all_edges.csv
    enrichment_up_edges.csv
    enrichment_down_edges.csv
    rf_importance_table.csv
    rf_predictions_table.csv
    rf_metrics_table.csv
    power_table.csv
    power_curves_long.csv
  plots/
    pca_plot.html
    umap_plot.html
    condition_only_volcano.html
    enrichment_all_barplot.html
    enrichment_up_barplot.html
    enrichment_down_barplot.html
  objects/
    seurat_object.rds
```

Some files are only created when the necessary metadata/results exist and the corresponding `--skip-*` flag is not used.

## Docker

CLI container:

```bash
docker build -f Dockerfile.cli -t jcap-scrna-cli .
docker run --rm -v "$PWD":/work jcap-scrna-cli \
  --counts /work/sample_data/demo_counts.csv \
  --metadata /work/sample_data/demo_metadata.csv \
  --outdir /work/results/scrna_demo \
  --skip-enrichment --skip-classifier --skip-power
```

API container:

```bash
docker build -f Dockerfile.api -t jcap-scrna-api .
docker run --rm -p 8000:8000 -v "$PWD/jobs":/app/jobs jcap-scrna-api
```

## Notes for dashboard integration

Use `/run` as the module execution route. The unified dashboard can show:

- uploaded counts + metadata names
- `run_summary.json`
- tables from `/jobs/{job_id}/tables`
- plot HTML files from the ZIP bundle or a later `/plots/{filename}` route
- one-click full bundle download from `/jobs/{job_id}/download`

The R CLI remains the source of truth; the API should stay thin.

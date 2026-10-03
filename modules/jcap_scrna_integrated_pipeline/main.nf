nextflow.enable.dsl=2

params.counts=null
params.metadata=null
params.context_json=null
params.outdir='results'

params.scrna_url='http://127.0.0.1:8003'
params.network_url='http://127.0.0.1:8001'
params.llm_url='http://127.0.0.1:8004'
params.literature_url='http://127.0.0.1:8006'
params.cell_state_api='http://127.0.0.1:8007'

params.condition_col='stim'
params.celltype_col='cell_type'
params.annotation_col='seurat_annotations'
params.sample_col='orig.ident'
params.pair_col=null
params.de_scope='target_celltype'

params.condition_a='stim'
params.condition_b='ctrl'
params.target_celltype='CD14+ Monocytes'
params.llm_mode='full'

params.cell_state_column='seurat_annotations'
params.cell_state_annotation_mode='provided'


workflow {

    if(!params.counts || !params.metadata || !params.context_json) {
        error 'counts, metadata, and context_json are required'
    }

    c = Channel.value(file(params.counts, checkIfExists:true))
    m = Channel.value(file(params.metadata, checkIfExists:true))
    x = Channel.value(file(params.context_json, checkIfExists:true))

    cell_state_runner = Channel.value(
        file("${projectDir}/scripts/run_cell_state_api.py", checkIfExists:true)
    )

    context_merger = Channel.value(
        file("${projectDir}/scripts/merge_cell_state_context.py", checkIfExists:true)
    )

    CHECK_APIS()

    CELL_STATE_ANALYSIS(
        m,
        cell_state_runner,
        CHECK_APIS.out.ready
    )

    MERGE_CELL_STATE_CONTEXT(
        x,
        CELL_STATE_ANALYSIS.out.handoff,
        context_merger
    )

    RUN_SCRNA(
        c,
        CELL_STATE_ANALYSIS.out.metadata,
        CHECK_APIS.out.ready
    )

    RUN_NETWORK(
        RUN_SCRNA.out.dir,
        MERGE_CELL_STATE_CONTEXT.out.context,
        CHECK_APIS.out.ready
    )

    RUN_LLM(
        RUN_SCRNA.out.dir,
        MERGE_CELL_STATE_CONTEXT.out.context,
        CHECK_APIS.out.ready
    )

    MAKE_LITERATURE_HANDOFF(
        RUN_LLM.out.dir,
        MERGE_CELL_STATE_CONTEXT.out.context
    )

    RUN_LITERATURE_REVIEW(
        RUN_SCRNA.out.dir,
        RUN_NETWORK.out.dir,
        MAKE_LITERATURE_HANDOFF.out.dir,
        MERGE_CELL_STATE_CONTEXT.out.context
    )

    ASSEMBLE_FINAL_REPORT(
        RUN_SCRNA.out.dir,
        RUN_NETWORK.out.dir,
        RUN_LLM.out.dir,
        RUN_LITERATURE_REVIEW.out.dir,
        MERGE_CELL_STATE_CONTEXT.out.context
    )

    BUNDLE(
        RUN_SCRNA.out.dir,
        RUN_NETWORK.out.dir,
        RUN_LLM.out.dir,
        RUN_LITERATURE_REVIEW.out.dir,
        ASSEMBLE_FINAL_REPORT.out.dir,
        CELL_STATE_ANALYSIS.out.metadata,
        MERGE_CELL_STATE_CONTEXT.out.context
    )
}


process CHECK_APIS {
    output:
    path 'ready.json', emit:ready

    script:
    """
    python3 '${projectDir}/scripts/check.py' \
      '${params.scrna_url}' \
      '${params.network_url}' \
      '${params.llm_url}' \
      '${params.literature_url}' \
      '${params.cell_state_api}' \
      > ready.json
    """
}


process CELL_STATE_ANALYSIS {
    tag "scRNA cell-state / lineage context"

    publishDir "${params.outdir}/cell_state", mode:'copy', overwrite:true

    input:
    path metadata
    path runner
    path ready

    output:
    path 'metadata_stateaware.csv', emit:metadata
    path 'cell_state_handoff.json', emit:handoff

    script:
    """
    python3 '${runner}' \
      --metadata '${metadata}' \
      --api-url '${params.cell_state_api}' \
      --assay scrna \
      --cell-id-column cell_id \
      --state-column '${params.cell_state_column}' \
      --condition-column '${params.condition_col}' \
      --annotation-mode '${params.cell_state_annotation_mode}' \
      --out-metadata metadata_stateaware.csv \
      --out-handoff cell_state_handoff.json
    """
}


process MERGE_CELL_STATE_CONTEXT {
    tag "Merge scRNA cell-state context"

    input:
    path context
    path handoff
    path merger

    output:
    path 'context_stateaware.json', emit:context

    script:
    """
    python3 '${merger}' \
      --context '${context}' \
      --cell-state-handoff '${handoff}' \
      --out context_stateaware.json
    """
}


process RUN_SCRNA {
    input:
    path c
    path m
    path z

    output:
    path '01_scrna', emit:dir

    script:
    """
    python3 '${projectDir}/scripts/scrna.py' \
      --api '${params.scrna_url}' \
      --counts '${c}' \
      --meta '${m}' \
      --condition-col '${params.condition_col}' \
      --celltype-col '${params.celltype_col}' \
      --annotation-col '${params.annotation_col}' \
      --sample-col '${params.sample_col}' \
      --pair-col '${params.pair_col}' \
      --de-scope '${params.de_scope}' \
      --a '${params.condition_a}' \
      --b '${params.condition_b}' \
      --target-celltype '${params.target_celltype}' \
      --out 01_scrna
    """
}


process RUN_NETWORK {
    input:
    path s
    path context
    path z

    output:
    path '02_network', emit:dir

    script:
    """
    python3 '${projectDir}/scripts/network.py' \
      --api '${params.network_url}' \
      --src '${s}' \
      --context '${context}' \
      --out 02_network
    """
}


process RUN_LLM {
    input:
    path s
    path x
    path z

    output:
    path '03_llm_triage', emit:dir

    script:
    """
    python3 '${projectDir}/scripts/llm.py' \
      --api '${params.llm_url}' \
      --src '${s}' \
      --context '${x}' \
      --mode '${params.llm_mode}' \
      --out 03_llm_triage
    """
}


process MAKE_LITERATURE_HANDOFF {
    input:
    path l
    path x

    output:
    path 'literature_handoff', emit:dir

    script:
    """
    mkdir -p literature_handoff

    python3 '${projectDir}/scripts/make_literature_handoff.py' \
      --llm-dir '${l}' \
      --context '${x}' \
      --out literature_handoff
    """
}


process RUN_LITERATURE_REVIEW {
    input:
    path s
    path n
    path l
    path x

    output:
    path '04_literature_review', emit:dir

    script:
    """
    python3 '${projectDir}/scripts/run_literature_api.py' \
      --api-url '${params.literature_url}' \
      --omics-dir '${s}' \
      --network-dir '${n}' \
      --llm-dir '${l}' \
      --context-json '${x}' \
      --outdir 04_literature_review
    """
}


process ASSEMBLE_FINAL_REPORT {
    tag "Assemble PI-facing integrated scRNA report"

    publishDir params.outdir, mode:'copy', overwrite:true

    input:
    path s
    path n
    path l
    path r
    path x

    output:
    path 'final_report', emit:dir

    script:
    """
    mkdir -p final_report

    python3 '${projectDir}/scripts/assemble_scrna_report.py' \
      --scrna-dir '${s}' \
      --network-dir '${n}' \
      --llm-dir '${l}' \
      --literature-dir '${r}' \
      --context-json '${x}' \
      --outdir final_report
    """
}


process BUNDLE {
    publishDir params.outdir, mode:'copy', overwrite:true

    input:
    path a
    path b
    path c
    path d
    path e
    path m
    path x

    output:
    path 'jcap_scrna_integrated_results'
    path 'jcap_scrna_integrated_results.zip'

    script:
    """
    python3 '${projectDir}/scripts/bundle.py' '${a}' '${b}' '${c}' '${m}' '${x}'

    rm -rf jcap_scrna_integrated_results/04_literature_review
    cp -RL '${d}' jcap_scrna_integrated_results/04_literature_review
    cp -RL '${e}' jcap_scrna_integrated_results/final_report

    mkdir -p jcap_scrna_integrated_results/cell_state
    cp '${m}' jcap_scrna_integrated_results/cell_state/metadata_stateaware.csv
    cp '${x}' jcap_scrna_integrated_results/cell_state/context_stateaware.json

    python3 - <<'PY'
from pathlib import Path
import zipfile

root = Path('jcap_scrna_integrated_results')
out = Path('jcap_scrna_integrated_results.zip')

if out.exists():
    out.unlink()

with zipfile.ZipFile(out, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
    for p in sorted(root.rglob('*')):
        if p.is_file():
            zf.write(p, p.as_posix())
PY
    """
}

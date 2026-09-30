nextflow.enable.dsl=2

params.counts=null; params.metadata=null; params.context_json=null; params.outdir='results'
params.scrna_url='http://127.0.0.1:8003'; params.network_url='http://127.0.0.1:8001'; params.llm_url='http://127.0.0.1:8004'; params.literature_url='http://127.0.0.1:8006'
params.condition_col='stim'; params.celltype_col='cell_type'; params.annotation_col='seurat_annotations'; params.sample_col='orig.ident'
params.condition_a='STIM'; params.condition_b='CTRL'; params.target_celltype='CD14+ Monocytes'; params.llm_mode='full'

workflow {

    if(!params.counts || !params.metadata || !params.context_json) {
        error 'counts, metadata, and context_json are required'
    }

    c=Channel.value(file(params.counts,checkIfExists:true))
    m=Channel.value(file(params.metadata,checkIfExists:true))
    x=Channel.value(file(params.context_json,checkIfExists:true))

    CHECK_APIS()
    RUN_SCRNA(c,m,CHECK_APIS.out.ready)
    RUN_NETWORK(RUN_SCRNA.out.dir,CHECK_APIS.out.ready)
    RUN_LLM(RUN_SCRNA.out.dir,x,CHECK_APIS.out.ready)

    RUN_LITERATURE_REVIEW(
        RUN_SCRNA.out.dir,
        RUN_NETWORK.out.dir,
        RUN_LLM.out.dir,
        x
    )

    BUNDLE(
        RUN_SCRNA.out.dir,
        RUN_NETWORK.out.dir,
        RUN_LLM.out.dir,
        RUN_LITERATURE_REVIEW.out.dir,
        m,
        x
    )
}

process CHECK_APIS {
    output:
    path 'ready.json', emit:ready

    script:
    """python3 '${projectDir}/scripts/check.py' '${params.scrna_url}' '${params.network_url}' '${params.llm_url}' > ready.json"""
}

process RUN_SCRNA {
    input:
    path c
    path m
    path z

    output:
    path '01_scrna', emit:dir

    script:
    """python3 '${projectDir}/scripts/scrna.py' --api '${params.scrna_url}' --counts '${c}' --meta '${m}' --condition-col '${params.condition_col}' --celltype-col '${params.celltype_col}' --annotation-col '${params.annotation_col}' --sample-col '${params.sample_col}' --pair-col '${params.pair_col}' --de-scope '${params.de_scope}' --a '${params.condition_a}' --b '${params.condition_b}' --target-celltype '${params.target_celltype}' --out 01_scrna"""
}

process RUN_NETWORK {
    input:
    path s
    path z

    output:
    path '02_network', emit:dir

    script:
    """python3 '${projectDir}/scripts/network.py' --api '${params.network_url}' --src '${s}' --out 02_network"""
}

process RUN_LLM {
    input:
    path s
    path x
    path z

    output:
    path '03_llm_triage', emit:dir

    script:
    """python3 '${projectDir}/scripts/llm.py' --api '${params.llm_url}' --src '${s}' --context '${x}' --mode '${params.llm_mode}' --out 03_llm_triage"""
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
    # cache-bust: scrna-network-literature-handoffs-v2
    python3 '${projectDir}/scripts/run_literature_api.py' \
      --api-url '${params.literature_url}' \
      --omics-dir '${s}' \
      --network-dir '${n}' \
      --llm-dir '${l}' \
      --context-json '${x}' \
      --outdir 04_literature_review
    """
}

process BUNDLE {
    publishDir params.outdir, mode:'copy', overwrite:true

    input:
    path a
    path b
    path c
    path d
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

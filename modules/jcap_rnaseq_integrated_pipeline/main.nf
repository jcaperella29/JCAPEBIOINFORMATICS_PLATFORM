nextflow.enable.dsl=2

params.counts              = null
params.phenotype           = null
params.phenotype_column    = null
params.context_json        = null
params.outdir              = "results"

params.rnaseq_url          = "http://127.0.0.1:8000"
params.network_url         = "http://127.0.0.1:8001"
params.llm_url             = "http://127.0.0.1:8004"
params.literature_url      = "http://127.0.0.1:8006"

params.species             = "hsapiens"
params.contrast            = ""
params.reference_group     = ""
params.fdr                 = 0.05
params.logfc_cutoff        = 1.0
params.enrich_backend      = "auto"
params.gprof_sources       = "GO:BP,GO:MF,GO:CC,KEGG,REAC"
params.enrichr_db          = "GO_Biological_Process_2023"
params.classifier_validation = "cv"
params.classifier_model    = "auto"

params.network_ranking_mode = "balanced"
params.network_projection_method = "jaccard"
params.network_top_n       = 50
params.network_candidate_top_n = 30

params.llm_mode            = "full"
params.make_pdf            = true
params.poll_seconds        = 5
params.timeout_minutes     = 240

workflow {

    if (!params.counts || !params.phenotype || !params.phenotype_column || !params.context_json) {
        error """
Required parameters:
  --counts PATH
  --phenotype PATH
  --phenotype_column NAME
  --context_json PATH

Example:
  nextflow run main.nf -params-file params.example.json
"""
    }

    counts_ch = Channel.value(file(params.counts, checkIfExists: true))
    phenotype_ch = Channel.value(file(params.phenotype, checkIfExists: true))
    context_ch = Channel.value(file(params.context_json, checkIfExists: true))

    CHECK_APIS()

    RUN_RNASEQ(counts_ch, phenotype_ch, CHECK_APIS.out.ready)
    RUN_NETWORK(RUN_RNASEQ.out.rnaseq_dir, CHECK_APIS.out.ready)
    RUN_LLM_TRIAGE(RUN_RNASEQ.out.rnaseq_dir, context_ch, CHECK_APIS.out.ready)

    RUN_LITERATURE_REVIEW(
        RUN_RNASEQ.out.rnaseq_dir,
        RUN_NETWORK.out.network_dir,
        RUN_LLM_TRIAGE.out.llm_dir,
        context_ch,
        CHECK_APIS.out.ready
    )

    BUILD_FINAL_REPORT(
        RUN_RNASEQ.out.rnaseq_dir,
        RUN_NETWORK.out.network_dir,
        RUN_LLM_TRIAGE.out.llm_dir,
        RUN_LITERATURE_REVIEW.out.literature_dir,
        context_ch
    )

    ASSEMBLE_BUNDLE(
        RUN_RNASEQ.out.rnaseq_dir,
        RUN_NETWORK.out.network_dir,
        RUN_LLM_TRIAGE.out.llm_dir,
        RUN_LITERATURE_REVIEW.out.literature_dir,
        BUILD_FINAL_REPORT.out.final_report_dir,
        counts_ch,
        phenotype_ch,
        context_ch
    )
}

process CHECK_APIS {
    tag "API health checks"

    output:
    path "apis_ready.txt", emit: ready

    script:
    """
    python3 ${projectDir}/scripts/check_apis.py \
      --rnaseq-url '${params.rnaseq_url}' \
      --network-url '${params.network_url}' \
      --llm-url '${params.llm_url}' \
      > apis_ready.txt
    """
}

process RUN_RNASEQ {
    tag "RNA-seq analysis"

    input:
    path counts
    path phenotype
    path api_ready

    output:
    path "01_rnaseq", emit: rnaseq_dir

    script:
    def contrastArg = params.contrast ? "--contrast '${params.contrast}'" : ""
    def refArg = params.reference_group ? "--reference-group '${params.reference_group}'" : ""
    """
    python3 ${projectDir}/scripts/run_rnaseq_api.py \
      --api-url '${params.rnaseq_url}' \
      --counts '${counts}' \
      --phenotype '${phenotype}' \
      --phenotype-column '${params.phenotype_column}' \
      --species '${params.species}' \
      --fdr '${params.fdr}' \
      --logfc-cutoff '${params.logfc_cutoff}' \
      --enrich-backend '${params.enrich_backend}' \
      --gprof-sources '${params.gprof_sources}' \
      --enrichr-db '${params.enrichr_db}' \
      --classifier-validation '${params.classifier_validation}' \
      --classifier-model '${params.classifier_model}' \
      --poll-seconds '${params.poll_seconds}' \
      --timeout-minutes '${params.timeout_minutes}' \
      ${contrastArg} \
      ${refArg} \
      --outdir 01_rnaseq
    """
}

process RUN_NETWORK {
    tag "Network analysis"

    input:
    path rnaseq_dir
    path api_ready

    output:
    path "02_network", emit: network_dir

    script:
    """
    python3 ${projectDir}/scripts/run_network_api.py \
      --api-url '${params.network_url}' \
      --rnaseq-dir '${rnaseq_dir}' \
      --ranking-mode '${params.network_ranking_mode}' \
      --projection-method '${params.network_projection_method}' \
      --top-n '${params.network_top_n}' \
      --candidate-top-n '${params.network_candidate_top_n}' \
      --outdir 02_network
    """
}

process RUN_LLM_TRIAGE {
    tag "LLM enrichment triage"

    input:
    path rnaseq_dir
    path context_json
    path api_ready

    output:
    path "03_llm_triage", emit: llm_dir

    script:
    """
    python3 ${projectDir}/scripts/run_llm_api.py \
      --api-url '${params.llm_url}' \
      --rnaseq-dir '${rnaseq_dir}' \
      --context-json '${context_json}' \
      --mode '${params.llm_mode}' \
      --make-pdf '${params.make_pdf}' \
      --outdir 03_llm_triage
    """
}

process RUN_LITERATURE_REVIEW {
    tag "Literature review + final synthesis"

    input:
    path rnaseq_dir
    path network_dir
    path llm_dir
    path context_json
    path api_ready

    output:
    path "04_literature_review", emit: literature_dir

    script:
    """
    echo "literature synthesis prompt v2" > /dev/null
    python3 ${projectDir}/scripts/run_literature_api.py \
      --api-url '${params.literature_url}' \
      --rnaseq-dir '${rnaseq_dir}' \
      --network-dir '${network_dir}' \
      --llm-dir '${llm_dir}' \
      --context-json '${context_json}' \
      --outdir 04_literature_review
    """
}


process BUILD_FINAL_REPORT {
    tag "Integrated RNA-seq final report"

    input:
    path rnaseq_dir
    path network_dir
    path llm_dir
    path literature_dir
    path context_json

    output:
    path "05_final_report", emit: final_report_dir

    script:
    """
    mkdir -p 05_final_report

    python3 ${projectDir}/scripts/assemble_rnaseq_report.py \
      --rnaseq-dir '${rnaseq_dir}' \
      --network-dir '${network_dir}' \
      --llm-dir '${llm_dir}' \
      --literature-dir '${literature_dir}' \
      --context-json '${context_json}' \
      --outdir 05_final_report
    """
}


process ASSEMBLE_BUNDLE {
    tag "Final results bundle"
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path rnaseq_dir
    path network_dir
    path llm_dir
    path literature_dir
    path final_report_dir
    path counts
    path phenotype
    path context_json

    output:
    path "jcap_integrated_results"
    path "jcap_integrated_results.zip", emit: final_zip

    script:
    """
    python3 ${projectDir}/scripts/assemble_bundle.py \
      --rnaseq-dir '${rnaseq_dir}' \
      --network-dir '${network_dir}' \
      --llm-dir '${llm_dir}' \
      --literature-dir '04_literature_review' \
      --counts '${counts}' \
      --phenotype '${phenotype}' \
      --context-json '${context_json}' \
      --outdir jcap_integrated_results \
      --zip-path jcap_integrated_results.zip

    cp '${final_report_dir}/rnaseq_final_report.md' jcap_integrated_results/
    cp '${final_report_dir}/rnaseq_final_report.json' jcap_integrated_results/

    rm -f jcap_integrated_results.zip
    python3 - <<'PYZIP'
import shutil
shutil.make_archive(
    "jcap_integrated_results",
    "zip",
    root_dir="jcap_integrated_results"
)
PYZIP
    """
}

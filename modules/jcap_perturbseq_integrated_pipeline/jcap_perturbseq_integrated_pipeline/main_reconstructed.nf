nextflow.enable.dsl=2

params.counts       = null
params.cells        = null
params.genes        = null
params.metadata     = null
params.context_json = null
params.targets      = null
params.outdir       = 'results/jcap_perturbseq_integrated_results'

params.network_url    = 'http://127.0.0.1:8001'
params.llm_url        = 'http://127.0.0.1:8004'
params.literature_url = 'http://127.0.0.1:8006'

params.species         = 'mmusculus'
params.nt_label        = 'CTRL'
params.neighbors       = 20
params.max_pcs         = 40
params.min_group_cells = 3
params.min_de_genes    = 3
params.iter_num        = 20

params.de_mode                   = 'auto'
params.pseudobulk_method         = 'limma_voom'
params.sample_col                = 'replicate'
params.pseudobulk_min_replicates = 2
params.de_fdr                    = 0.05
params.logfc_threshold           = 0
params.min_pct                   = 0.1
params.test_use                  = 'wilcox'

params.enrich_backend = 'auto'
params.gprof_sources  = 'GO:BP,GO:MF,GO:CC,REAC,KEGG'
params.enrichr_dbs    = 'GO_Biological_Process_2021,GO_Molecular_Function_2021,GO_Cellular_Component_2021,Reactome_2022,KEGG_2021_Human'

params.network_ranking_mode      = 'balanced'
params.network_projection_method = 'jaccard'
params.llm_mode                  = 'full'

params.max_findings               = 20
params.max_queries_per_finding    = 4
params.max_articles_per_query     = 5
params.max_articles_for_synthesis = 20


process CHECK_APIS {
    tag 'Perturb-seq downstream API health checks'

    output:
    path 'apis_ready.txt', emit: ready

    script:
    """
    set -euo pipefail

    python3 - <<'PY' > apis_ready.txt
import urllib.request

urls = {
    "network": "${params.network_url}",
    "llm": "${params.llm_url}",
    "literature": "${params.literature_url}",
}

for name, base in urls.items():
    ok = False
    last = None
    for suffix in ("/health", "/docs", "/openapi.json", ""):
        try:
            with urllib.request.urlopen(base.rstrip("/") + suffix, timeout=10) as r:
                if 200 <= r.status < 500:
                    print(f"{name}: reachable ({r.status}) at {base.rstrip('/') + suffix}")
                    ok = True
                    break
        except Exception as e:
            last = e
    if not ok:
        raise SystemExit(f"{name} API is not reachable at {base}: {last}")
PY
    """
}


process CLASSIFY_MIXSCAPE {
    tag 'Whole-dataset Mixscape classification'

    cpus 4
    memory '16 GB'
    time '12h'

    input:
    path counts
    path cells
    path genes
    path metadata
    path api_ready

    output:
    path '00_classification', emit: classification_dir
    path '00_classification/mixscape_classified.rds', emit: classified_rds

    script:
    """
    Rscript '${projectDir}/R/crispr_mixscape_giladi_classify.R'       --counts '${counts}'       --cells '${cells}'       --genes '${genes}'       --metadata '${metadata}'       --outdir 00_classification       --species '${params.species}'       --nt_label '${params.nt_label}'       --neighbors '${params.neighbors}'       --max_pcs '${params.max_pcs}'       --min_group_cells '${params.min_group_cells}'       --min_de_genes '${params.min_de_genes}'       --iter_num '${params.iter_num}'
    """
}


process SELECT_AVAILABLE_TARGETS {
    tag 'Select Mixscape-supported perturbations'

    input:
    path classification_dir

    output:
    path 'available_targets.txt', emit: available
    path 'skipped_targets.tsv', emit: skipped
    path 'target_selection.tsv', emit: selection

    script:
    """
    Rscript '${projectDir}/scripts/select_available_targets.R'       --classification-dir '${classification_dir}'       --targets '${params.targets}'       --available available_targets.txt       --skipped skipped_targets.tsv       --selection target_selection.tsv
    """
}


process RUN_TARGET_MIXSCAPE {
    tag { "Perturb-seq target: ${target}" }

    cpus 4
    memory '16 GB'
    time '12h'

    input:
    tuple val(target), path(classified_rds)

    output:
    tuple val(target), path('01_mixscape'), emit: results

    script:
    """
    Rscript '${projectDir}/R/crispr_mixscape_from_rds.R'       --rds '${classified_rds}'       --outdir 01_mixscape       --ko_label '${target}'       --species '${params.species}'       --nt_label '${params.nt_label}'       --de_mode '${params.de_mode}'       --pseudobulk_method '${params.pseudobulk_method}'       --sample_col '${params.sample_col}'       --pseudobulk_min_replicates '${params.pseudobulk_min_replicates}'       --de_fdr '${params.de_fdr}'       --logfc_threshold '${params.logfc_threshold}'       --min_pct '${params.min_pct}'       --test_use '${params.test_use}'       --enrich_backend '${params.enrich_backend}'       --gprof_sources '${params.gprof_sources}'       --enrichr_dbs '${params.enrichr_dbs}'
    """
}


process RUN_NETWORK {
    tag { "Network: ${target}" }

    input:
    tuple val(target), path(mixscape_dir)
    path context_json
    path api_ready

    output:
    tuple val(target), path('02_network'), emit: results

    script:
    """
    python3 '${projectDir}/scripts/network.py'       --api '${params.network_url}'       --src '${mixscape_dir}'       --context '${context_json}'       --ranking-mode '${params.network_ranking_mode}'       --projection-method '${params.network_projection_method}'       --out 02_network
    """
}


process RUN_LLM {
    tag { "LLM triage: ${target}" }

    input:
    tuple val(target), path(mixscape_dir)
    path context_json
    path api_ready

    output:
    tuple val(target), path('03_llm_triage'), emit: results

    script:
    """
    python3 '${projectDir}/scripts/llm.py'       --api '${params.llm_url}'       --src '${mixscape_dir}'       --context '${context_json}'       --mode '${params.llm_mode}'       --out 03_llm_triage
    """
}


process RUN_LITERATURE {
    tag { "Literature synthesis: ${target}" }

    input:
    tuple val(target), path(mixscape_dir), path(network_dir), path(llm_dir)
    path context_json
    path api_ready

    output:
    tuple val(target), path('04_literature_review'), emit: results

    script:
    """
    python3 '${projectDir}/scripts/run_literature_api.py'       --api-url '${params.literature_url}'       --omics-dir '${mixscape_dir}'       --network-dir '${network_dir}'       --llm-dir '${llm_dir}'       --context-json '${context_json}'       --outdir 04_literature_review       --max-findings '${params.max_findings}'       --max-queries-per-finding '${params.max_queries_per_finding}'       --max-articles-per-query '${params.max_articles_per_query}'       --max-articles-for-synthesis '${params.max_articles_for_synthesis}'
    """
}


process ASSEMBLE_TARGET {
    tag { "Bundle target: ${target}" }

    publishDir { "${params.outdir}/targets/${target}" }, mode: 'copy', overwrite: true

    input:
    tuple val(target), path(mixscape_dir), path(network_dir), path(llm_dir), path(literature_dir)

    output:
    tuple val(target), path("target_${target}"), emit: bundle

    script:
    """
    set -euo pipefail

    mkdir -p 'target_${target}'
    cp -RL '${mixscape_dir}'   'target_${target}/01_mixscape'
    cp -RL '${network_dir}'    'target_${target}/02_network'
    cp -RL '${llm_dir}'        'target_${target}/03_llm_triage'
    cp -RL '${literature_dir}' 'target_${target}/04_literature_review'
    """
}


process PUBLISH_SELECTION {
    tag 'Publish target selection'

    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path selection
    path skipped
    path classification_dir

    output:
    path 'target_selection.tsv'
    path 'skipped_targets.tsv'
    path '00_classification'

    script:
    """
    cp '${selection}' target_selection.tsv
    cp '${skipped}' skipped_targets.tsv
    cp -RL '${classification_dir}' 00_classification
    """
}


workflow {

    if (!params.counts || !params.cells || !params.genes || !params.metadata ||
        !params.context_json || !params.targets) {
        error """
Required parameters:
  --counts PATH
  --cells PATH
  --genes PATH
  --metadata PATH
  --context_json PATH
  --targets CEBPA,CSF1R,FCGR3,...

Example:
  nextflow run main.nf -params-file params.giladi.json
"""
    }

    counts_ch   = Channel.value(file(params.counts, checkIfExists: true))
    cells_ch    = Channel.value(file(params.cells, checkIfExists: true))
    genes_ch    = Channel.value(file(params.genes, checkIfExists: true))
    metadata_ch = Channel.value(file(params.metadata, checkIfExists: true))
    context_ch  = Channel.value(file(params.context_json, checkIfExists: true))

    CHECK_APIS()
    api_ready_ch = CHECK_APIS.out.ready.collect()

    CLASSIFY_MIXSCAPE(
        counts_ch,
        cells_ch,
        genes_ch,
        metadata_ch,
        api_ready_ch
    )

    SELECT_AVAILABLE_TARGETS(
        CLASSIFY_MIXSCAPE.out.classification_dir
    )

    PUBLISH_SELECTION(
        SELECT_AVAILABLE_TARGETS.out.selection,
        SELECT_AVAILABLE_TARGETS.out.skipped,
        CLASSIFY_MIXSCAPE.out.classification_dir
    )

    targets_ch = SELECT_AVAILABLE_TARGETS.out.available
        .splitText()
        .map { it.trim() }
        .filter { it }

    target_jobs = targets_ch.combine(CLASSIFY_MIXSCAPE.out.classified_rds)

    RUN_TARGET_MIXSCAPE(target_jobs)

    RUN_NETWORK(
        RUN_TARGET_MIXSCAPE.out.results,
        context_ch,
        api_ready_ch
    )

    RUN_LLM(
        RUN_TARGET_MIXSCAPE.out.results,
        context_ch,
        api_ready_ch
    )

    mix_net = RUN_TARGET_MIXSCAPE.out.results
        .join(RUN_NETWORK.out.results, by: 0)

    mix_net_llm = mix_net
        .join(RUN_LLM.out.results, by: 0)

    RUN_LITERATURE(
        mix_net_llm,
        context_ch,
        api_ready_ch
    )

    complete = mix_net_llm
        .join(RUN_LITERATURE.out.results, by: 0)

    ASSEMBLE_TARGET(complete)
}

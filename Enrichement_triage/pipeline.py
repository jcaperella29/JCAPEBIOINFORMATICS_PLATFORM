import pandas as pd
from triage import triage_enrichment_table
from program_summarizer import summarize_programs

from reasoner import gpt5_reason_simple
from playbook_retriever import retrieve_playbook_context
from pubmed_client import fetch_pubmed_context
from input_validation import validate_enrichment_df
from schemas import build_run_metadata, normalize_claims


def run_enrichment_pipeline(
    df: pd.DataFrame,
    *,
    phenotype: str,
    context: dict,
):
    validation = validate_enrichment_df(df)
    if not validation.get("ok"):
        raise ValueError(validation.get("error") or "Input validation failed.")

    # 1) stats + biofit + gene overlap clustering
    tri = triage_enrichment_table(
        df,
        phenotype=phenotype,
        context=context,
    )

    # 2) collapse into biological programs
    programs = summarize_programs(
        tri["rows"],
        phenotype=phenotype,
    )

    # 3) retrieve literature context from PubMed / NCBI
    pubmed_context = fetch_pubmed_context(
        phenotype=phenotype,
        context=context,
        triage=tri,
        programs=programs,
    )

    # 4) Retrieve only the relevant playbook guidance for this specific run.
    # This avoids injecting full markdown playbooks into every GPT prompt.
    playbook_context = retrieve_playbook_context(
        phenotype=phenotype,
        context=context,
        triage=tri,
        programs=programs,
    )

    # 5) GPT-5 + compact retrieved playbook guidance + PubMed evidence.
    gpt = gpt5_reason_simple(
        phenotype=phenotype,
        context=context,
        triage=tri,
        programs=programs,
        pubmed_context=pubmed_context,
        playbook_context=playbook_context,
    )

    claims = normalize_claims(gpt.get("parsed", {}) or {})

    metadata = build_run_metadata(
        df_columns=list(df.columns),
        n_rows=len(df),
        phenotype=phenotype,
        context=context,
        model_name=gpt.get("model", "gpt-5"),
    )

    return {
        "metadata": metadata,
        "input_validation": validation,
        "triage": tri,
        "programs": programs,
        "pubmed": pubmed_context,
        "playbook_context": playbook_context,
        "gpt": gpt,
        "claims": claims,
    }

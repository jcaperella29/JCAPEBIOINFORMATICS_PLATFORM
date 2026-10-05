#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def jload(path):
    p = Path(path)
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else {}


def csvload(path):
    p = Path(path)
    if not p.exists():
        return []
    with p.open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def find(root, name):
    hits = list(Path(root).rglob(name))
    return hits[0] if hits else None


def find_csv(root, exact=(), contains=()):
    root = Path(root)
    for name in exact:
        p = find(root, name)
        if p:
            return p
    for p in sorted(root.rglob('*.csv')):
        low = p.name.lower()
        if any(token.lower() in low for token in contains):
            return p
    return None


def first(row, keys):
    for k in keys:
        v = row.get(k)
        if v not in (None, ''):
            return str(v).strip()
    return None


def number(v):
    try:
        return float(v)
    except Exception:
        return None


def finding_labels(obj):
    labels = {}
    def walk(x):
        if isinstance(x, dict):
            fid = x.get('finding_id') or x.get('id')
            if isinstance(fid, str) and fid.startswith('finding_'):
                for key in ('claim','finding','finding_text','question','summary','description','text','title','hypothesis'):
                    val = x.get(key)
                    if isinstance(val, str) and val.strip():
                        labels[fid] = val.strip()
                        break
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(obj)
    return labels


def literature_markdown(request_obj, response_obj):
    assessments = response_obj.get('assessments') or []
    if not isinstance(assessments, list) or not assessments:
        return ['_No literature assessments were returned._']

    labels = finding_labels(request_obj)
    grouped = {}
    for a in assessments:
        if isinstance(a, dict):
            grouped.setdefault(str(a.get('finding_id') or 'unassigned'), []).append(a)

    def fkey(fid):
        try:
            return int(fid.rsplit('_', 1)[-1])
        except Exception:
            return 10**9

    def pkey(a):
        rank = {'direct':0,'support':1,'supporting':1,'indirect':2,'context':3,'exclude':9}
        level = str(a.get('level') or '').lower()
        try:
            rel = float(a.get('relevance_score'))
        except Exception:
            rel = 0.0
        return (rank.get(level, 5), -rel)

    out = []
    for fid in sorted(grouped, key=fkey):
        rows = grouped[fid]
        out += [f"### {labels.get(fid) or fid.replace('_',' ').title()}", '']
        if fid in labels:
            out += [f"*Literature target: `{fid}`*", '']

        counts = {}
        for row in rows:
            level = str(row.get('level') or 'unclassified')
            counts[level] = counts.get(level, 0) + 1
        count_text = ', '.join(f'**{k}**: {v}' for k,v in sorted(counts.items()))
        out += [f"**Evidence screen:** {len(rows)} paper{'s' if len(rows) != 1 else ''}" + (f" ({count_text})" if count_text else ''), '']

        for a in sorted(rows, key=pkey):
            paper = a.get('paper') or {}
            title = paper.get('title') or 'Untitled paper'
            meta = [str(x) for x in (paper.get('year'), paper.get('journal')) if x]
            if paper.get('pmid'):
                meta.append(f"PMID {paper['pmid']}")
            out.append('- **' + title + '**' + ((' — ' + '; '.join(meta)) if meta else ''))
            detail = f"  - Evidence level: **{a.get('level') or 'unclassified'}**"
            if a.get('relevance_score') is not None:
                try:
                    detail += f"; relevance score: **{float(a['relevance_score']):.2f}**"
                except Exception:
                    detail += f"; relevance score: **{a['relevance_score']}**"
            out.append(detail)
            for reason in a.get('reasons') or []:
                if reason:
                    out.append(f'  - {reason}')
            if a.get('caution'):
                out.append(f"  - Caution: {a['caution']}")
        out.append('')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rnaseq-dir', required=True)
    ap.add_argument('--network-dir', required=True)
    ap.add_argument('--llm-dir', required=True)
    ap.add_argument('--literature-dir', required=True)
    ap.add_argument('--context-json', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ctx = jload(a.context_json)
    context = ctx.get('context', ctx) if isinstance(ctx, dict) else {}
    if not isinstance(context, dict):
        context = {}

    run_summary_path = find(a.rnaseq_dir, 'run_summary.json')
    run_summary = jload(run_summary_path) if run_summary_path else {}

    de_path = find_csv(a.rnaseq_dir,
        exact=('DE_results.csv','differential_expression.csv','de_table.csv'),
        contains=('de_results_','limma_all_results_','differential'))
    de_rows = csvload(de_path) if de_path else []

    enrich_path = find_csv(a.rnaseq_dir, exact=('enrichment_all.csv',), contains=('enrichment_all',))
    enrich_rows = csvload(enrich_path) if enrich_path else []

    network_path = (find(a.network_dir, 'consensus_candidates.csv') or
                    find(a.network_dir, 'projection_top_candidates.csv') or
                    find(a.network_dir, 'bipartite_top_candidates.csv'))
    network_rows = csvload(network_path) if network_path else []

    claims_path = find(a.llm_dir, 'claims.csv')
    claims_rows = csvload(claims_path) if claims_path else []
    programs_path = find(a.llm_dir, 'programs.csv')
    programs_rows = csvload(programs_path) if programs_path else []
    llm_report_path = find(a.llm_dir, 'report.md')
    llm_report = llm_report_path.read_text(encoding='utf-8') if llm_report_path else ''

    lit_request_path = find(a.literature_dir, 'literature_review_request.json')
    lit_request = jload(lit_request_path) if lit_request_path else {}
    lit_path = find(a.literature_dir, 'literature_review_response.json')
    lit = jload(lit_path) if lit_path else {}

    dataset = context.get('dataset') or ctx.get('dataset') or 'bulk RNA-seq dataset'
    organism = context.get('organism') or ctx.get('organism') or ''
    assay = context.get('assay') or ctx.get('assay') or 'bulk RNA-seq'
    study = context.get('study_context') or context.get('phenotype') or ''
    tissue = context.get('cell_type') or context.get('tissue') or ''
    perturbation = context.get('perturbation') or context.get('treatment') or ''

    claims = []
    for row in claims_rows:
        val = first(row, ['claim','claim_text','summary','finding','description','interpretation','rationale'])
        if val and val not in claims:
            claims.append(val)

    programs = []
    for row in programs_rows:
        name = first(row, ['program','pathway','term','name','label'])
        summary = first(row, ['summary','description','interpretation','rationale'])
        if name or summary:
            programs.append((name, summary))

    def gene(row):
        return first(row, ['gene','Gene','symbol','Symbol','feature','Feature','Ensembl_IDs','ensembl_gene_id'])
    def fc(row):
        for k in ('logFC','avg_log2FC','log2FoldChange'):
            x = number(row.get(k))
            if x is not None:
                return x
    def padj(row):
        for k in ('adj.P.Val','padj','FDR','p_val_adj'):
            x = number(row.get(k))
            if x is not None:
                return x

    de_top = []
    for row in sorted(de_rows, key=lambda r: padj(r) if padj(r) is not None else 1e99):
        g = gene(row)
        if g:
            de_top.append((g, fc(row), padj(row)))
        if len(de_top) >= 15:
            break

    enrich_top = []
    for row in enrich_rows[:15]:
        term = first(row, ['Term','term','pathway','Pathway','Description','name'])
        p = None
        for k in ('Adjusted.P.value','adjusted_pvalue','padj','FDR','p.adjust','qvalue'):
            p = number(row.get(k))
            if p is not None:
                break
        source = first(row, ['Source','source','database','DB'])
        if term:
            enrich_top.append((term, p, source))

    network_top = []
    for row in network_rows[:10]:
        name = first(row, ['gene','Gene','candidate','Candidate','node','Node','name','Name','symbol','Symbol'])
        if not name:
            name = next((str(v) for v in row.values() if v not in (None,'')), 'unknown')
        score = first(row, ['consensus_score','score','Score','rank_score','diffusion_score','projection_score'])
        network_top.append((name, score))

    assessments = lit.get('assessments') or []
    n_samples = run_summary.get('matched_samples') or run_summary.get('n_samples')
    tested = run_summary.get('tested_genes')
    sig = run_summary.get('significant_genes')

    md = ['# Integrated Bulk RNA-seq Analysis', '', '## Executive Summary', '']
    md.append(f'This analysis integrates differential expression, enrichment, network prioritization, LLM biological triage, and literature review for **{dataset}**.')
    if n_samples:
        md.append(f'The RNA-seq analysis included **{int(n_samples):,} matched samples**.')
    if tested is not None:
        s = f'The differential-expression model tested **{int(tested):,} genes**'
        s += f' and identified **{int(sig):,} significant genes**.' if sig is not None else '.'
        md.append(s)
    if claims:
        md.append(f'LLM triage produced **{len(claims)} structured biological claims**.')
    if isinstance(assessments, list):
        md.append(f'The literature stage produced **{len(assessments)} evidence assessments**.')

    md += ['', '### PI-facing Biological Synthesis', '']
    if claims:
        md += [f'- {x}' for x in claims[:8]]
    elif llm_report:
        md.append(llm_report[:5000])
    else:
        md.append('_No LLM synthesis was available._')

    md += ['', '## Experimental Context', '', f'- **Dataset:** {dataset}', f'- **Organism:** {organism}', f'- **Assay:** {assay}']
    if study: md.append(f'- **Study context:** {study}')
    if tissue: md.append(f'- **Cell/tissue context:** {tissue}')
    if perturbation: md.append(f'- **Perturbation/treatment:** {perturbation}')

    md += ['', '## Differential Expression', '']
    if de_top:
        md += ['| Gene | log2FC | adjusted p |', '|---|---:|---:|']
        for g, f, p in de_top:
            md.append(f"| {g} | {'' if f is None else f'{f:.3f}'} | {'' if p is None else f'{p:.3g}'} |")
    else:
        md.append('_No differential-expression table available._')

    md += ['', '## Enrichment / Biological Programs', '']
    if enrich_top:
        md += ['| Term | adjusted p | Source |', '|---|---:|---|']
        for term, p, source in enrich_top:
            md.append(f"| {term} | {'' if p is None else f'{p:.3g}'} | {source or ''} |")
    elif programs:
        for name, summary in programs[:12]:
            md.append(f"- **{name}** — {summary}" if name and summary else f"- **{name}**" if name else f'- {summary}')
    else:
        md.append('_No enrichment/program table available._')

    md += ['', '## LLM Triage', '']
    if claims:
        md += [f'- {x}' for x in claims[:12]]
    elif llm_report:
        md.append(llm_report[:5000])
    else:
        md.append('_No LLM triage output available._')

    md += ['', '## Network Prioritization', '']
    if network_top:
        md += ['| Candidate | Score |', '|---|---:|']
        md += [f'| {name} | {score or ""} |' for name, score in network_top]
    else:
        md.append('_No network candidate table available._')

    md += ['', '## Literature Evidence', '']
    md += literature_markdown(lit_request, lit)
    md += ['', '## Interpretation Notes', '',
           '- Differential-expression statistics are the primary quantitative evidence from the RNA-seq experiment.',
           '- Enrichment and network outputs prioritize programs and candidates; they do not independently establish causality.',
           '- LLM triage is an interpretation layer and should be checked against the numerical source tables.',
           '- Literature retrieval contextualizes the experiment; primary papers should be reviewed directly before publication or experimental decisions.',
           '', '## Supporting Material', '',
           '- `01_rnaseq/` — bulk RNA-seq DE, enrichment, classifier, QC, and plots',
           '- `02_network/` — network and candidate prioritization',
           '- `03_llm_triage/` — LLM claims, programs, and report',
           '- `04_literature_review/` — literature request and evidence response',
           '', 'Numerical source tables remain the authoritative record for individual statistics.']

    report = '\n'.join(md).rstrip() + '\n'
    md_path = out / 'rnaseq_final_report.md'
    md_path.write_text(report, encoding='utf-8')

    payload = {
        'dataset': dataset,
        'context': ctx,
        'run_summary': run_summary,
        'leading_de': [{'gene':g,'log2fc':f,'adjusted_p':p} for g,f,p in de_top],
        'leading_enrichment': [{'term':t,'adjusted_p':p,'source':s} for t,p,s in enrich_top],
        'network_candidates': network_top,
        'llm_claims': claims,
        'llm_programs': programs,
        'literature_summary_counts': lit.get('summary_counts') or {},
        'literature_assessments': assessments,
    }
    json_path = out / 'rnaseq_final_report.json'
    json_path.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')

    print(md_path)
    print(json_path)


if __name__ == '__main__':
    main()

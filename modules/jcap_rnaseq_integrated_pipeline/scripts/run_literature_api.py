#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, zipfile
from pathlib import Path
import requests

KNOWN = {
    'rnaseq_literature_handoff.json',
    'network_literature_handoff.json',
    'literature_handoff.json',
    'enrichment_triage_literature_handoff.json',
}


def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))


def find_handoffs(*roots):
    found = []
    for root in map(Path, roots):
        if not root.exists(): continue
        candidates = [root] if root.is_file() else [p for n in KNOWN for p in root.rglob(n)] + list(root.rglob('*literature*handoff*.json'))
        for p in candidates:
            p = p.resolve()
            if p not in found: found.append(p)
    return found


def recover_network(root):
    root = Path(root)
    for zp in sorted(root.rglob('*.zip')):
        try:
            with zipfile.ZipFile(zp) as z:
                for n in z.namelist():
                    if Path(n).name == 'network_literature_handoff.json':
                        out = root/'network_literature_handoff.json'
                        out.write_bytes(z.read(n))
                        return out
        except Exception:
            pass


def first(row, names):
    for name in names:
        v = row.get(name)
        if v not in (None, ''): return str(v).strip()
    return None


def fallback_targets(rnaseq_dir, llm_dir, limit=20):
    targets = []
    llm = Path(llm_dir)
    claims = next(llm.rglob('claims.csv'), None)
    programs = next(llm.rglob('programs.csv'), None)
    for path, kind in [(claims, 'claim'), (programs, 'program')]:
        if not path: continue
        with path.open(newline='', encoding='utf-8-sig') as fh:
            for row in csv.DictReader(fh):
                summary = first(row, ['claim','claim_text','summary','finding','description','interpretation','rationale'])
                gene = first(row, ['gene','Gene','symbol','Symbol','target_gene'])
                pathway = first(row, ['pathway','program','term','name','label'])
                if summary or gene or pathway:
                    item = {}
                    if summary: item['summary'] = summary
                    if gene: item['genes'] = [gene]
                    if pathway: item['pathways'] = [pathway]
                    targets.append(item)
                if len(targets) >= limit: return targets

    # Last fallback: leading DE genes from common bulk RNA-seq tables.
    csvs = list(Path(rnaseq_dir).rglob('*.csv'))
    de = next((p for p in csvs if p.name.lower() in {'de_results.csv','differential_expression.csv','de_table.csv'}), None) \
         or next((p for p in csvs if 'differential' in p.name.lower() or p.name.lower().startswith('de_')), None)
    if de:
        with de.open(newline='', encoding='utf-8-sig') as fh:
            for row in csv.DictReader(fh):
                gene = first(row, ['gene','Gene','symbol','Symbol','feature','Feature'])
                if gene:
                    targets.append({'genes':[gene], 'summary':f'{gene} identified in bulk RNA-seq differential expression.'})
                if len(targets) >= limit: break
    return targets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--api-url', required=True)
    ap.add_argument('--rnaseq-dir', required=True)
    ap.add_argument('--network-dir', required=True)
    ap.add_argument('--llm-dir', required=True)
    ap.add_argument('--context-json')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--max-findings', type=int, default=20)
    ap.add_argument('--max-articles-per-query', type=int, default=5)
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    context = load(a.context_json) if a.context_json and Path(a.context_json).exists() else {}

    network_skip = next(Path(a.network_dir).rglob('network_skipped.txt'), None)
    llm_skip = next(Path(a.llm_dir).rglob('llm_skipped.txt'), None)
    if network_skip and llm_skip:
        reason = network_skip.read_text().strip() + ' | ' + llm_skip.read_text().strip()
        (out/'literature_skipped.txt').write_text(f'Literature review skipped because both network and LLM stages were skipped. {reason}\n')
        print('Literature review skipped: no interpretive upstream outputs available.')
        return

    handoff_paths = find_handoffs(a.rnaseq_dir, a.network_dir, a.llm_dir)
    if not any('network_literature_handoff' in p.name for p in handoff_paths):
        recover_network(a.network_dir)
        handoff_paths = find_handoffs(a.rnaseq_dir, a.network_dir, a.llm_dir)

    handoffs = []
    for p in handoff_paths:
        try:
            d = load(p)
            if isinstance(d, dict): handoffs.append(d)
        except Exception as e:
            print(f'WARNING: could not load {p}: {e}')

    merged_context = context.get('context', context) if isinstance(context, dict) else {}
    if not isinstance(merged_context, dict): merged_context = {}
    for h in handoffs:
        hc = h.get('context') or {}
        if isinstance(hc, dict):
            for k,v in hc.items():
                if v not in (None, ''): merged_context[k] = v

    findings = []
    counter = 0
    for h in handoffs:
        existing = h.get('findings')
        if isinstance(existing, list):
            for f in existing:
                if isinstance(f, dict) and f.get('summary'):
                    counter += 1
                    x = dict(f)
                    x.setdefault('finding_id', f'finding_{counter:03d}')
                    x.setdefault('assay_type', merged_context.get('assay') or 'bulk RNA-seq')
                    findings.append(x)
        for item in h.get('literature_targets') or []:
            if len(findings) >= a.max_findings: break
            counter += 1
            if isinstance(item, str):
                item = {'summary': item}
            if not isinstance(item, dict): continue
            genes = item.get('genes') or ([item['gene']] if item.get('gene') else [])
            pathways = item.get('pathways') or ([item['pathway']] if item.get('pathway') else [])
            phenotypes = item.get('phenotypes') or ([item['phenotype']] if item.get('phenotype') else [])
            findings.append({
                'finding_id': f'finding_{counter:03d}',
                'assay_type': merged_context.get('assay') or 'bulk RNA-seq',
                'target': item.get('target') or item.get('perturbation'),
                'genes': genes,
                'pathways': pathways,
                'phenotypes': phenotypes,
                'organism': merged_context.get('organism'),
                'tissue_or_cell_type': merged_context.get('cell_type') or merged_context.get('tissue'),
                'disease_or_context': merged_context.get('phenotype') or merged_context.get('study_context'),
                'summary': item.get('summary') or item.get('finding') or item.get('description') or 'RNA-seq literature-review target.',
            })

    if not findings:
        for item in fallback_targets(a.rnaseq_dir, a.llm_dir, a.max_findings):
            counter += 1
            findings.append({
                'finding_id': f'finding_{counter:03d}',
                'assay_type': merged_context.get('assay') or 'bulk RNA-seq',
                'target': merged_context.get('perturbation'),
                'genes': item.get('genes', []),
                'pathways': item.get('pathways', []),
                'phenotypes': item.get('phenotypes', []),
                'organism': merged_context.get('organism'),
                'tissue_or_cell_type': merged_context.get('cell_type') or merged_context.get('tissue'),
                'disease_or_context': merged_context.get('phenotype') or merged_context.get('study_context'),
                'summary': item.get('summary', 'Bulk RNA-seq result selected for literature review.'),
            })

    if not findings:
        (out/'literature_skipped.txt').write_text('Literature review skipped: no literature targets could be derived from RNA-seq, network, or LLM outputs.\n')
        print('Literature review skipped: no findings available.')
        return

    pipeline_handoff = {
        'run_id': 'rnaseq_integrated',
        'analysis_type': 'Bulk RNA-seq integrated analysis',
        'findings': findings[:a.max_findings],
        'quantitative_results': {'rnaseq_handoffs':[h for h in handoffs if 'rnaseq' in str(h.get('source_api','')).lower()]},
        'network_results': {'handoffs':[h for h in handoffs if 'network' in str(h.get('source_api','')).lower()]},
        'llm_triage': {'handoffs':[h for h in handoffs if any(x in str(h.get('source_api','')).lower() for x in ('llm','triage'))]},
        'metadata': {'context': merged_context, 'source_handoff_count': len(handoffs)},
    }
    payload = {'handoff': pipeline_handoff, 'max_papers_per_finding': a.max_articles_per_query}
    (out/'literature_review_request.json').write_text(json.dumps(payload, indent=2) + '\n')

    r = requests.post(
        a.api_url.rstrip('/') + '/review',
        json=payload,
        timeout=1800
    )

    if r.status_code != 200:
        raise SystemExit(
            f'Literature API HTTP {r.status_code}\n'
            f'{r.text[:4000]}'
        )

    result = r.json()

    (out / 'literature_review_response.json').write_text(
        json.dumps(result, indent=2) + '\n'
    )

    print(
        f"Literature review complete: "
        f"{len(result.get('assessments', []))} evidence assessments"
    )


if __name__ == '__main__': main()

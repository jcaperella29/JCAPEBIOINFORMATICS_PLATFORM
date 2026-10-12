#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,sys,zipfile
from pathlib import Path
import requests
KNOWN={'rnaseq_literature_handoff.json','network_literature_handoff.json','literature_handoff.json','enrichment_triage_literature_handoff.json','scrna_literature_handoff.json','crispr_mixscape_literature_handoff.json'}
def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def find_handoffs(*roots):
    found=[]
    for root in map(Path,roots):
        if not root.exists(): continue
        cs=[root] if root.is_file() else [p for n in KNOWN for p in root.rglob(n)]+list(root.rglob('*literature*handoff*.json'))
        for p in cs:
            p=p.resolve()
            if p not in found: found.append(p)
    return found
def recover_network(root):
    root=Path(root)
    for zp in sorted(root.rglob('*.zip')):
        try:
            with zipfile.ZipFile(zp) as z:
                for n in z.namelist():
                    if Path(n).name=='network_literature_handoff.json':
                        out=root/'network_literature_handoff.json'; out.write_bytes(z.read(n)); return out
        except Exception: pass
def augment(h,c):
    if not c:return h
    existing=h.get('context') if isinstance(h.get('context'),dict) else {}; cand=c.get('context') if isinstance(c.get('context'),dict) else c
    for k in ['phenotype','organism','assay','tissue','cell_type','perturbation','selected_perturbation']:
        if not existing.get(k) and cand.get(k) is not None: existing[k]=cand[k]
    h['context']=existing; return h
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--api-url',required=True);ap.add_argument('--omics-dir',required=True);ap.add_argument('--network-dir',required=True);ap.add_argument('--llm-dir',required=True);ap.add_argument('--context-json');ap.add_argument('--outdir',required=True);ap.add_argument('--max-findings',type=int,default=20);ap.add_argument('--max-queries-per-finding',type=int,default=4);ap.add_argument('--max-articles-per-query',type=int,default=5);ap.add_argument('--max-articles-for-synthesis',type=int,default=20);a=ap.parse_args()
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True);ctx=load(a.context_json) if a.context_json and Path(a.context_json).exists() else None

    network_skip = next(Path(a.network_dir).rglob('network_skipped.txt'), None)
    llm_skip = next(Path(a.llm_dir).rglob('llm_skipped.txt'), None)

    if network_skip or llm_skip:
        reasons = []

        if network_skip:
            reasons.append(
                network_skip.read_text(encoding='utf-8').strip()
            )

        if llm_skip:
            reasons.append(
                llm_skip.read_text(encoding='utf-8').strip()
            )

        reason = ' | '.join(reasons)

        (out/'literature_skipped.txt').write_text(
            f'Literature review skipped because upstream analysis was skipped. {reason}\n',
            encoding='utf-8'
        )

        print(f'Literature review skipped: {reason}')
        return

    hp=find_handoffs(a.omics_dir,a.network_dir,a.llm_dir)
    if not any('network_literature_handoff' in p.name for p in hp): recover_network(a.network_dir); hp=find_handoffs(a.omics_dir,a.network_dir,a.llm_dir)
    hand=[]
    for p in hp:
        try:
            d=load(p)
            if isinstance(d,dict): hand.append(augment(d,ctx)); print(f'+ {p} [source_api={d.get("source_api","unknown")}]')
        except Exception as e: print(f'WARNING {p}: {e}',file=sys.stderr)
    if not hand:
        (out/'literature_skipped.txt').write_text(
            'Literature review skipped: no JCAP literature handoff JSON files were produced.\n',
            encoding='utf-8'
        )
        print('Literature review skipped: no JCAP literature handoff JSON files were produced.')
        return
    # Canonical perturbation provenance: trust Mixscape's actual DE target.
    canonical = None

    for h in hand:
        if h.get('source_api') == 'crispr_mixscape':
            ctx = h.get('context') or {}
            canonical = ctx.get('selected_perturbation')
            if canonical:
                break

    if canonical:
        canonical = str(canonical).strip()
        print(f'Canonical perturbation from Mixscape: {canonical}')

        for h in hand:
            ctx = h.get('context') or {}
            ctx['selected_perturbation'] = canonical
            h['context'] = ctx

    # ---- Adapt accumulated JCAP handoffs to current Literature API ----

    # Merge useful context from the handoffs.
    merged_context = {}
    for h in hand:
        ctx = h.get('context') or {}
        if isinstance(ctx, dict):
            merged_context.update({k: v for k, v in ctx.items() if v not in (None, '')})

    if canonical:
        merged_context['selected_perturbation'] = canonical

    assay = (
        merged_context.get('assay')
        or 'Perturb-seq / scRNA-seq'
    )
    organism = merged_context.get('organism')
    tissue_or_cell_type = (
        merged_context.get('cell_type')
        or merged_context.get('tissue')
    )
    disease_or_context = (
        merged_context.get('phenotype')
        or merged_context.get('study_context')
        or merged_context.get('perturbation')
    )

    findings = []

    # Preserve any already-structured StudyFinding-like entries.
    for h in hand:
        existing = h.get('findings')
        if isinstance(existing, list):
            for f in existing:
                if not isinstance(f, dict):
                    continue
                if 'finding_id' in f and 'assay_type' in f and 'summary' in f:
                    findings.append(f)

    # Convert explicit literature targets from the handoffs.
    finding_counter = len(findings)

    for h in hand:
        targets = h.get('literature_targets') or []

        if not isinstance(targets, list):
            continue

        for item in targets:
            finding_counter += 1

            genes = []
            pathways = []
            phenotypes = []
            target = canonical
            summary = None

            if isinstance(item, str):
                genes = [item]
                summary = f'{item} identified as a literature-review target.'

            elif isinstance(item, dict):
                raw_gene = (
                    item.get('gene')
                    or item.get('symbol')
                    or item.get('target_gene')
                )

                if raw_gene:
                    genes = [str(raw_gene)]

                raw_genes = item.get('genes')
                if isinstance(raw_genes, list):
                    genes.extend(str(x) for x in raw_genes if x)

                raw_pathways = item.get('pathways') or item.get('pathway')
                if isinstance(raw_pathways, list):
                    pathways = [str(x) for x in raw_pathways if x]
                elif raw_pathways:
                    pathways = [str(raw_pathways)]

                raw_pheno = item.get('phenotypes') or item.get('phenotype')
                if isinstance(raw_pheno, list):
                    phenotypes = [str(x) for x in raw_pheno if x]
                elif raw_pheno:
                    phenotypes = [str(raw_pheno)]

                target = (
                    item.get('target')
                    or item.get('perturbation')
                    or canonical
                )

                summary = (
                    item.get('summary')
                    or item.get('finding')
                    or item.get('description')
                )

                if not summary:
                    bits = genes + pathways + phenotypes
                    summary = (
                        'Literature-review target: ' + ', '.join(bits)
                        if bits
                        else f'Literature-review finding for {target or "Perturb-seq experiment"}.'
                    )

            else:
                continue

            # Preserve order while removing duplicates.
            genes = list(dict.fromkeys(genes))
            pathways = list(dict.fromkeys(pathways))
            phenotypes = list(dict.fromkeys(phenotypes))

            findings.append({
                'finding_id': f'finding_{finding_counter:03d}',
                'assay_type': assay,
                'target': target,
                'genes': genes,
                'pathways': pathways,
                'phenotypes': phenotypes,
                'organism': organism,
                'tissue_or_cell_type': tissue_or_cell_type,
                'disease_or_context': disease_or_context,
                'summary': summary,
            })

    # Fallback: make one target-level finding if no explicit literature
    # targets were supplied.
    if not findings:
        top_genes = []

        for h in hand:
            de = h.get('top_perturbation_de') or []
            if not isinstance(de, list):
                continue

            for row in de[:10]:
                if isinstance(row, dict):
                    gene = (
                        row.get('gene')
                        or row.get('symbol')
                        or row.get('feature')
                    )
                    if gene:
                        top_genes.append(str(gene))

        top_genes = list(dict.fromkeys(top_genes))

        findings.append({
            'finding_id': 'finding_001',
            'assay_type': assay,
            'target': canonical,
            'genes': top_genes,
            'pathways': [],
            'phenotypes': [],
            'organism': organism,
            'tissue_or_cell_type': tissue_or_cell_type,
            'disease_or_context': disease_or_context,
            'summary': (
                f'Perturb-seq transcriptional response associated with {canonical}.'
                if canonical
                else 'Perturb-seq transcriptional response.'
            ),
        })

    crispr_handoffs = [
        h for h in hand
        if str(h.get('source_api', '')).lower() == 'crispr_mixscape'
    ]

    network_handoffs = [
        h for h in hand
        if 'network' in str(h.get('source_api', '')).lower()
    ]

    llm_handoffs = [
        h for h in hand
        if any(
            x in str(h.get('source_api', '')).lower()
            for x in ('llm', 'triage')
        )
    ]

    pipeline_handoff = {
        'run_id': f'perturbseq_{canonical or "unknown"}',
        'analysis_type': 'Perturb-seq / Mixscape integrated analysis',
        'findings': findings,
        'quantitative_results': {
            'crispr_mixscape_handoffs': crispr_handoffs
        },
        'network_results': {
            'handoffs': network_handoffs
        },
        'llm_triage': {
            'handoffs': llm_handoffs
        },
        'metadata': {
            'context': merged_context,
            'selected_perturbation': canonical,
            'source_handoff_count': len(hand),
        },
    }

    payload = {
        'handoff': pipeline_handoff,
        'max_papers_per_finding': a.max_articles_per_query,
    }

    (out/'literature_review_request.json').write_text(
        json.dumps(payload, indent=2),
        encoding='utf-8'
    )

    r = requests.post(
        a.api_url.rstrip('/') + '/review',
        json=payload,
        timeout=1800
    )

    if r.status_code != 200:
        raise SystemExit(
            f'Literature API HTTP {r.status_code}\n{r.text[:4000]}'
        )

    result = r.json()

    (out/'literature_review_response.json').write_text(
        json.dumps(result, indent=2),
        encoding='utf-8'
    )

    print(
        f"Literature review complete: "
        f"{len(result.get('assessments', []))} evidence assessments"
    )


if __name__=='__main__':main()

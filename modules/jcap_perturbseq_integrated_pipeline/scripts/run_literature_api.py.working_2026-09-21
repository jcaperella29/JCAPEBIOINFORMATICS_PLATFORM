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

    payload={'handoffs':hand,'max_findings':a.max_findings,'max_queries_per_finding':a.max_queries_per_finding,'max_articles_per_query':a.max_articles_per_query,'max_articles_for_synthesis':a.max_articles_for_synthesis,'use_pubmed':True,'use_europe_pmc':True,'use_openalex':True,'use_semantic_scholar':False,'run_gpt_synthesis':True}
    (out/'literature_review_request.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
    r=requests.post(a.api_url.rstrip('/')+'/export/report-bundle',json=payload,timeout=1800)
    if r.status_code!=200: raise SystemExit(f'Literature API HTTP {r.status_code}\n{r.text[:4000]}')
    zp=out/'jcap_literature_review.zip';zp.write_bytes(r.content)
    if not zipfile.is_zipfile(zp):raise SystemExit('Literature API did not return ZIP')
    with zipfile.ZipFile(zp) as z:z.extractall(out)
if __name__=='__main__':main()

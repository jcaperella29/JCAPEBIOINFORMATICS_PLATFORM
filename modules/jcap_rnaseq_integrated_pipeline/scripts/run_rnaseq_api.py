#!/usr/bin/env python3
from __future__ import annotations
import argparse, io, json, shutil, time, zipfile
from pathlib import Path
import requests


def as_file_tuple(path: Path):
    return (path.name, path.open('rb'), 'text/csv')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--api-url', required=True)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--phenotype', required=True)
    ap.add_argument('--phenotype-column', required=True)
    ap.add_argument('--species', default='hsapiens')
    ap.add_argument('--contrast', default='')
    ap.add_argument('--reference-group', default='')
    ap.add_argument('--fdr', type=float, default=0.05)
    ap.add_argument('--logfc-cutoff', type=float, default=1.0)
    ap.add_argument('--enrich-backend', default='auto')
    ap.add_argument('--gprof-sources', default='GO:BP,GO:MF,GO:CC,KEGG,REAC')
    ap.add_argument('--enrichr-db', default='GO_Biological_Process_2023')
    ap.add_argument('--classifier-validation', default='cv')
    ap.add_argument('--classifier-model', default='auto')
    ap.add_argument('--poll-seconds', type=int, default=5)
    ap.add_argument('--timeout-minutes', type=int, default=240)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    if a.contrast and a.reference_group:
        raise SystemExit('Use either --contrast or --reference-group, not both.')

    api = a.api_url.rstrip('/')
    counts = Path(a.counts)
    phenotype = Path(a.phenotype)
    out = Path(a.outdir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    form = {
        'phenotype_column': a.phenotype_column,
        'species': a.species,
        'fdr': str(a.fdr),
        'logfc_cutoff': str(a.logfc_cutoff),
        'enrich_backend': a.enrich_backend,
        'gprof_sources': a.gprof_sources,
        'enrichr_db': a.enrichr_db,
        'classifier_validation': a.classifier_validation,
        'classifier_model': a.classifier_model,
        'zip_outputs': 'true',
    }
    if a.contrast:
        form['contrast'] = a.contrast
    if a.reference_group:
        form['reference_group'] = a.reference_group

    print('Submitting RNA-seq job...')
    print(json.dumps(form, indent=2))

    with counts.open('rb') as fc, phenotype.open('rb') as fp:
        r = requests.post(
            f'{api}/rnaseq/jobs',
            files={
                'counts': (counts.name, fc, 'text/csv'),
                'phenotype': (phenotype.name, fp, 'text/csv'),
            },
            data=form,
            timeout=3600,
        )
    if r.status_code not in (200, 201, 202):
        raise RuntimeError(f'RNA-seq submission failed: HTTP {r.status_code}\n{r.text}')

    submit = r.json()
    (out / 'submission.json').write_text(json.dumps(submit, indent=2) + '\n')
    job_id = submit['job_id']
    print(f'Submitted RNA-seq job: {job_id}')

    deadline = time.time() + a.timeout_minutes * 60
    status = submit
    while True:
        r = requests.get(f'{api}/rnaseq/jobs/{job_id}', timeout=60)
        r.raise_for_status()
        status = r.json()
        (out / 'job_status.json').write_text(json.dumps(status, indent=2) + '\n')
        state = str(status.get('status', '')).lower()
        print(f'RNA-seq job status: {state}')
        if state in {'completed', 'complete', 'success', 'succeeded'}:
            break
        if state in {'failed', 'error', 'deleted'}:
            raise RuntimeError('RNA-seq job failed:\n' + json.dumps(status, indent=2))
        if time.time() >= deadline:
            raise TimeoutError(f'RNA-seq job timed out after {a.timeout_minutes} minutes')
        time.sleep(max(1, a.poll_seconds))

    url = f'{api}/rnaseq/jobs/{job_id}/download'
    print(f'Downloading RNA-seq bundle from {url}')
    r = requests.get(url, timeout=3600)
    r.raise_for_status()
    data = r.content
    bio = io.BytesIO(data)
    if not zipfile.is_zipfile(bio):
        fallback = out / 'download_response.bin'
        fallback.write_bytes(data)
        raise RuntimeError(f'RNA-seq download was not a ZIP bundle; saved response to {fallback}')

    bio.seek(0)
    with zipfile.ZipFile(bio) as zf:
        zf.extractall(out)

    print(f'RNA-seq results ready: {out}')


if __name__ == '__main__':
    main()

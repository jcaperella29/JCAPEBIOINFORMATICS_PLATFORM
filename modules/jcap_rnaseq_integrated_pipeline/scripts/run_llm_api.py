#!/usr/bin/env python3
import argparse, json, shutil, subprocess
from pathlib import Path


def bool_text(x):
    return 'true' if str(x).strip().lower() in {'1','true','yes','y','on'} else 'false'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--api-url', required=True)
    p.add_argument('--rnaseq-dir', required=True)
    p.add_argument('--context-json', required=True)
    p.add_argument('--mode', required=True)
    p.add_argument('--make-pdf', default='true')
    p.add_argument('--outdir', required=True)
    a = p.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    src = Path(a.rnaseq_dir)

    csvs = list(src.rglob('*.csv'))
    f = next((x for x in csvs if x.name.lower() == 'enrichment_all.csv'), None) \
        or next((x for x in csvs if 'enrichment' in x.name.lower() and 'edge' not in x.name.lower()), None)
    if not f:
        reason = 'LLM skipped: no usable enrichment CSV was produced by RNA-seq.'
        (out/'llm_skipped.txt').write_text(reason + '\n')
        print(reason); return

    c = json.loads(Path(a.context_json).read_text())
    phenotype = str(c.pop('phenotype', 'RNA-seq differential expression / enrichment response'))
    std = {k: str(c.pop(k, '')) for k in ['organism','assay','tissue','cell_type','perturbation','timepoint']}
    z = out/'llm_bundle.zip'
    cmd = ['curl','--fail-with-body','--silent','--show-error','-X','POST',a.api_url.rstrip('/')+'/analyze-bundle','-F',f'file=@{f}','--form-string',f'phenotype={phenotype}']
    for k,v in std.items(): cmd += ['--form-string', f'{k}={v}']
    cmd += ['--form-string', f'extra_context_json={json.dumps(c)}', '--form-string','input_format=auto', '--form-string',f'mode={a.mode}', '--form-string',f'make_pdf={bool_text(a.make_pdf)}', '-o', str(z)]
    subprocess.run(cmd, check=True)
    results = out/'results'; results.mkdir()
    shutil.unpack_archive(z, results)


if __name__ == '__main__': main()

#!/usr/bin/env python3
import argparse, json, shutil, zipfile
from pathlib import Path


def copy_any(src, dst):
    src, dst = Path(src), Path(dst)
    if src.is_dir(): shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rnaseq-dir', required=True)
    ap.add_argument('--network-dir', required=True)
    ap.add_argument('--llm-dir', required=True)
    ap.add_argument('--literature-dir', required=True)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--phenotype', required=True)
    ap.add_argument('--context-json', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--zip-path', required=True)
    a = ap.parse_args()

    out = Path(a.outdir)
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)

    copy_any(a.rnaseq_dir, out/'01_rnaseq')
    copy_any(a.network_dir, out/'02_network')
    copy_any(a.llm_dir, out/'03_llm_triage')
    copy_any(a.literature_dir, out/'04_literature_review')

    inputs = out/'inputs'; inputs.mkdir()
    counts_name = Path(a.counts).name
    phenotype_name = Path(a.phenotype).name
    copy_any(a.counts, inputs/counts_name)
    copy_any(a.phenotype, inputs/phenotype_name)
    copy_any(a.context_json, inputs/'context.json')

    manifest = {
        'pipeline': 'jcap_rnaseq_integrated_pipeline',
        'rnaseq_results': '01_rnaseq',
        'network_results': '02_network',
        'llm_triage': '03_llm_triage',
        'literature_review': '04_literature_review',
        'inputs': {
            'counts': f'inputs/{counts_name}',
            'phenotype': f'inputs/{phenotype_name}',
            'context': 'inputs/context.json',
        },
    }
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

    zp = Path(a.zip_path)
    if zp.exists(): zp.unlink()
    with zipfile.ZipFile(zp, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for p in out.rglob('*'):
            if p.is_file(): z.write(p, arcname=str(p.relative_to(out.parent)))
    print(f'Wrote {out}')
    print(f'Wrote {zp}')


if __name__ == '__main__': main()

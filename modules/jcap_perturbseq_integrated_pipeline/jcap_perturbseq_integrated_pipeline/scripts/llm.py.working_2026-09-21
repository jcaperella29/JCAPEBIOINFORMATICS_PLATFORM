#!/usr/bin/env python3
import argparse,json,subprocess,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--api',required=True);p.add_argument('--src',required=True);p.add_argument('--context',required=True);p.add_argument('--mode',required=True);p.add_argument('--out',required=True);a=p.parse_args();o=Path(a.out);o.mkdir()

src=Path(a.src)

network_skip=next(src.rglob('network_skipped.txt'),None)
if network_skip:
    reason=network_skip.read_text(encoding='utf-8').strip()
    (o/'llm_skipped.txt').write_text(
        f'LLM skipped because upstream network was skipped. {reason}\n',
        encoding='utf-8'
    )
    print(f'LLM skipped: {reason}')
    raise SystemExit(0)

cs=list(src.rglob('*.csv'))
f=next((x for x in cs if x.name.lower()=='enrichment_all.csv'),None) or next((x for x in cs if 'enrichment' in x.name.lower() and 'network_edges' not in x.name.lower()),None) or next((x for x in cs if 'network_edges' in x.name.lower()),None)
if not f:
    (o/'llm_skipped.txt').write_text(
        'LLM skipped: no usable network/enrichment CSV was produced.\n',
        encoding='utf-8'
    )
    print('LLM skipped: no usable network/enrichment CSV was produced.')
    raise SystemExit(0)
c=json.loads(Path(a.context).read_text()); ph=c.pop('phenotype','CRISPR perturbation transcriptional response'); std={k:str(c.pop(k,'')) for k in ['organism','assay','tissue','cell_type','perturbation','timepoint']}; fmt='long_edges' if 'network_edges' in f.name.lower() else 'auto'; z=o/'llm_bundle.zip'
cmd=['curl','--fail-with-body','--silent','--show-error','-X','POST',a.api.rstrip('/')+'/analyze-bundle','-F',f'file=@{f}','--form-string',f'phenotype={ph}']
for k,v in std.items():cmd+=['--form-string',f'{k}={v}']
cmd += ['--form-string',f'extra_context_json={json.dumps(c)}','--form-string',f'input_format={fmt}','--form-string',f'mode={a.mode}','--form-string','make_pdf=true','-o',str(z)]
subprocess.run(cmd,check=True); (o/'results').mkdir(); shutil.unpack_archive(z,o/'results')


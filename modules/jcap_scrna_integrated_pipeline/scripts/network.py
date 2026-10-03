#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, io, json, re, subprocess, urllib.request, zipfile
from pathlib import Path

def post_json(url,payload,timeout=900):
    req=urllib.request.Request(url,data=json.dumps(payload).encode(),headers={'Content-Type':'application/json'},method='POST')
    with urllib.request.urlopen(req,timeout=timeout) as r: return json.loads(r.read())
def write_json(p,x): Path(p).write_text(json.dumps(x,indent=2),encoding='utf-8')
def write_text(p,x): Path(p).write_text(x or '',encoding='utf-8')
def find_enrichment(root):
    cs=sorted(Path(root).rglob('*.csv'))
    for n in ['enrichment_all.csv','enrichment_up.csv','enrichment_down.csv']:
        h=[p for p in cs if p.name.lower()==n];
        if h:return h[0]
    h=[p for p in cs if 'enrichment' in p.name.lower() and 'edge' not in p.name.lower()]
    if h:return h[0]
    raise SystemExit('No enrichment CSV found.\n'+'\n'.join(map(str,cs)))
def make_edges(src,dst):
    with open(src,encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    if not rows: raise SystemExit('Enrichment table is empty')
    cols=rows[0].keys(); term=next((c for c in cols if c.lower() in {'term','pathway','description'}),None); genes=next((c for c in cols if c.lower() in {'genes','gene_list','core_enrichment'}),None); padj=next((c for c in cols if c.lower().replace('.','_') in {'adjusted_p_value','padj','fdr','qvalue'}),None)
    if not term or not genes:
        msg = rows[0].get('message', '') if rows else ''
        reason = msg or f'Need Term and Genes columns; found {list(cols)}'

        skip = Path(dst).parent / 'network_skipped.txt'
        skip.write_text(f'Network skipped: {reason}\n', encoding='utf-8')

        print(f'Network skipped: {reason}')
        return None
    out=[]
    for r in rows:
        for g in re.split(r'[,;/|]+',r.get(genes,'') or ''):
            g=g.strip()
            if g: out.append({'gene':g,'term':r.get(term,''),'adjusted_pvalue':r.get(padj,'') if padj else ''})
    with open(dst,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['gene','term','adjusted_pvalue']); w.writeheader(); w.writerows(out)
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--api',required=True);ap.add_argument('--src',required=True);ap.add_argument('--context',required=True);ap.add_argument('--ranking-mode',default='balanced');ap.add_argument('--projection-method',default='jaccard');ap.add_argument('--out',required=True);a=ap.parse_args()
    api=a.api.rstrip('/'); out=Path(a.out);out.mkdir(parents=True,exist_ok=True); src=find_enrichment(a.src); edges=out/'enrichment_all_edges.csv'
    made = make_edges(src,edges)
    write_text(out/'selected_input.txt',str(src))
    if made is None:
        write_json(out/'network_status.json', {
            'status': 'skipped',
            'reason': (out/'network_skipped.txt').read_text(encoding='utf-8').strip()
        })
        return
    opts={'apply_preset':False,'item_col':'gene','group_col':'term','weight_col':'adjusted_pvalue'}
    cmd=['curl','--fail-with-body','--silent','--show-error','-X','POST',f'{api}/network/build','-F',f'file=@{edges}','-F',f'options_json={json.dumps(opts)}']
    r=subprocess.run(cmd,text=True,capture_output=True); 
    if r.returncode: raise SystemExit(r.stdout+'\n'+r.stderr)
    built=json.loads(r.stdout); write_json(out/'network_build.json',built); graph=built['graph']; write_text(out/'main_nodes.csv',built.get('exports',{}).get('nodes_csv')); write_text(out/'main_edges.csv',built.get('exports',{}).get('edges_csv'))
    bip=post_json(f'{api}/diffusion/bipartite',{'graph':graph,'seed_node':None,'alpha':0.85,'ranking_mode':a.ranking_mode,'top_n':50,'candidate_top_n':30,'candidate_node_type':'group'}); write_json(out/'bipartite_diffusion.json',bip);write_text(out/'bipartite_diffusion_results.csv',bip.get('csv'));write_text(out/'bipartite_top_candidates.csv',bip.get('candidate_csv'))
    proj=post_json(f'{api}/projection/build',{'graph':graph,'method':a.projection_method,'return_figure':True,'show_labels':True});write_json(out/'projection_build.json',proj);write_text(out/'projection_nodes.csv',proj.get('exports',{}).get('nodes_csv'));write_text(out/'projection_edges.csv',proj.get('exports',{}).get('edges_csv'))
    pd=post_json(f'{api}/diffusion/projection',{'projection_graph':proj['projection_graph'],'seed_node':None,'alpha':0.85,'ranking_mode':a.ranking_mode,'top_n':50,'candidate_top_n':30});write_json(out/'projection_diffusion.json',pd);write_text(out/'projection_diffusion_results.csv',pd.get('csv'));write_text(out/'projection_top_candidates.csv',pd.get('candidate_csv'))
    con=post_json(f'{api}/consensus',{'bipartite_results':bip,'projection_results':pd,'top_n':30});write_json(out/'consensus.json',con);write_text(out/'consensus_candidates.csv',con.get('csv'))
    ctx=json.loads(Path(a.context).read_text())
    payload={'graph':graph,'projection_graph':proj['projection_graph'],'bipartite_results':bip,'projection_results':pd,'consensus_results':con,'mapped_columns':{'item_col':'gene','group_col':'term','weight_col':'adjusted_pvalue'},'settings':{'ranking_mode':a.ranking_mode,'projection_method':a.projection_method},'context':ctx.get('context',ctx),'input_source':str(src)}
    req=urllib.request.Request(f'{api}/export/report-bundle',data=json.dumps(payload).encode(),headers={'Content-Type':'application/json'},method='POST')
    with urllib.request.urlopen(req,timeout=1200) as q: blob=q.read()
    zp=out/'network_report_bundle.zip'; zp.write_bytes(blob)
    with zipfile.ZipFile(zp) as z:z.extractall(out)
if __name__=='__main__':main()

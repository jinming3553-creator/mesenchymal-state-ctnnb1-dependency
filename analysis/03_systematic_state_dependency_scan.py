#!/usr/bin/env python
from pathlib import Path
import sys, numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.stats_utils import residualize, pearson_pair, fixed_effect_meta, bh_fdr
P=ROOT/'data/processed'; R=ROOT/'results/source_data'; R.mkdir(parents=True,exist_ok=True)

# The exact 20-state list used in the submitted analysis is frozen here when the repository is finalized.
# If `config/systematic_states.txt` exists, it is used; otherwise all programmes are intersected and the
# user is asked to create the frozen list before a publication run.
state_file=ROOT/'config/systematic_states.txt'
if not state_file.exists():
    raise SystemExit('Create config/systematic_states.txt containing the exact 20 prespecified state names, one per line.')
states=[x.strip() for x in state_file.read_text().splitlines() if x.strip()]
if len(states)!=20: raise SystemExit(f'Expected 20 states, found {len(states)}')

bs=pd.read_csv(P/'broad_program_scores.csv',index_col=0); bg=pd.read_pickle(P/'broad_gene_effect.pkl'); bm=pd.read_csv(P/'broad_metadata.csv',index_col=0)
ss=pd.read_csv(P/'sanger_program_scores.csv',index_col=0); sg=pd.read_pickle(P/'sanger_gene_effect.pkl'); sm=pd.read_csv(P/'sanger_metadata.csv',index_col=0)
genes=sorted(set(bg.columns)&set(sg.columns))
rows=[]
for state in states:
    bx=residualize(bs[state],bm['lineage']); sx=residualize(ss[state],sm['lineage'])
    for gene in genes:
        br,bp,bn=pearson_pair(bx,residualize(bg[gene],bm['lineage']))
        sr,sp,sn=pearson_pair(sx,residualize(sg[gene],sm['lineage']))
        meta=fixed_effect_meta([br,sr],[bn,sn])
        rows.append([state,gene,br,bp,bn,sr,sp,sn,meta['r'],meta['p']])
out=pd.DataFrame(rows,columns=['State','Gene','Broad_r','Broad_P','Broad_N','Sanger_r','Sanger_P','Sanger_N','Meta_r','Meta_P'])
out['Meta_FDR']=bh_fdr(out['Meta_P'])
out['Same_direction']=(out['Broad_r'] * out['Sanger_r']) > 0
out.to_csv(R/'crosscohort_state_dependency_meta.csv',index=False)
print('Tests:',len(out),'FDR<0.05:',int((out.Meta_FDR<0.05).sum()))
print(out.sort_values('Meta_FDR').head(10).to_string(index=False))

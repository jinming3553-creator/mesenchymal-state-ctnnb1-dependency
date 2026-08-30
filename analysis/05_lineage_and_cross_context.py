#!/usr/bin/env python
from pathlib import Path
import sys, pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.stats_utils import pearson_pair, fixed_effect_meta, bh_fdr
P=ROOT/'data/processed'; R=ROOT/'results/source_data'; R.mkdir(parents=True,exist_ok=True)
rows=[]
for prefix in ['broad','sanger']:
    s=pd.read_csv(P/f'{prefix}_program_scores.csv',index_col=0); g=pd.read_pickle(P/f'{prefix}_gene_effect.pkl'); m=pd.read_csv(P/f'{prefix}_metadata.csv',index_col=0)
    d=pd.DataFrame({'EMT-I':s['EMT-I'],'CTNNB1':g['CTNNB1'],'lineage':m['lineage']}).dropna()
    for line,x in d.groupby('lineage'):
        r,p,n=pearson_pair(x['EMT-I'],x['CTNNB1']); rows.append([prefix,line,n,r,p])
pd.DataFrame(rows,columns=['Cohort','Lineage','N','r','P']).to_csv(R/'lineage_specific_EMTI_CTNNB1.csv',index=False)

# Secondary selected-context architecture: Broad Esophagus/Stomach vs Sanger COLO.
bs=pd.read_csv(P/'broad_program_scores.csv',index_col=0); bg=pd.read_pickle(P/'broad_gene_effect.pkl'); bm=pd.read_csv(P/'broad_metadata.csv',index_col=0)
ss=pd.read_csv(P/'sanger_program_scores.csv',index_col=0); sg=pd.read_pickle(P/'sanger_gene_effect.pkl'); sm=pd.read_csv(P/'sanger_metadata.csv',index_col=0)
programs=sorted(set(bs.columns)&set(ss.columns)); out=[]
bi=bm.index[bm.lineage.eq('Esophagus/Stomach')]; si=sm.index[sm.lineage.eq('COLO')]
for pr in programs:
    br,bp,bn=pearson_pair(bs.loc[bi,pr],bg.loc[bi,'CTNNB1']); sr,sp,sn=pearson_pair(ss.loc[si,pr],sg.loc[si,'CTNNB1'])
    mm=fixed_effect_meta([br,sr],[bn,sn]); out.append([pr,br,bp,sr,sp,mm['r'],mm['p']])
out=pd.DataFrame(out,columns=['Program','Broad_ES_r','Broad_ES_P','Sanger_COLO_r','Sanger_COLO_P','Meta_r','Meta_P'])
out['Meta_FDR']=bh_fdr(out.Meta_P); out['Same_direction']=(out.Broad_ES_r*out.Sanger_COLO_r)>0
out.to_csv(R/'cross_lineage_program_architecture.csv',index=False)

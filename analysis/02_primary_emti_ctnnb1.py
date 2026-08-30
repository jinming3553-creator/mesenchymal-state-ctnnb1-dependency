#!/usr/bin/env python
from pathlib import Path
import sys, numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.stats_utils import residualize, pearson_pair, spearman_pair, fisher_ci, fixed_effect_meta
OUT=ROOT/'data/processed'; RES=ROOT/'results/source_data'; RES.mkdir(parents=True,exist_ok=True)

def cohort(prefix,id_name):
    scores=pd.read_csv(OUT/f'{prefix}_program_scores.csv',index_col=0)
    ge=pd.read_pickle(OUT/f'{prefix}_gene_effect.pkl')
    meta=pd.read_csv(OUT/f'{prefix}_metadata.csv',index_col=0)
    d=pd.DataFrame({'EMT-I':scores['EMT-I'],'CTNNB1':ge['CTNNB1'],'lineage':meta['lineage']}).dropna()
    d['EMT-I adjusted']=residualize(d['EMT-I'],d['lineage'])
    d['CTNNB1 adjusted']=residualize(d['CTNNB1'],d['lineage'])
    r,p,n=pearson_pair(d['EMT-I adjusted'],d['CTNNB1 adjusted'])
    rs,ps,_=spearman_pair(d['EMT-I adjusted'],d['CTNNB1 adjusted'])
    lo,hi=fisher_ci(r,n)
    d.index.name=id_name; d.to_csv(RES/f'{prefix}_EMTI_CTNNB1_model_level.csv')
    return {'cohort':prefix,'n':n,'r':r,'p':p,'ci_low':lo,'ci_high':hi,'spearman_rho':rs,'spearman_p':ps}

b=cohort('broad','ModelID'); s=cohort('sanger','sample_ID')
m=fixed_effect_meta([b['r'],s['r']],[b['n'],s['n']])
pd.DataFrame([b,s,{"cohort":"fixed-effect meta","n":b['n']+s['n'],**m}]).to_csv(RES/'EMTI_CTNNB1_primary_summary.csv',index=False)
print(pd.DataFrame([b,s])); print('META',m)

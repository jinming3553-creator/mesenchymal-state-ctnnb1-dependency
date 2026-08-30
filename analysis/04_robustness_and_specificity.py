#!/usr/bin/env python
from pathlib import Path
import sys, numpy as np, pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.io_utils import load_genesets
from src.scoring import score_programs
from src.stats_utils import residualize, pearson_pair, fixed_effect_meta
P=ROOT/'data/processed'; RAW=ROOT/'data/raw'; R=ROOT/'results/source_data'; R.mkdir(parents=True,exist_ok=True)
STROMA=['COL1A2','SPARC','COL1A1','COL3A1','TAGLN','FN1','MGP','BGN','COL6A2','COL4A2','COL6A1','MMP2','DCN','COL4A1']
WNT=['CTNNB1','APC','TCF7L2','LGR4','MESD','PORCN','WLS','FZD5','LRP6']
gs=load_genesets(RAW/'geneset_table_updated.csv'); emti=gs['EMT-I']
rows=[]; loo=[]
for prefix in ['broad','sanger']:
    expr=pd.read_pickle(P/f'{prefix}_expression.pkl'); ge=pd.read_pickle(P/f'{prefix}_gene_effect.pkl'); meta=pd.read_csv(P/f'{prefix}_metadata.csv',index_col=0)
    signatures={'Full EMT-I':emti,'Stroma-light EMT-I':[g for g in emti if g not in STROMA]}
    sc,_=score_programs(expr,signatures,min_genes=5)
    for label in sc:
        r,p,n=pearson_pair(residualize(sc[label],meta.lineage),residualize(ge.CTNNB1,meta.lineage)); rows.append([prefix,label,n,r,p])
    for omit in [g for g in emti if g in expr.columns]:
        score,_=score_programs(expr,{'EMT-I':[g for g in emti if g!=omit]},min_genes=5)
        r,p,n=pearson_pair(residualize(score['EMT-I'],meta.lineage),residualize(ge.CTNNB1,meta.lineage)); loo.append([prefix,omit,n,r,p])
pd.DataFrame(rows,columns=['Cohort','Score','N','r','P']).to_csv(R/'stroma_light_EMTI.csv',index=False)
pd.DataFrame(loo,columns=['Cohort','Omitted_gene','N','r','P']).to_csv(R/'EMTI_leave_one_gene_out.csv',index=False)

# Representative WNT specificity using same primary lineage adjustment.
w=[]
for gene in WNT:
    vals=[]
    for prefix in ['broad','sanger']:
        scores=pd.read_csv(P/f'{prefix}_program_scores.csv',index_col=0); ge=pd.read_pickle(P/f'{prefix}_gene_effect.pkl'); meta=pd.read_csv(P/f'{prefix}_metadata.csv',index_col=0)
        if gene not in ge: vals.append((np.nan,np.nan,0)); continue
        vals.append(pearson_pair(residualize(scores['EMT-I'],meta.lineage),residualize(ge[gene],meta.lineage)))
    m=fixed_effect_meta([vals[0][0],vals[1][0]],[vals[0][2],vals[1][2]])
    w.append([gene,vals[0][0],vals[0][1],vals[1][0],vals[1][1],m['r'],m['p']])
pd.DataFrame(w,columns=['Gene','Broad_r','Broad_P','Sanger_r','Sanger_P','Meta_r','Meta_P']).to_csv(R/'WNT_specificity.csv',index=False)

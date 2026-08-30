#!/usr/bin/env python
from pathlib import Path
import sys, re, numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.io_utils import load_genesets
from src.stats_utils import pearson_pair, spearman_pair
import hdf5plugin  # registers HDF5 compression filters used by the public H5AD files
import anndata as ad

RAW=ROOT/'data/raw'; R=ROOT/'results/source_data'; R.mkdir(parents=True,exist_ok=True)
EMT=load_genesets(RAW/'geneset_table_updated.csv')['EMT-I']
EPI=['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT20','MUC1','KRT17','KRT23']
PAIRS={
'COAD_0143':('AG001_HCM-CSHL-0143-C20_gt.h5ad','AG002_HCM-CSHL-0143-C20_gt.h5ad','COAD'),
'COAD_0322':('AG001_HCM-CSHL-0322-C20_gt.h5ad','AG002_HCM-CSHL-0322-C20_gt.h5ad','COAD'),
'COAD_BROD0001':('AG004_HCM-BROD-0001-C18_gt.h5ad','AG005_HCM-BROD-0001-C18_gt.h5ad','COAD'),
'PAAD_BROD0110':('AG008_HCM-BROD-0110-C25_gt.h5ad','AG009_HCM-BROD-0110-C25_gt.h5ad','PAAD'),
'PAAD_0073':('AG008_HCM-CSHL-0073-C25_gt.h5ad','AG009_HCM-CSHL-0073-C25_gt.h5ad','PAAD'),
'PAAD_0089_primary':('AG010_HCM-CSHL-0089-C25_gt.h5ad','AG011_HCM-CSHL-0089-C25_gt.h5ad','PAAD'),
'PAAD_0089_meta':('AG012_HCM-CSHL-0089-C25_gt.h5ad','AG013_HCM-CSHL-0089-C25_gt.h5ad','PAAD'),
'PAAD_0078_primary':('AG010_HCM-CSHL-0078-C25_gt.h5ad','AG011_HCM-CSHL-0078-C25_gt.h5ad','PAAD'),
'PAAD_0078_meta':('AG012_HCM-CSHL-0078-C25_gt.h5ad','AG013_HCM-CSHL-0078-C25_gt.h5ad','PAAD'),
}

def pseudobulk_score(path, epithelial=False):
    a=ad.read_h5ad(path)
    genes=pd.Index(a.var_names.astype(str))
    mask=np.ones(a.n_obs,dtype=bool)
    if epithelial:
        epi=[g for g in EPI if g in genes]
        if not epi: raise ValueError('No epithelial markers found')
        X=a[:,epi].X
        mask=np.asarray((X>0).sum(axis=1)).ravel()>0
    aa=a[mask]
    sig=[g for g in EMT if g in genes]
    X=aa[:,sig].X
    counts=np.asarray(X.sum(axis=0)).ravel()
    total=float(np.asarray(aa.X.sum()).ravel()[0] if hasattr(np.asarray(aa.X.sum()),'ravel') else aa.X.sum())
    cpm=counts/total*1e6
    return float(np.log2(cpm+1).mean()), int(mask.sum()), len(sig)

rows=[]
for key,(tf,mf,line) in PAIRS.items():
    for mode in ['whole','epithelial_enriched']:
        epi=mode!='whole'
        ts,tn,ng=pseudobulk_score(RAW/tf,epi); ms,mn,_=pseudobulk_score(RAW/mf,epi)
        rows.append([key,line,mode,ts,ms,ms-ts,tn,mn,ng])
out=pd.DataFrame(rows,columns=['Pair','Lineage','Mode','Tumour_EMTI','Model_EMTI','Delta_model_minus_tumour','Tumour_nuclei','Model_nuclei','EMTI_genes'])
out.to_csv(R/'HCMI_EMTI_pair_scores.csv',index=False)

stats=[]
for mode,d in out.groupby('Mode'):
    d=d.copy(); d['Tumour_adj']=d.Tumour_EMTI-d.groupby('Lineage').Tumour_EMTI.transform('mean'); d['Model_adj']=d.Model_EMTI-d.groupby('Lineage').Model_EMTI.transform('mean')
    r,p,n=pearson_pair(d.Tumour_adj,d.Model_adj); rs,ps,_=spearman_pair(d.Tumour_adj,d.Model_adj)
    stats.append([mode,n,r,p,rs,ps,d.Delta_model_minus_tumour.mean()])
pd.DataFrame(stats,columns=['Mode','N_pairs','Pearson_r','Pearson_P','Spearman_rho','Spearman_P','Mean_delta']).to_csv(R/'HCMI_EMTI_preservation_stats.csv',index=False)

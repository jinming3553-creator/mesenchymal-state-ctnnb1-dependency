#!/usr/bin/env python
from pathlib import Path
import sys, pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.io_utils import load_broad_expression, load_broad_gene_effect, load_sanger_expression, load_sanger_minlib, load_genesets
from src.scoring import score_programs

RAW=ROOT/'data/raw'; OUT=ROOT/'data/processed'; OUT.mkdir(parents=True,exist_ok=True)
BROAD_LINES=['Bowel','Esophagus/Stomach','Pancreas','Ampulla of Vater','Biliary Tract']
SANGER_TISSUES=['COLO','OESO','PANC']

# Broad
expr=load_broad_expression(RAW/'next_gen_expression.csv')
meta=pd.read_csv(RAW/'model_metadata.csv').set_index('ModelID')
ge, screen_meta=load_broad_gene_effect(RAW/'screen_gene_effect.zip', RAW/'screen_metadata.csv')
gsets=load_genesets(RAW/'geneset_table_updated.csv')
ids=expr.index.intersection(ge.index).intersection(meta.index)
ids=ids[meta.loc[ids,'OncotreeLineage'].isin(BROAD_LINES)]
expr_b=expr.loc[ids]; ge_b=ge.loc[ids]
meta_b=meta.loc[ids].copy(); meta_b['lineage']=meta_b['OncotreeLineage']
# Freeze programme standardization within the GI analysis cohort.
# This is the implementation used by the systematic discovery scan.
scores_b,cov_b=score_programs(expr_b,gsets,min_genes=5)
expr_b.to_pickle(OUT/'broad_expression.pkl'); ge_b.to_pickle(OUT/'broad_gene_effect.pkl')
meta_b.to_csv(OUT/'broad_metadata.csv'); scores_b.to_csv(OUT/'broad_program_scores.csv'); cov_b.to_csv(OUT/'broad_program_coverage.csv',index=False)

# Sanger
sexpr=load_sanger_expression(RAW/'data_RNAseq_ALL-organoids.csv')
sge=load_sanger_minlib(RAW/'41586_2026_10830_MOESM3_ESM.zip')
scov=pd.read_csv(RAW/'data_covariates_ALL-organoids.csv').set_index('sample_ID')
sids=sexpr.index.intersection(sge.index).intersection(scov.index)
sids=sids[scov.loc[sids,'tissue'].isin(SANGER_TISSUES)]
sexpr=sexpr.loc[sids]; sge=sge.loc[sids]; scov=scov.loc[sids].copy(); scov['lineage']=scov['tissue']
# Use the same within-GI standardization rule as Broad.
sscores,scovg=score_programs(sexpr,gsets,min_genes=5)
sexpr.to_pickle(OUT/'sanger_expression.pkl'); sge.to_pickle(OUT/'sanger_gene_effect.pkl')
scov.to_csv(OUT/'sanger_metadata.csv'); sscores.to_csv(OUT/'sanger_program_scores.csv'); scovg.to_csv(OUT/'sanger_program_coverage.csv',index=False)

print('Broad matched GI models:',len(ids))
print('Sanger matched GI organoids:',len(sids))

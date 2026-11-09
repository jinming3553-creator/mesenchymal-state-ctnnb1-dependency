import csv,re,numpy as np
from scipy import stats
import statsmodels.api as sm
from pathlib import Path
P=Path('/mnt/data')
sets={}
with open(P/'geneset_table_updated(2).csv') as f:
 for x in csv.DictReader(f):
  if x['Geneset'] in ['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT']:
   sets.setdefault(x['Geneset'],set()).add(x['GeneSymbol'])
print('Gene sets', {k:len(v) for k,v in sets.items()},flush=True)
outputs=[]
for cohort,fn,expr,idfield in [('Broad','WNT_model_level_broad_reconstructed.csv','next_gen_expression(2).csv',''),('Sanger','WNT_model_level_sanger_reconstructed.csv','data_RNAseq_ALL-organoids(1).csv','sample_ID')]:
 with open(P/fn) as f: rows=list(csv.DictReader(f))
 ids=[r[idfield] for r in rows]; idx={k:i for i,k in enumerate(ids)}; lineage=[r['lineage'] for r in rows]
 score={k:np.zeros((len(ids),),float) for k in sets}; count={k:0 for k in sets}; missing={k:0 for k in sets}
 with open(P/expr) as f:
  rd=csv.reader(f); head=next(rd)
  symbol=[re.sub(r' \(\d+\)$','',x) if cohort=='Broad' else (x[:-4] if x.endswith('_exp') else x) for x in head]
  columns={k:[j for j,s in enumerate(symbol) if s in genes] for k,genes in sets.items()}
  print(cohort,'available genes',{k:len(v) for k,v in columns.items()},flush=True)
  union=sorted(set(sum(columns.values(),[]))); union_map={j:i for i,j in enumerate(union)}
  subcols={k:[union_map[j] for j in v] for k,v in columns.items()}
  matrix=np.full((len(ids),len(union)),np.nan)
  for r in rd:
   i=idx.get(r[0]);
   if i is None:continue
   matrix[i]=[float(r[j]) if r[j] else np.nan for j in union]
 print('matched',np.isfinite(matrix).all(axis=1).sum(),'of',len(ids),flush=True)
 # Within-gene z-score across matched models, then average per set; avoids high-expression genes dominating mean
 mu=np.nanmean(matrix,axis=0); sd=np.nanstd(matrix,axis=0,ddof=0); z=(matrix-mu)/np.where(sd>1e-10,sd,np.nan)
 for k,ix in subcols.items():score[k]=np.nanmean(z[:,ix],axis=1)
 y=np.array([float(r['CTNNB1']) for r in rows]); emt=np.array([float(r['EMTI']) for r in rows]);
 L=np.column_stack([np.ones(len(rows))]+[(np.array(lineage)==t).astype(float) for t in sorted(set(lineage))[1:]])
 e=score['HALLMARK_E2F_TARGETS'];g=score['HALLMARK_G2M_CHECKPOINT'];
 def res(v,x):return sm.OLS(v,x).fit().resid
 def cor(x,y,ctrl):
  a=res(x,ctrl);b=res(y,ctrl);r=np.corrcoef(a,b)[0,1];df=len(x)-np.linalg.matrix_rank(ctrl)-1;p=2*stats.t.sf(abs(r*np.sqrt(df/(1-r*r))),df);return r,p,df
 for name,extra in [('lineage only',[]),('lineage + E2F',[e]),('lineage + G2M',[g]),('lineage + E2F + G2M',[e,g])]:
  C=np.column_stack([L]+extra) if extra else L
  r,p,df=cor(emt,y,C)
  # HC3 coefficient for z-scaled EMT and scores with lineage
  xx=np.column_stack([L,(emt-emt.mean())/emt.std()]+[(v-v.mean())/v.std() for v in extra]);fit=sm.OLS(y,xx).fit(cov_type='HC3')
  beta=fit.params[L.shape[1]];p_hc3=fit.pvalues[L.shape[1]];ci=fit.conf_int()[L.shape[1]]
  outputs.append([cohort,len(rows),name,r,p,df,beta,p_hc3,*ci]);print(cohort,name,'r',round(r,4),'p',p,'beta',beta,'HC3 p',p_hc3,flush=True)
 for k,s in score.items():print(cohort,k,'corr EMT',stats.pearsonr(s,emt).statistic,'corr dependency',stats.pearsonr(s,y).statistic)
 print(cohort,'E2F-G2M corr',stats.pearsonr(e,g).statistic,flush=True)
 with open(P/f'Comment8_{cohort}_model_level.csv','w') as f:
  w=csv.writer(f);w.writerow(['Model_ID','Lineage','EMTI','CTNNB1_dependency','E2F_score','G2M_score']);w.writerows(zip(ids,lineage,emt,y,e,g))
with open(P/'Comment8_cellcycle_results.csv','w') as f:
 w=csv.writer(f);w.writerow(['Cohort','N','Adjustment','partial_r','partial_p','df','EMTI_standardized_beta_HC3','EMTI_HC3_p','EMTI_HC3_CI_low','EMTI_HC3_CI_high']);w.writerows(outputs)

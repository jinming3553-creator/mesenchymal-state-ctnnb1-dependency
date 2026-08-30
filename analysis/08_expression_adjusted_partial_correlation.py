#!/usr/bin/env python
"""Reproduce the expression-adjusted EMT-I–CTNNB1 partial correlations used in Fig. 3c.

For each cohort, EMT-I score and CTNNB1 dependency are residualized separately
against standardized CTNNB1 expression and reference-coded tumour-lineage
indicator variables. The reported partial r is the Pearson correlation between
the two residual vectors. Statistical significance uses the standard partial-
correlation t statistic with n - k - 2 degrees of freedom, where k is the
number of covariates in the residualization model. Fisher-transformed 95%
confidence intervals use standard error 1/sqrt(n - k - 3).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t as t_dist, norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.stats_utils import pearson_pair

P = ROOT / 'data' / 'processed'
R = ROOT / 'results' / 'source_data'
R.mkdir(parents=True, exist_ok=True)

COHORTS = {
    'Broad NextGen GI': 'broad_multivariable_design.csv',
    'Sanger GI organoids': 'sanger_multivariable_design.csv',
}

def partial_corr_inference(r, n, k, alpha=0.05):
    df = n - k - 2
    if df <= 0:
        raise ValueError('Insufficient degrees of freedom for partial-correlation test')
    t_stat = r * np.sqrt(df / (1.0 - r**2))
    p = 2.0 * t_dist.sf(abs(t_stat), df)
    if n - k - 3 <= 0:
        lo = hi = np.nan
    else:
        z = np.arctanh(r)
        se = 1.0 / np.sqrt(n - k - 3)
        crit = norm.ppf(1.0 - alpha / 2.0)
        lo, hi = np.tanh([z - crit * se, z + crit * se])
    return t_stat, df, p, lo, hi

rows = []
for cohort, filename in COHORTS.items():
    d = pd.read_csv(P / filename)
    required = ['EMTI_z', 'CTNNB1_dependency_z', 'CTNNB1_expression_z', 'cov_lineage']
    d = d[required].dropna().copy()

    cov = pd.get_dummies(
        d[['CTNNB1_expression_z', 'cov_lineage']],
        columns=['cov_lineage'], drop_first=True, dtype=float
    )
    X = sm.add_constant(cov, has_constant='add')
    emti_resid = sm.OLS(d['EMTI_z'], X).fit().resid
    dep_resid = sm.OLS(d['CTNNB1_dependency_z'], X).fit().resid

    r, _, n = pearson_pair(emti_resid, dep_resid)
    k = cov.shape[1]  # CTNNB1 expression + reference-coded lineage indicators
    t_stat, df, p, lo, hi = partial_corr_inference(r, n, k)
    rows.append({
        'Cohort': cohort,
        'N': n,
        'K_covariates': k,
        'df': df,
        'partial_r': r,
        'P': p,
        'CI_low': lo,
        'CI_high': hi,
    })

out = pd.DataFrame(rows)
out.to_csv(R / 'Figure3c_expression_independent_source_data.csv', index=False)
print(out.to_string(index=False))

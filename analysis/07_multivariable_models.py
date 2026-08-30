#!/usr/bin/env python
"""Reproduce the fully adjusted EMT-I–CTNNB1 regression models used in Fig. 3.

The cohort-specific processed design tables expose the exact covariate coding used
in the manuscript. Continuous predictors enter in standardized form; binary event
indicators enter as 0/1 variables; categorical variables are reference-coded by
``hc3_ols``.
"""
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.stats_utils import hc3_ols

P = ROOT / 'data' / 'processed'
R = ROOT / 'results' / 'source_data'
R.mkdir(parents=True, exist_ok=True)

MODELS = {
    'broad': {
        'file': 'broad_multivariable_design.csv',
        'predictors': [
            'EMTI_z', 'CTNNB1_expression_z',
            'cov_mut_APC', 'cov_mut_CTNNB1', 'cov_mut_KRAS', 'cov_mut_BRAF', 'cov_mut_RNF43',
            'cov_cn_CTNNB1_z', 'cov_cn_APC_z',
            'cov_lineage', 'cov_screen_type'
        ],
        'categorical': ['cov_lineage', 'cov_screen_type'],
    },
    'sanger': {
        'file': 'sanger_multivariable_design.csv',
        'predictors': [
            'EMTI_z', 'CTNNB1_expression_z',
            'cov_mut_APC', 'cov_mut_CTNNB1', 'cov_mut_KRAS', 'cov_mut_BRAF', 'cov_mut_RNF43',
            'cov_CTNNB1_amp', 'cov_APC_del', 'cov_KRAS_amp', 'cov_MSI', 'cov_ploidy_z',
            'cov_lineage'
        ],
        'categorical': ['cov_lineage'],
    },
}

for cohort, spec in MODELS.items():
    path = P / spec['file']
    if not path.exists():
        print(f"SKIP {spec['file']}: frozen design table not found")
        continue
    d = pd.read_csv(path)
    fit = hc3_ols(d, 'CTNNB1_dependency_z', spec['predictors'], categorical=spec['categorical'])
    row = 'EMTI_z'
    ci = fit.conf_int().loc[row]
    out = pd.DataFrame([{
        'cohort': cohort,
        'N': int(fit.nobs),
        'beta': fit.params[row],
        'HC3_P': fit.pvalues[row],
        'CI_low': ci.iloc[0],
        'CI_high': ci.iloc[1],
    }])
    out.to_csv(R / spec['file'].replace('_design.csv', '_HC3_result.csv'), index=False)
    print(out.to_string(index=False))

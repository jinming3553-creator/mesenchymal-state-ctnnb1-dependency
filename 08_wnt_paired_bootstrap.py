#!/usr/bin/env python3
"""Reviewer 1 Comment 5: paired, lineage-adjusted WNT correlation contrasts.

Run from the repository root after analysis/01_prepare_cohorts.py:
    python analysis/08_wnt_paired_bootstrap.py

The original Figure 3d correlations are preserved; this script adds direct
comparisons against eight prespecified WNT-pathway dependencies.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
GENES = ('CTNNB1', 'APC', 'TCF7L2', 'LGR4', 'MESD', 'PORCN', 'WLS', 'FZD5', 'LRP6')
EXPECTED = {'Broad': (78, 0.4109262561313052), 'Sanger': (144, 0.4031670129478884)}


def load_cohort(prefix, processed):
    scores = pd.read_csv(processed / f'{prefix}_program_scores.csv', index_col=0)
    ge = pd.read_pickle(processed / f'{prefix}_gene_effect.pkl')
    meta = pd.read_csv(processed / f'{prefix}_metadata.csv', index_col=0)
    ids = scores.index.intersection(ge.index).intersection(meta.index, sort=False)
    missing = set(GENES).difference(ge.columns)
    if missing:
        raise ValueError(f'{prefix}: missing dependencies: {sorted(missing)}')
    df = pd.concat([scores.loc[ids, 'EMT-I'].rename('EMTI'),
                    ge.loc[ids, list(GENES)],
                    meta.loc[ids, 'lineage']], axis=1)
    df = df.rename(columns={'lineage': 'lineage'}).dropna()
    return df


def residualized_correlations(df):
    vals = df[['EMTI', *GENES]].to_numpy(dtype=float).copy()
    groups = pd.factorize(df['lineage'])[0]
    for group in np.unique(groups):
        selected = groups == group
        vals[selected] -= vals[selected].mean(axis=0)
    x = vals[:, 0]
    y = vals[:, 1:]
    numerator = np.sum(x[:, None] * y, axis=0)
    denominator = np.sqrt(np.sum(x*x) * np.sum(y*y, axis=0))
    if np.any(denominator <= 0):
        raise ValueError('Zero residual variance in a cohort')
    return numerator / denominator


def paired_bootstrap(df, n_boot, seed):
    rng = np.random.default_rng(seed)
    vals = df[['EMTI', *GENES]].to_numpy(dtype=float)
    groups = pd.factorize(df['lineage'])[0]
    n = len(df)
    out = np.full((n_boot, len(GENES)-1), np.nan)
    for i in range(n_boot):
        indices = rng.integers(0, n, size=n)
        sample = vals[indices].copy()
        group_sample = groups[indices]
        for group in np.unique(group_sample):
            mask = group_sample == group
            sample[mask] -= sample[mask].mean(axis=0)
        x, y = sample[:, 0], sample[:, 1:]
        numerator = np.sum(x[:, None]*y, axis=0)
        denominator = np.sqrt(np.sum(x*x)*np.sum(y*y, axis=0))
        with np.errstate(divide='ignore', invalid='ignore'):
            r = numerator / denominator
        out[i] = r[0] - r[1:]
    return out


def summary(bootstrap_values):
    a = bootstrap_values[np.isfinite(bootstrap_values)]
    if not len(a):
        return np.nan, np.nan, np.nan, 0
    # Tail-proportion bootstrap significance, matching the revision source data.
    # This is not a null-centered/permutation test.
    p = min(1., 2*min((np.count_nonzero(a <= 0)+1)/(len(a)+1),
                      (np.count_nonzero(a >= 0)+1)/(len(a)+1)))
    return np.quantile(a, .025), np.quantile(a, .975), p, len(a)


def run(processed, output, n_boot=10000, seed=20261109, verify=True):
    output.mkdir(parents=True, exist_ok=True)
    cohorts, draws, rows = {}, {}, []
    for label, prefix, cohort_seed in [('Broad','broad',seed), ('Sanger','sanger',seed+1)]:
        df = load_cohort(prefix, processed)
        r = residualized_correlations(df)
        if verify:
            expected_n, expected_r = EXPECTED[label]
            if len(df) != expected_n or not np.isclose(r[0], expected_r, atol=1e-8, rtol=0):
                raise RuntimeError(f'{label}: failed original Figure 3d validation: n={len(df)}, CTNNB1 r={r[0]:.12f}')
        cohorts[label] = (df, r)
        draws[label] = paired_bootstrap(df, n_boot, cohort_seed)
        for j, gene in enumerate(GENES[1:]):
            lo, hi, p, valid = summary(draws[label][:, j])
            rows.append({'Cohort':label, 'Gene':gene, 'N':len(df),
                         'CTNNB1_r':r[0], 'Comparator_r':r[j+1],
                         'Delta_r':r[0]-r[j+1], 'CI_low':lo, 'CI_high':hi,
                         'P_bootstrap':p, 'Bootstrap_valid':valid})
    by = pd.DataFrame(rows)
    for label in cohorts:
        mask = by['Cohort'].eq(label)
        by.loc[mask, 'BH_FDR_within_cohort'] = multipletests(by.loc[mask,'P_bootstrap'], method='fdr_bh')[1]
    by.to_csv(output/'WNT_CTNNB1_paired_bootstrap_by_cohort.csv', index=False)
    weights = np.array([len(cohorts[label][0])-3 for label in ('Broad','Sanger')], dtype=float)
    weights /= weights.sum()
    pooled = []
    for j, gene in enumerate(GENES[1:]):
        a = weights[0]*draws['Broad'][:, j] + weights[1]*draws['Sanger'][:, j]
        lo, hi, p, valid = summary(a)
        delta = sum(weights[k]*(cohorts[label][1][0]-cohorts[label][1][j+1])
                    for k,label in enumerate(('Broad','Sanger')))
        pooled.append({'Gene':gene, 'Weighted_Delta_r':delta, 'CI_low':lo,
                       'CI_high':hi, 'P_bootstrap':p, 'Bootstrap_valid':valid})
    cross = pd.DataFrame(pooled)
    cross['BH_FDR_8_comparisons'] = multipletests(cross['P_bootstrap'],method='fdr_bh')[1]
    cross.to_csv(output/'WNT_CTNNB1_paired_bootstrap_cross_cohort.csv', index=False)
    print('Reproduced correlations:', {label:round(cohorts[label][1][0],9) for label in cohorts})
    print('Saved results to', output)
    return by, cross


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--processed', type=Path, default=ROOT/'data/processed')
    parser.add_argument('--output', type=Path, default=ROOT/'results/source_data')
    parser.add_argument('--bootstrap', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20261109)
    parser.add_argument('--skip-reference-check', action='store_true')
    args = parser.parse_args()
    run(args.processed,args.output,args.bootstrap,args.seed,not args.skip_reference_check)

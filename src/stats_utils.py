from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm


def residualize(values: pd.Series, groups: pd.Series) -> pd.Series:
    """Subtract group-specific means; preserves index and missing values."""
    frame = pd.DataFrame({"y": values, "g": groups})
    return frame["y"] - frame.groupby("g", observed=True)["y"].transform("mean")


def pearson_pair(x, y):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(d) < 3 or d.x.nunique() < 2 or d.y.nunique() < 2:
        return np.nan, np.nan, len(d)
    r, p = stats.pearsonr(d.x, d.y)
    return float(r), float(p), len(d)


def spearman_pair(x, y):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(d) < 3:
        return np.nan, np.nan, len(d)
    r, p = stats.spearmanr(d.x, d.y)
    return float(r), float(p), len(d)


def fisher_ci(r: float, n: int, alpha: float = 0.05):
    if not np.isfinite(r) or n <= 3 or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    q = stats.norm.ppf(1 - alpha / 2)
    return tuple(np.tanh([z - q * se, z + q * se]))


def fixed_effect_meta(correlations, ns):
    """Fixed-effect meta-analysis of Pearson r using Fisher r-to-z."""
    arr = [(r, n) for r, n in zip(correlations, ns) if np.isfinite(r) and n > 3 and abs(r) < 1]
    if not arr:
        return {"r": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    z = np.array([np.arctanh(r) for r, _ in arr])
    w = np.array([n - 3 for _, n in arr], dtype=float)
    zm = np.sum(w * z) / np.sum(w)
    se = 1 / np.sqrt(np.sum(w))
    zstat = zm / se
    p = 2 * stats.norm.sf(abs(zstat))
    q = stats.norm.ppf(0.975)
    lo, hi = np.tanh([zm - q * se, zm + q * se])
    return {"r": float(np.tanh(zm)), "p": float(p), "ci_low": float(lo), "ci_high": float(hi)}


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def hc3_ols(df: pd.DataFrame, outcome: str, predictors: list[str], categorical: list[str] | None = None):
    """OLS with HC3 robust SE and reference-coded categorical covariates.

    HC3 is computed explicitly to avoid a statsmodels-version-dependent 0/0
    artefact when a singleton indicator has leverage numerically equal to one
    and residual numerically equal to zero.
    """
    from types import SimpleNamespace
    categorical = categorical or []
    d = df[[outcome] + predictors].dropna().copy()
    X = d[predictors].copy()
    cats = [c for c in categorical if c in X.columns]
    if cats:
        X = pd.get_dummies(X, columns=cats, drop_first=True, dtype=float)
    X = X.apply(pd.to_numeric, errors="coerce")
    ok = X.notna().all(axis=1)
    X = X.loc[ok]
    y = pd.to_numeric(d.loc[ok, outcome], errors="coerce")
    ok2 = y.notna()
    X, y = X.loc[ok2], y.loc[ok2]
    X = sm.add_constant(X, has_constant="add")
    fit = sm.OLS(y, X).fit()
    h = fit.get_influence().hat_matrix_diag
    denom = np.maximum(1.0 - h, np.sqrt(np.finfo(float).eps))
    scaled_sq_resid = (np.asarray(fit.resid) / denom) ** 2
    x = np.asarray(X, dtype=float)
    bread = np.linalg.pinv(x.T @ x)
    meat = x.T @ (x * scaled_sq_resid[:, None])
    cov = bread @ meat @ bread
    se = np.sqrt(np.diag(cov))
    params = pd.Series(np.asarray(fit.params), index=X.columns)
    z = params.to_numpy() / se
    pvalues = pd.Series(2 * stats.norm.sf(np.abs(z)), index=X.columns)
    q = stats.norm.ppf(0.975)
    ci = pd.DataFrame({0: params.to_numpy() - q * se, 1: params.to_numpy() + q * se}, index=X.columns)
    return SimpleNamespace(params=params, pvalues=pvalues, bse=pd.Series(se,index=X.columns), nobs=fit.nobs,
                           conf_int=lambda: ci, model=fit.model, cov_params=lambda: pd.DataFrame(cov,index=X.columns,columns=X.columns))


def zscore_series(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd

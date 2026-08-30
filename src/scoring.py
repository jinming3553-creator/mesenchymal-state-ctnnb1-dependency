from __future__ import annotations
import numpy as np
import pandas as pd


def zscore_expression(expr: pd.DataFrame) -> pd.DataFrame:
    """Gene-wise z-score across models using population SD (ddof=0)."""
    mu = expr.mean(axis=0)
    sd = expr.std(axis=0, ddof=0).replace(0, np.nan)
    return (expr - mu) / sd


def score_programs(expr: pd.DataFrame, genesets: dict[str, list[str]], min_genes=1):
    z = zscore_expression(expr)
    scores, coverage = {}, []
    for name, genes in genesets.items():
        avail = [g for g in genes if g in z.columns]
        if len(avail) < min_genes:
            continue
        scores[name] = z[avail].mean(axis=1)
        coverage.append({"Program": name, "n_signature": len(set(genes)), "n_available": len(avail)})
    return pd.DataFrame(scores, index=expr.index), pd.DataFrame(coverage)

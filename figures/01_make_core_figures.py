#!/usr/bin/env python
"""Minimal reproducible plotting layer.

Figure placement is deliberately decoupled from analysis. Moving a panel between a main figure and
Extended Data changes only output names/layout, not statistical computation.
"""
from pathlib import Path
import pandas as pd, matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'results/source_data'; OUT=ROOT/'results/figures'; OUT.mkdir(parents=True,exist_ok=True)

# Broad and Sanger primary scatter plots
for prefix,label in [('broad','Broad NextGen GI'),('sanger','Sanger GI organoids')]:
    f=SRC/f'{prefix}_EMTI_CTNNB1_model_level.csv'
    if not f.exists(): continue
    d=pd.read_csv(f)
    fig,ax=plt.subplots(figsize=(4,4))
    ax.scatter(d['EMT-I adjusted'],d['CTNNB1 adjusted'],s=24,alpha=.8)
    ax.set_xlabel('Lineage-adjusted EMT-I score')
    ax.set_ylabel('Lineage-adjusted CTNNB1 dependency')
    ax.set_title(label)
    fig.tight_layout(); fig.savefig(OUT/f'{prefix}_EMTI_CTNNB1.pdf'); fig.savefig(OUT/f'{prefix}_EMTI_CTNNB1.png',dpi=300); plt.close(fig)

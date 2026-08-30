from __future__ import annotations
from pathlib import Path
import io, re, zipfile
import numpy as np
import pandas as pd

GENE_ID_RE = re.compile(r"\s+\(\d+\)$")


def strip_gene_id(name: str) -> str:
    return GENE_ID_RE.sub("", str(name))


def load_broad_expression(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [strip_gene_id(c) for c in df.columns]
    df.index.name = "ModelID"
    return df


def load_broad_gene_effect(zip_path, screen_metadata_path, aggregate="mean"):
    """Read QC-passing Broad NextGen CRISPR screens and aggregate screens to model level."""
    sm = pd.read_csv(screen_metadata_path)
    keep = sm.loc[(sm.PassesQC == True) & (sm.IsNextGen == True), ["ScreenID","ModelID","ScreenType","OncotreeLineage"]]
    with zipfile.ZipFile(zip_path) as z:
        with z.open("screen_gene_effect.csv") as f:
            ge = pd.read_csv(f, index_col=0)
    ge.columns = [strip_gene_id(c) for c in ge.columns]
    ge.index.name = "ScreenID"
    x = keep.merge(ge, left_on="ScreenID", right_index=True, how="inner")
    gene_cols = [c for c in ge.columns if c in x.columns]
    agg = x.groupby("ModelID", observed=True)[gene_cols].agg(aggregate)
    meta = x.groupby("ModelID", observed=True).agg(
        ScreenType=("ScreenType", lambda s: ";".join(sorted(set(map(str,s))))),
        lineage=("OncotreeLineage","first"),
    )
    return agg, meta


def load_genesets(path):
    x = pd.read_csv(path)
    return x.groupby("Geneset", observed=True)["GeneSymbol"].apply(lambda s: sorted(set(s.dropna().astype(str)))).to_dict()


def load_sanger_expression(path):
    x = pd.read_csv(path, index_col="sample_ID")
    x.columns = [c[:-4] if c.endswith("_exp") else c for c in x.columns]
    return x


def load_sanger_minlib(zip_path):
    """Read Supplementary Table 6 directly from the Sanger supplementary ZIP."""
    with zipfile.ZipFile(zip_path) as z:
        raw = z.read("supplementary_table_6_revision.csv")
    x = pd.read_csv(io.BytesIO(raw), sep=r"\s+")
    x = x.loc[x["library"].astype(str).str.lower().eq("minlib")]
    return x.pivot_table(index="sample_ID", columns="gene", values="LFC", aggfunc="mean")


def load_sanger_supp_sheet(zip_path, sheet):
    with zipfile.ZipFile(zip_path) as z:
        raw = z.read("supplementary_table_3_revision.xlsx")
    return pd.read_excel(io.BytesIO(raw), sheet_name=sheet)

from pathlib import Path
import pandas as pd
from scipy.stats import pearsonr

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'reference_outputs'

def test_broad_figure2b_source_values():
    d=pd.read_csv(REF/'Figure2b_Broad_EMTI_CTNNB1_source_data.csv')
    r,p=pearsonr(d['EMT-I adjusted'],d['CTNNB1 adjusted'])
    assert abs(r-0.4109262561313052)<1e-12

def test_sanger_figure2c_source_values():
    d=pd.read_csv(REF/'Figure2c_Sanger_EMTI_CTNNB1_source_data.csv')
    r,p=pearsonr(d['EMT-I adjusted'],d['CTNNB1 adjusted'])
    assert abs(r-0.4031670129478884)<1e-12

def test_meta_summary():
    d=pd.read_csv(REF/'Figure2d_EMTI_CTNNB1_meta_source_data.csv')
    meta=d[d['Cohort'].eq('Fixed-effect meta')].iloc[0]
    assert abs(meta['r']-0.40586784858857816)<1e-12
    assert abs(meta['P']-2.462797448441407e-10)<1e-20

def test_fully_adjusted_effects():
    d=pd.read_csv(REF/'Figure3b_fully_adjusted_source_data.csv')
    b=d[d['Cohort'].str.contains('Broad')].iloc[0]
    s=d[d['Cohort'].str.contains('Sanger')].iloc[0]
    assert abs(b['EMT-I standardized beta']-0.4756974755596273)<1e-12
    assert abs(s['EMT-I standardized beta']-0.3399402368871638)<1e-12


def test_expression_adjusted_partial_correlations():
    d=pd.read_csv(REF/'Figure3c_expression_independent_source_data.csv')
    b=d[d['Cohort'].str.contains('Broad')].iloc[0]
    s=d[d['Cohort'].str.contains('Sanger')].iloc[0]
    assert abs(b['partial_r']-0.4640334312562166)<1e-12
    assert abs(b['P']-3.549408077658038e-05)<1e-16
    assert int(b['df']) == 71
    assert abs(b['CI_low']-0.26192945696151615)<1e-12
    assert abs(b['CI_high']-0.6271471621580395)<1e-12
    assert abs(s['partial_r']-0.38961614549250956)<1e-12
    assert abs(s['P']-1.7940357125388227e-06)<1e-17
    assert int(s['df']) == 139
    assert abs(s['CI_low']-0.23974560668838368)<1e-12
    assert abs(s['CI_high']-0.521349088640739)<1e-12

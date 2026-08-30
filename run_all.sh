#!/usr/bin/env bash
set -euo pipefail
python analysis/01_prepare_cohorts.py
python analysis/02_primary_emti_ctnnb1.py
python analysis/04_robustness_and_specificity.py
python analysis/05_lineage_and_cross_context.py
python analysis/06_hcmi_emti_preservation.py
python analysis/07_multivariable_models.py
python analysis/08_expression_adjusted_partial_correlation.py
# Systematic scan is intentionally separate because it is computationally expensive and requires
# the frozen exact 20-state list in config/systematic_states.txt.
# python analysis/03_systematic_state_dependency_scan.py
python figures/01_make_core_figures.py

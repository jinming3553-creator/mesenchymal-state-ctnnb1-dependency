# Tumour-state plasticity and CTNNB1 dependency in gastrointestinal cancer models

Reproducible Python code for the manuscript analysis linking an EMT/ECM-associated transcriptional state to reduced CTNNB1 dependency across independent gastrointestinal patient-derived cancer model collections.

## Analysis overview

The repository separates **statistical analysis** from **figure placement**. Therefore, moving a panel from a main figure to Extended Data does not change the underlying analysis or the code required to reproduce it.

1. Prepare matched Broad NextGen and Sanger GI cohorts and calculate transcriptional programme scores.
2. Reproduce the primary lineage-adjusted EMT-I–CTNNB1 associations and fixed-effect meta-analysis.
3. Run the systematic cross-cohort state–dependency screen (20 states × shared CRISPR dependencies).
4. Test robustness to stromal/ECM-heavy EMT-I genes and leave-one-gene-out scoring; compare representative WNT dependencies.
5. Quantify lineage heterogeneity and secondary cross-context programme convergence.
6. Analyse matched HCMI tumour/model snRNA-seq for relative EMT-I-state preservation and ex vivo attenuation.
7. Fit cohort-specific HC3 multivariable models using frozen design tables.
8. Reproduce the CTNNB1-expression- and lineage-adjusted residual correlations shown in Fig. 3c.
8. Generate figures from source-data tables.

## Installation

```bash
conda env create -f environment.yml
conda activate ctnnb1-state-dependency
```

or

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
```

## Data

Raw public datasets are not redistributed because of repository size and upstream data-use considerations. Download the files listed in `data/raw/README.md` and place them under `data/raw/`.

Broad NextGen data: https://doi.org/10.6084/m9.figshare.29472362

Sanger tumour-derived organoid data: https://doi.org/10.6084/m9.figshare.28339340

HCMI snRNA-seq: GEO GSE297148

## Reproduce core analyses

```bash
bash run_all.sh
```

The systematic 326,360-test scan is kept separate because it is comparatively expensive and should only be run after freezing the exact 20-programme list:

```bash
python analysis/03_systematic_state_dependency_scan.py
```

## Important statistical conventions

- CRISPR scores are oriented such that more negative values indicate stronger dependency; therefore, a positive EMT-I–CTNNB1 correlation indicates reduced dependency at higher EMT-I activity.
- Primary cross-cohort correlations are adjusted for lineage by within-lineage mean residualization.
- Broad and Sanger Pearson correlations are combined by fixed-effect meta-analysis after Fisher transformation.
- Multiple-testing correction for the systematic analysis is performed across the complete cross-cohort test universe.
- HCMI epithelial-enriched analysis is a permissive epithelial-marker sensitivity analysis and is not presented as exact inferCNV-based malignant-cell classification.
- Cross-lineage programme analysis is secondary and selected on the basis of lineage-stratified results; its FDR values should not be described as independent discovery evidence.

## Repository status

This repository is structured for manuscript reproducibility. Before public release, two items should be frozen against the final manuscript version:

1. the exact 20 transcriptional programmes in `config/systematic_states.txt`;
2. the exact cohort-specific genomic-event coding in the multivariable design tables described in `docs/MULTIVARIABLE_DESIGN.md`.

These two explicit freeze points prevent exploratory choices from being silently changed after manuscript finalization.

## Citation

Please cite the associated manuscript and the three public cancer-model resources used in the analysis.

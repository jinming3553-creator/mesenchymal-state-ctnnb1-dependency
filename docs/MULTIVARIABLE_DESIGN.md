# Multivariable design tables

The final manuscript analysis used cohort-specific models because mutation and copy-number annotations differ between Broad and Sanger.

Broad fully adjusted model:
- EMT-I score
- CTNNB1 expression
- APC, CTNNB1, KRAS, BRAF and RNF43 hotspot/damaging mutation indicators
- CTNNB1 and APC copy number
- tumour lineage
- CRISPR screen type

Sanger fully adjusted model:
- EMT-I score
- CTNNB1 expression
- APC, CTNNB1, KRAS, BRAF and RNF43 driver-event indicators
- CTNNB1/APC/KRAS copy-number events
- MSI
- ploidy
- tumour lineage

Before public release, generate `data/processed/broad_multivariable_design.csv` and `data/processed/sanger_multivariable_design.csv` with one row per analysed model and covariates prefixed `cov_`. This deliberately exposes the exact event coding rather than hiding it inside a long script.

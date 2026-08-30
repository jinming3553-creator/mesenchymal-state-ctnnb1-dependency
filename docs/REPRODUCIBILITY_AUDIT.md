# Reproducibility audit

## Resolved: EMT-I–CTNNB1 correlation discrepancy

The small Broad discrepancy was traced to the **expression standardization universe**, not to CRISPR aggregation, lineage residualization, rounding, or sample inclusion.

- When EMT-I genes were z-scored across the complete Broad NextGen expression collection before selecting GI models, the lineage-adjusted Broad association was r = 0.409482 (P = 1.971e-4).
- When EMT-I genes were z-scored within the 78-model GI analysis cohort, as in the systematic discovery scan, the association was r = 0.410926 (P = 1.861e-4).

The same issue was detectable in Sanger:

- complete Sanger expression collection: r = 0.404033 (P = 5.115e-7);
- 144-model GI analysis cohort: r = 0.403167 (P = 5.439e-7).

The repository now freezes a single rule for the primary manuscript analysis: **gene-wise z-scoring is performed within each cohort's GI analysis set before programme averaging**. This matches the Broad systematic discovery scan and avoids mixing standardization universes between cohorts/panels.

Under the unified rule, EMT-I–CTNNB1 is:

- Broad: r = 0.410926, P = 1.861e-4;
- Sanger: r = 0.403167, P = 5.439e-7;
- fixed-effect meta-analysis: r = 0.405868, P = 2.463e-10.

Recalculation of the cross-cohort multiple-testing table under the unified Sanger GI-standardization rule retains EMT-I–CTNNB1 as the only association with global FDR < 0.05 (approximately 8.0e-5).


## Frozen manuscript values after full Fig. 2–4 audit
- Broad EMT-I–CTNNB1: r = 0.410926, P = 1.861e-4.
- Sanger EMT-I–CTNNB1: r = 0.403167, P = 5.439e-7.
- Fixed-effect meta-analysis: r = 0.405868, 95% CI 0.288836–0.510949, P = 2.463e-10; global BH FDR = 8.038e-5.
- Fully adjusted Broad EMT-I coefficient: beta = 0.475697, HC3 P = 0.002186.
- Fully adjusted Sanger EMT-I coefficient: beta = 0.339940, HC3 P = 0.001361.
- Figure 5 and Extended Data stroma-light/leave-one-gene-out outputs were already generated using the frozen GI-cohort scoring rule and therefore did not require numerical revision.

## Fig. 3 multivariable-model implementation check
The public Fig. 3 script uses explicit cohort-specific predictor lists rather than selecting every `cov_` column. This prevents simultaneous inclusion of raw and standardized copies of continuous covariates and ensures that standardized CTNNB1 expression is included. Binary genomic event variables are coded 0/1; lineage and Broad screen type are reference-coded categorical variables. The frozen design tables contain no missing values in the variables used by the final models (Broad n=78; Sanger n=144). The HC3 helper includes a numerical guard for a singleton Sanger amplification indicator whose leverage is floating-point-equivalent to one; this reproduces the finite limiting HC3 covariance and avoids version-dependent 0/0 behavior.

## Systematic scan sample sizes

The systematic scan computes the effective sample size separately for every state–dependency association after pairwise removal of missing values. The submitted 326,360-row source table was audited against the cohort-specific Pearson P values and fixed-effect meta-analysis: all retained shared dependencies have complete measurements in the frozen GI analysis cohorts, so the effective sample sizes are Broad n=78 and Sanger n=144 for every submitted association. The meta-analysis weights are therefore 75 and 141, respectively. The public scan script nevertheless retains pair-specific N calculation so that the implementation remains correct if applied to inputs containing missing values.

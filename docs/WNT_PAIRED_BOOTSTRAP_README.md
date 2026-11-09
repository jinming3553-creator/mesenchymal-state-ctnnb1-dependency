## Reviewer 1 Comment 5: paired WNT-specificity comparisons

Run after `python analysis/01_prepare_cohorts.py`:

```bash
python analysis/08_wnt_paired_bootstrap.py
```

The script uses the existing `data/processed/{broad,sanger}_program_scores.csv`,
`{broad,sanger}_gene_effect.pkl`, and `{broad,sanger}_metadata.csv` inputs.
It verifies the original CTNNB1 correlations (Broad n=78, r=0.4109262561;
Sanger n=144, r=0.4031670129), and does not alter Figure 3d or
`results/source_data/WNT_specificity.csv`.

Outputs in `results/source_data/`:
- `WNT_CTNNB1_paired_bootstrap_by_cohort.csv`
- `WNT_CTNNB1_paired_bootstrap_cross_cohort.csv`

Within each cohort, all model-level EMT-I and WNT-dependency values are resampled
together (10,000 resamples; seeds 20261109 and 20261110). Lineage means are
re-estimated within each resample, and the difference in lineage-residualized
Pearson correlations is computed as r(CTNNB1) minus r(comparator). Cross-cohort
contrasts are pooled with fixed weights proportional to n−3 (75/216 and
141/216). Confidence intervals are bootstrap percentile intervals. Reported
P values use twice the smaller empirical bootstrap tail proportion relative
to zero with an add-one correction, followed by Benjamini–Hochberg correction
across eight comparisons separately in each cohort and in the pooled analysis.

**Statistical interpretation:** These empirical tail-proportion P values are
bootstrap sign probabilities, not exact null-hypothesis permutation P values.
The paired comparisons support relative specificity among eight selected WNT
genes; they do not establish causality or exclude effects involving other genes.

Raw public expression and gene-effect datasets are not committed to GitHub.

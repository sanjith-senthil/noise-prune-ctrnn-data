# Statistical Reporting Guide

Use this guide and the `paper_relevant_*` CSVs for paper-facing statistical
reporting. The older `revised_scope_task_only_significance_*` files remain
available as a frozen selected-comparison analysis, but they use raw sequence
accuracy and selected one-sided alternatives. Do not mix their p-values into
the reporting-grade analysis below.

## Primary Test

The reporting-grade export uses paired two-sided Wilcoxon signed-rank tests
with `zero_method="wilcox"` and SciPy `method="auto"`. The significance level
is `alpha=0.05`. Holm correction is applied within each named family across all
method contrasts and sparsities in that family. Report the actual
`wilcoxon_p` and `holm_p_within_family` values for every cited comparison.

The paired non-parametric test is appropriate because each method is evaluated
on the same trained networks and held-out task batches, while task retention is
bounded and normality is not assumed. A t-test is not used, so a normality test
is not required. The Wilcoxon signed-rank test does assume that paired
differences are symmetrically distributed. Each row therefore also exports an
exact two-sided binomial sign test and its Holm-adjusted p-value as a robustness
check that does not require paired-difference symmetry.

## Analysis Units

- `task_network_cell`: pruning-seed technical replicates are averaged within
  each trained-network, method, and sparsity cell before inference. Each
  per-sparsity comparison has `n=24` paired trained-network cells.
- `task_mean`: the three trained networks are then averaged within each task.
  Each per-sparsity sensitivity comparison has `n=8` paired task means.

Use the `task_network_cell` result as the trained-network analysis and report
the corresponding `task_mean` result as a clustered sensitivity analysis. The
small-`n` task view is accompanied by exact sign-test p-values and range-based
descriptive statistics. Do not treat the `72` stochastic pruning rows as
independent observations.

## Files

- `paper_relevant_h512_comparison_tests.csv`: `768` explicit H=512
  comparisons. Every row records the dataset, metric, methods, sparsity, test,
  tails, alpha, analysis unit, `n`, raw p-value, Holm-adjusted p-value, exact
  sign-test robustness p-value, and paired-difference summaries.
- `paper_relevant_h512_descriptive_statistics.csv`: `576` method-level
  summaries with `n`, mean, median, sample SD, SEM, minimum, maximum, and
  range.
- `paper_relevant_amplification_descriptive_statistics.csv`: descriptive
  summaries for the candidate-edge `1 / p_ij` amplification distributions.
- `../capped_rescale/*_q10_to_q90_3seed_vs_uncapped_tests.csv`:
  per-sparsity q10--q90 capped-rescale comparisons against q100 uncapped
  rescale. These rows retain the four sparsity levels separately.
- `../capped_rescale/*_q10_to_q90_3seed_vs_q100_plot_level_holm_tests.csv`:
  plot-level capped-rescale comparisons against q100 uncapped rescale. These
  rows match the cap-percentile figure: pruning seeds are averaged within each
  trained-network/sparsity cell, the four sparsities are averaged within each
  trained network, and inference is then run across `n=24` trained networks.
  Holm correction is applied across the full 18-comparison L-NP/S-NP q10--q90
  capped-curve family.

## Comparison Coverage

The comprehensive export includes per-sparsity pairwise comparisons and
comparisons against the unpruned baseline for:

- main H=512 task retention;
- main H=512 recurrent spectral-radius ratio;
- exploratory one-pruning-seed S-NP q50/q75/q90 cap probe retention and
  spectral-radius ratio;
- matched three-pruning-seed L-NP/S-NP q50 cap retention and spectral-radius
  ratio.
- three-pruning-seed L-NP/S-NP q10--q90 cap-percentile curve retention at the
  plot level, compared against q100 uncapped rescale.

Task-retention tests use baseline-normalized `sequence_retention`, matching the
paper-facing metric. Spectral-radius tests use
`post_rec_linear_rho_ratio_to_unpruned`.

## Claims That Remain Descriptive

- “Comparable,” “matches,” and “close to” are equivalence claims. Failure to
  reject a difference does not establish equivalence. Use descriptive wording
  unless an equivalence margin is prespecified and an equivalence analysis is
  added.
- “Near chance” requires a task-specific chance benchmark. The current export
  supports wording such as “low retention,” not an inferential chance-level
  claim.
- Candidate-edge amplification tails are a descriptive mechanistic result.
  Report quantiles, maxima, and ranges rather than a significance test.
- The archived H=256/H=512 scaling bundle is excluded from current paper
  claims and figures. The full suites use unmatched task batteries, and the
  exact overlap contains only two tasks. Do not claim a formal hidden-size
  effect without a matched H=256 rerun and a prespecified scaling analysis.

## Error Bars And Plus-Minus Notation

Main H=512 and matched-q50 figures should use `mean +/- SEM` after averaging
technical pruning-seed replicates within each of the `24` trained networks.
State explicitly that error bars are SEM. For the `n=8` task-mean sensitivity
view, include the exported range where space permits.

## Regeneration And Validation

```bash
python3 scripts/analyze_paper_relevant_statistics.py
python3 scripts/validate_paper_relevant_statistics.py
```

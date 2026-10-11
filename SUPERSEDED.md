# Superseded rows

The sampling rule's probability normalisation was corrected on 2026-10-09.
Before the fix it scaled its unnormalised scores by `K = M / sum(q)` and then
clipped them at 1, discarding the probability mass above 1, so it retained fewer
connections than requested -- at a requested 10% sparsity it removed 33% of the
recurrent weights. It now solves `sum(min(1, K q)) = M`, so the expected retained
count equals the budget.

Six strategies draw a Bernoulli mask and were therefore affected:

  - `simulation_noise_prune_rescale`
  - `simulation_noise_prune_capped_rescale`
  - `simulation_noise_prune_magnitude_rescale`
  - `simulation_noise_prune_capped_magnitude_rescale`
  - `simulation_noise_prune_shuffled_rescale`
  - `simulation_noise_prune_uniform_rescale`

They were re-run over 10-80% at both network sizes. **The re-run is the data the
paper reports**, and it is deposited under `results/renormalized_2026_10/`.

Four strategies are deterministic and were NOT affected; their numbers are
bit-identical before and after: `l1_unstructured`, `obs_compensated`,
`random_unstructured`, `simulation_noise_prune_mask_only`.

## Files that still contain pre-correction rows

These tables mix affected and unaffected rows, so they are kept whole rather than
rewritten. Filter out the six strategies above and take them from the re-run
instead; the remaining rows are unchanged and still valid.

| file | pre-correction rows | total rows |
|---|---|---|
| `results/amplification_factors/rescale_value_histogram_snp_lnp_p80_by_run.csv` | 7517 | 15071 |
| `results/amplification_factors/rescale_value_histogram_snp_lnp_p80_by_run_logbin0p20.csv` | 2099 | 4205 |
| `results/amplification_factors/rescale_value_histogram_snp_lnp_p80_overall.csv` | 136 | 267 |
| `results/amplification_factors/rescale_value_histogram_snp_lnp_p80_overall_logbin0p20.csv` | 38 | 77 |
| `results/amplification_factors/rescale_value_run_summary_snp_lnp_p80.csv` | 72 | 144 |
| `results/capped_rescale/capped_rescale_q10_to_q100_3seed_curve_rows.csv` | 2880 | 5760 |
| `results/capped_rescale/capped_rescale_q10_to_q100_3seed_spectral_abscissa_summary.csv` | 40 | 80 |
| `results/capped_rescale/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_lnp_snp_capped_rescale_q10_to_q90_3seed_curve_rows.csv` | 2592 | 5184 |
| `results/capped_rescale/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_lnp_snp_capped_rescale_q10_to_q90_3seed_curve_summary.csv` | 36 | 72 |
| `results/renormalized_2026_10/capq50_h512_p10_80.csv` | 1152 | 1176 |
| `results/renormalized_2026_10/h1024_p10_80.csv` | 576 | 600 |
| `results/renormalized_2026_10/h512_p10_80.csv` | 3456 | 3480 |
| `results/revision_2026/covariance_ablation/reviewer1_controls_combined.csv` | 16 | 52 |
| `results/revision_2026/covariance_ablation/reviewer1_controls_tests.csv` | 12 | 36 |
| `results/revision_2026/covariance_ablation/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_covariance_ablation_full_p50_80_raw.csv` | 288 | 768 |
| `results/revision_2026/covariance_ablation/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_covariance_ablation_full_p50_80_summary_by_sparsity.csv` | 4 | 16 |
| `results/revision_2026/low_sparsity/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_low_sparsity_p10_40_raw.csv` | 864 | 1728 |
| `results/revision_2026/low_sparsity/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_low_sparsity_p10_40_summary.csv` | 24 | 64 |
| `results/revision_2026/reference_tables/floor_corrected_retention_summary.csv` | 8 | 28 |
| `results/revision_2026/reference_tables/recurrent_weight_retention.csv` | 8 | 28 |
| `results/revision_2026/reference_tables/retention_above_one_counts.csv` | 8 | 28 |
| `results/revision_2026/score_controls/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_full_p50_80_holm_tests.csv` | 30 | 75 |
| `results/revision_2026/score_controls/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_full_p50_80_plot_points.csv` | 3 | 11 |
| `results/revision_2026/score_controls/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_full_p50_80_raw.csv` | 864 | 2592 |
| `results/revision_2026/score_controls/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_full_p50_80_summary_by_sparsity.csv` | 12 | 44 |
| `results/revision_2026/score_controls_round2/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_round2_p50_80_gain_audit.csv` | 64 | 320 |
| `results/revision_2026/score_controls_round2/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_round2_p50_80_raw.csv` | 64 | 320 |
| `results/revision_2026/score_controls_round2/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_round2_p50_80_regression_check.csv` | 32 | 96 |
| `results/revision_2026/score_controls_round2/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_score_control_round2_p50_80_summary_by_sparsity.csv` | 8 | 40 |
| `results/revision_2026/statistics/paper_relevant_h512_descriptive_statistics.csv` | 256 | 576 |
| `results/revision_2026/statistics/per_task_retention_curves.csv` | 64 | 224 |
| `results/statistics/paper_relevant_h512_descriptive_statistics.csv` | 256 | 576 |
| `results/task_retention/revised_task_preservation_h512_24net_p50_80.csv` | 576 | 1464 |
| `results/task_retention/revised_task_preservation_h512_24net_p50_80_spectral_abscissa_summary.csv` | 8 | 28 |
| `results/task_retention/revised_task_preservation_h512_24net_p50_80_summary.csv` | 8 | 28 |
| `results/weight_distributions/recurrent_weight_distribution_summary_p80_by_method.csv` | 1 | 3 |
| `results/weight_distributions/recurrent_weight_distribution_summary_p80_by_run.csv` | 72 | 168 |
| `results/weight_distributions/recurrent_weight_existing_metric_summary_p80_all_methods.csv` | 144 | 360 |
| `results/weight_distributions/recurrent_weight_log10_abs_histogram_p80_by_run.csv` | 5760 | 13440 |
| `results/weight_distributions/recurrent_weight_log10_abs_histogram_p80_overall.csv` | 80 | 240 |
| `results/weight_distributions/recurrent_weight_signed_histogram_p80_by_run.csv` | 14400 | 33600 |
| `results/weight_distributions/recurrent_weight_signed_histogram_p80_overall.csv` | 200 | 600 |

## Affected but NOT re-run

These are reviewer-control ablations of the sampling rule -- a sigma sweep and a
zero-task-input condition -- that the re-run did not cover. Every row in them
predates the correction and nothing replaces them, so they are kept and labelled
rather than deleted. Read them as evidence about the OLD rule.

  - `results/revision_2026/noise_necessity/`
  - `results/revision_2026/noise_only/`

## Files removed

  (none -- see the note above)

These statistical outputs are superseded by the change of PROCEDURE, not only of
data: the paper now reports a hierarchical bootstrap over tasks with an exact
sign-flip test, in place of Wilcoxon signed-rank with Holm correction. The
replacements are in `results/statistics_2026_10/`:

  - `results/statistics/main_task_retention_84_wilcoxon_holm_sign_tests.csv`
  - `results/statistics/paper_relevant_h512_comparison_tests.csv`
  - `results/statistics/revised_scope_task_only_significance_tests.csv`
  - `results/statistics/revised_scope_task_only_significance_means.csv`

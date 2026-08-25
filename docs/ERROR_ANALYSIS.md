# Error Analysis

Status: EXECUTED AND VERIFIED for the fast FireProtDB Ridge baseline.

Input:

- `results/fireprotdb_fast_baselines/fast_baseline_out_of_fold_predictions.csv`

Generated outputs:

- `results/fireprotdb_fast_baselines/ridge_fast_largest_residuals.csv`
- `results/fireprotdb_fast_baselines/ridge_fast_metrics_by_mutant_aa.csv`
- `results/fireprotdb_fast_baselines/ridge_fast_cys_out_of_fold_predictions.csv`
- `results/fireprotdb_fast_baselines/figures/ridge_fast_predicted_vs_observed.png`
- `results/fireprotdb_fast_baselines/figures/ridge_fast_residual_distribution.png`
- `results/fireprotdb_fast_baselines/figures/ridge_fast_mae_by_mutant_aa.png`

## Largest Residuals

The largest absolute residuals are extreme DDG measurements that the simple physicochemical Ridge model predicts near the central range. Examples include:

- Cytosolic beta-glucosidase `N391A`: observed `23.21`, predicted `-0.74`, absolute residual `23.95` kcal/mol.
- Divalent-cation tolerance protein CutA `S11V`: observed `-22.39`, predicted `-0.74`, absolute residual `21.65` kcal/mol.
- Tail spike protein `R383S`: observed `17.40`, predicted `-1.19`, absolute residual `18.59` kcal/mol.

These are not necessarily "bad data"; they show that a mutation-only linear model cannot capture all structural, protein-specific, and experimental-condition effects.

## Cys-Specific Residuals

For X->Cys mutations:

- `n = 25,026`
- MAE `0.731 kcal/mol`
- RMSE `0.968 kcal/mol`
- R2 `0.110`
- Pearson `0.336`
- Spearman `0.301`

This is better than the Dummy baseline but still modest. It should be described as an initial baseline, not a mature protein-engineering predictor.

## By Mutant Amino Acid

The best MAE among mutant residue categories was for Cys (`0.731 kcal/mol`). The worst was Pro (`0.973 kcal/mol`). Proline substitutions are often structurally special because they constrain backbone geometry, but this analysis is associative, not causal.

## Remaining Analyses

Still pending:

- structural residual analysis for buried/exposed residues;
- protein-family or homology-cluster residuals;
- condition-aware residual inspection;
- Random Forest/Gradient Boosting residual comparison;
- external benchmark residual analysis.

No biological causal claims should be made from residuals alone.

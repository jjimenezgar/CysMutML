# Structural Error Analysis

Status: PARTIALLY EXECUTED on 114 mapped structural rows.

The ablation used identical rows and frozen protein-grouped folds for physicochemical-only and physicochemical-plus-structure feature sets.

Largest residuals are available in:

`results/structural_ablation/out_of_fold_predictions.csv`

The paired comparison shows where structural features improved or worsened predictions fold-by-fold:

`results/structural_ablation/paired_ablation_comparison.csv`

On this small subset, adding structure worsened average MAE for both Ridge and HistGradientBoosting. This should be interpreted as a preliminary dataset/mapping-stage result, not a biological claim that structure is unimportant.

Pending analyses require a larger high-confidence structural dataset:

- buried versus exposed mutations;
- Pro/Gly mutations;
- charge reversals;
- Cys-specific effects beyond 6 observations;
- proteins with systematic high error;
- mutations where structure substantially improves or worsens prediction.

# CysMutML v1.0 Model Status Report

Date: 2026-08-25

## Executive Summary

CysMutML v1.0 has a real, deployed physicochemical model for mutation-associated destabilization. The model is scientifically modest, executable, and evaluated on the FireProtDB mutation-level median-aggregated dataset using protein-grouped cross-validation.

Current deployed model:

```text
Ridge regression using aggregated physicochemical mutation features
```

This is the production ML component of the v1.0 hybrid architecture. Structural features are used downstream in the engineering heuristic, not as ML training features.

## Current Model Status

Status:

```text
IMPLEMENTED AND EXECUTED
```

Model artifact:

```text
models/cysmutml_model.joblib
```

Metadata:

```text
models/model_metadata.json
```

Model class:

```text
Ridge
```

Training dataset:

```text
FireProtDB v2.0 API CSV export, aggregated to one row per protein/mutation using median DDG
```

Feature configuration:

```text
physicochemical
```

Target:

```text
destabilization_ddg_kcal_mol
```

Target convention:

```text
larger positive values = greater destabilization
```

## Dataset Used

FireProtDB was downloaded through the public documented API and processed locally.

Processed dataset summary:

| Quantity | Value |
|---|---:|
| Raw FireProtDB rows | 613,208 |
| Processed valid single substitutions | 555,932 |
| Unique proteins | 542 |
| Unique protein/mutation pairs | 352,005 |
| Measurement-level X->Cys rows | 25,026 |
| Aggregated X->Cys rows | 16,236 |
| Rows with missing PDB field | 547,955 |
| Duplicate protein/mutation measurements flagged | 203,927 |

Data files:

```text
data/raw/fireprotdb.csv
data/processed/fireprotdb_mutations.csv
data/processed/fireprotdb_features.csv
data/processed/fireprotdb_mutations_aggregated.csv
data/processed/fireprotdb_aggregated_features.csv
```

## Current Features

The deployed v1.0 model uses mutation-level physicochemical features only.

Feature examples:

```text
wt_aa
mut_aa
wt_hydrophobicity
wt_volume
wt_mass
wt_charge
wt_polarity
wt_aromatic
mut_hydrophobicity
mut_volume
mut_mass
mut_charge
mut_polarity
mut_aromatic
delta_hydrophobicity
delta_volume
delta_mass
delta_charge
delta_polarity
delta_aromatic
blosum62
```

Delta convention:

```text
delta_property = mutant_property - WT_property
```

Structural features are implemented in code but are intentionally not part of the deployed FireProtDB model. The project now uses structure only on the target PDB as a ranking heuristic.

## Aggregated Physicochemical Model Comparison

Executed command:

```bash
.venv/bin/cysmutml compare-physchem \
  --features data/processed/fireprotdb_aggregated_features.csv \
  --results-dir results/physchem_model_comparison \
  --models dummy_mean,ridge,hist_gradient_boosting
```

Important audit note:

An initial comparison produced unrealistically good results because target-derived aggregate columns were still eligible as features. Those leakage columns were excluded and the comparison was rerun. The valid results are below.

Overall mean metrics:

| Model | MAE | RMSE | R2 | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| Dummy mean | 0.800 | 1.049 | -0.002 | undefined | undefined |
| Ridge | 0.684 | 0.930 | 0.213 | 0.462 | 0.450 |
| HistGradientBoosting | 0.669 | 0.917 | 0.234 | 0.485 | 0.475 |

Cys-specific mean metrics:

| Model | MAE | RMSE | R2 | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| Dummy mean | 0.739 | 0.923 | -0.172 | undefined | undefined |
| Ridge | 0.587 | 0.803 | 0.115 | 0.348 | 0.345 |
| HistGradientBoosting | 0.579 | 0.795 | 0.131 | 0.364 | 0.362 |

Decision:

Ridge remains deployed. HGB is slightly better by grouped CV, but the improvement is modest and Ridge is simpler, faster, and more interpretable.

## Current Validation Strategy

Validation uses protein-grouped cross-validation:

```text
3-fold GroupKFold by protein_id
```

This prevents mutations from the same protein from being split between training and validation folds.

Important:

```text
No random mutation-level split is used as the primary validation strategy.
```

## Executed Model Evaluation

Executed command:

```bash
.venv/bin/cysmutml evaluate-fast \
  --features data/processed/fireprotdb_features.csv \
  --results-dir results/fireprotdb_fast_baselines
```

Then rerun with:

```bash
.venv/bin/cysmutml evaluate-fast \
  --features data/processed/fireprotdb_features.csv \
  --results-dir results/fireprotdb_fast_baselines \
  --save-oof
```

Evaluated models:

```text
dummy_mean
ridge_fast
```

`ridge_fast` is a fast full-dataset Ridge-style linear baseline used for efficient grouped cross-validation.

## Overall Regression Metrics

Mean across 3 protein-grouped folds:

| Model | MAE | RMSE | R2 | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| Dummy mean | 0.949 | 1.193 | -0.003 | null | null |
| Ridge fast | 0.852 | 1.099 | 0.150 | 0.392 | 0.380 |

Interpretation:

- Ridge improves over the Dummy baseline.
- The improvement is real but modest.
- The model captures some mutation-level signal from physicochemical descriptors.
- The remaining error is substantial, which is expected for a linear model without verified structural context.

## Cys-Specific Evaluation

The downstream application focuses on X->Cys mutations, so Cys-only performance was evaluated separately.

Mean across 3 protein-grouped folds:

| Model | MAE | RMSE | R2 | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| Dummy mean | 0.906 | 1.102 | -0.156 | null | null |
| Ridge fast | 0.731 | 0.968 | 0.110 | 0.337 | 0.306 |

Interpretation:

- Ridge also improves over Dummy for X->Cys mutations.
- The Cys-only subset contains 25,026 rows, so this is a meaningful first baseline.
- Performance is still modest and should not be presented as a mature engineering predictor.

## Out-of-Fold Predictions

Out-of-fold predictions were generated for real FireProtDB rows:

```text
results/fireprotdb_fast_baselines/fast_baseline_out_of_fold_predictions.csv
```

This file includes:

```text
model
fold
protein_id
mutation
wt_aa
mut_aa
observed
predicted
residual
absolute_residual
```

These predictions support honest residual and error analysis.

## Error Analysis

Executed command:

```bash
.venv/bin/cysmutml error-analysis \
  --predictions results/fireprotdb_fast_baselines/fast_baseline_out_of_fold_predictions.csv \
  --output-dir results/fireprotdb_fast_baselines \
  --model ridge_fast
```

Generated files:

```text
results/fireprotdb_fast_baselines/ridge_fast_largest_residuals.csv
results/fireprotdb_fast_baselines/ridge_fast_metrics_by_mutant_aa.csv
results/fireprotdb_fast_baselines/ridge_fast_cys_out_of_fold_predictions.csv
```

Generated figures:

```text
results/fireprotdb_fast_baselines/figures/ridge_fast_predicted_vs_observed.png
results/fireprotdb_fast_baselines/figures/ridge_fast_residual_distribution.png
results/fireprotdb_fast_baselines/figures/ridge_fast_mae_by_mutant_aa.png
```

Largest observed residual examples:

| Protein | Mutation | Observed | Predicted | Absolute residual |
|---|---:|---:|---:|---:|
| Cytosolic beta-glucosidase | N391A | 23.21 | -0.74 | 23.95 |
| Divalent-cation tolerance protein CutA | S11V | -22.39 | -0.74 | 21.65 |
| Tail spike protein | R383S | 17.40 | -1.19 | 18.59 |

Interpretation:

- The largest errors occur for extreme DDG values.
- A simple mutation-property linear model predicts these near the central range.
- These errors highlight the need for structural context, protein identity/homology controls, and careful condition metadata analysis.
- These residuals do not prove that the data are wrong.

## Mutant Amino-Acid Error Pattern

From `ridge_fast_metrics_by_mutant_aa.csv`:

- Best MAE among mutant categories: Cys, approximately `0.731 kcal/mol`.
- Worst MAE among mutant categories: Pro, approximately `0.973 kcal/mol`.

This is plausible because proline mutations often have special backbone effects, but this is only an association in model error, not a causal conclusion.

## Real Example Prediction

The real Ridge model was used to run X->Cys prediction on the example PDB:

Input:

```text
examples/tiny_protein.pdb
```

Outputs:

```text
examples/output_real/mutation_predictions.csv
examples/output_real/residue_ranking.csv
examples/output_real/visualize_rankings.pml
```

Important:

The example PDB is tiny and artificial. It verifies inference mechanics with the real model, but it is not a biologically meaningful case study.

## What the Model Can Do Now

The current real model can:

- load a PDB;
- generate X->Cys candidates;
- calculate mutation physicochemical features;
- calculate structural features for the input PDB;
- use the trained Ridge model to predict destabilization;
- produce `mutation_predictions.csv`;
- produce a separate engineering ranking using the heuristic ranking layer.

## What the Model Cannot Do Yet

The current real model cannot yet:

- claim structural-feature improvement over the physicochemical baseline;
- use verified FireProtDB structural context at scale;
- provide calibrated uncertainty;
- predict immobilization success;
- guarantee activity retention;
- replace experimental validation;
- claim state-of-the-art protein stability prediction.

## Main Scientific Limitations

The current model is limited by:

- heterogeneous experimental conditions;
- database bias;
- duplicate and repeated measurements;
- protein-family imbalance;
- missing reliable bulk residue-to-structure mapping;
- no verified structure features in the real FireProtDB model yet;
- static structures not representing conformational ensembles;
- physicochemical descriptors being too simple for many structural effects;
- stability not being equivalent to activity;
- stability not being equivalent to immobilization suitability.

## Current Engineering Ranking Status

The engineering ranking is implemented and executable.

It is deliberately separate from the ML prediction.

Ranking output:

```text
residue_ranking.csv
```

Ranking components:

```text
stability_component
exposure_component
protected_site_component
existing_cys_proximity_warning
cys_suitability_score
```

The ranking weights are heuristic and configurable. They are not experimentally optimized.

## Verification Status

Final checks executed:

```text
pytest -q
ruff check .
python -m json.tool reports/run_summary.json
```

Status:

```text
pytest: 8 passed
ruff: all checks passed
run_summary.json: valid JSON
```

## Current Bottom Line

CysMutML now has a real, evaluated, reproducible baseline model:

```text
FireProtDB -> physicochemical features -> GroupKFold Ridge baseline -> real metrics -> real OOF error analysis -> real model artifact -> X->Cys inference
```

This is a scientifically honest first model. The next major scientific step is to add verified structural features at scale and test whether they improve over this baseline.

## Recommended Next Steps

1. Build a structure-mapping subset from FireProtDB rows with reliable PDB and chain metadata.
2. Extract structural features only where WT identity can be verified.
3. Train/evaluate the same Ridge baseline on:
   - physicochemical-only features;
   - physicochemical plus structural features.
4. Run a real ablation study.
5. Add permutation importance for the Ridge/linear baseline and any tree-based model that completes.
6. Optimize Random Forest/Gradient Boosting evaluation or use a protein-grouped controlled subset.
7. Investigate S669 only after overlap checking.

# CysMutML Model Card

## Model overview

CysMutML v1.0 deploys an interpretable Ridge regression baseline that predicts
mutation-associated destabilization. It is one component of a hybrid system:
the learned stability estimate and the deterministic target-PDB ranking heuristic
remain deliberately separate.

| Item | Value |
|---|---|
| Model | Ridge regression |
| Training source | FireProtDB v2.0 API CSV export |
| Training rows | 352,005 median-aggregated mutation records |
| Protein groups | 542 |
| X→Cys rows | 16,236 |
| Target | `destabilization_ddg_kcal_mol` |
| Sign convention | Larger positive values mean greater destabilization |
| Validation | 3-fold GroupKFold by `protein_id` |
| Artifact | `models/cysmutml_model.joblib` |

## Intended use

The model provides a lightweight, reproducible mutation-tolerance signal for
prioritizing candidate X→Cys substitutions. It is suitable for ranking and
hypothesis generation, not for replacing experimental measurements.

## Features

The deployed model uses amino-acid physicochemical descriptors:

- wild-type and mutant hydrophobicity, volume, mass, charge, polarity, and aromaticity;
- mutant-minus-wild-type property deltas;
- BLOSUM62 substitution score;
- one-hot encoded wild-type and mutant amino-acid identities.

Structural descriptors are not training features in v1.0. SASA, B-factor-derived
flexibility, protected-site distance, local exposed Lys context, and native-Cys
context belong to the separate engineering heuristic.

## Performance

Mean metrics across protein-grouped folds:

| Population | Model | MAE | RMSE | R² | Pearson | Spearman |
|---|---|---:|---:|---:|---:|---:|
| All mutations | Dummy mean | 0.800 | 1.049 | -0.002 | — | — |
| All mutations | Ridge | 0.684 | 0.930 | 0.213 | 0.462 | 0.450 |
| All mutations | HistGradientBoosting | 0.669 | 0.917 | 0.234 | 0.485 | 0.475 |
| X→Cys | Dummy mean | 0.739 | 0.923 | -0.172 | — | — |
| X→Cys | Ridge | 0.587 | 0.803 | 0.115 | 0.348 | 0.345 |
| X→Cys | HistGradientBoosting | 0.579 | 0.795 | 0.131 | 0.364 | 0.362 |

HGB performs slightly better. Ridge remains deployed because the gain is small
relative to the interpretability and operational simplicity of the linear model.

## Evaluation design

Primary evaluation groups rows by protein, so mutations from the same protein do
not appear in both train and test folds. Random mutation-level splitting is not
used as the headline result.

This design does not guarantee separation of homologous proteins. CysMutML v1.2
adds an MMseqs2-based sequence-cluster split to estimate performance on less-related
protein families. The infrastructure is CI-tested; numerical results are not claimed
until the FireProtDB tables are regenerated and the experiment is executed.

## Known limitations

- The same substitution can receive the same ML prediction at different positions
  because v1.0 has no position-specific sequence or structural context.
- FireProtDB combines heterogeneous proteins, assays, temperatures, pH values, and methods.
- Stability is not activity, cysteine reactivity, immobilization yield, or retained activity.
- The model does not provide calibrated predictive uncertainty.
- The final engineering score uses documented heuristic weights and is not a probability.
- Performance is modest and must not be described as state of the art.

## Ethical and scientific use

CysMutML is a research prioritization tool. Predictions should be accompanied by
domain review and experimental validation. Negative or positive scores must not be
presented as proof that a mutation will fail or succeed.

## Reproducibility

The repository contains:

- exact package dependencies and a command-line interface;
- grouped cross-validation results;
- model metadata and a serialized artifact;
- self-contained tests;
- a CI workflow for lint, tests, package build, and portfolio-notebook execution;
- a retrospective validation with explicitly frozen heuristic settings.

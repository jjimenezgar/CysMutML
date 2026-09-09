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
| Training rows | 351,487 median-aggregated mutation records |
| Evaluation groups | 543 |
| X→Cys rows | 16,208 |
| Target | `destabilization_ddg_kcal_mol` |
| Sign convention | Larger positive values mean greater destabilization |
| Validation | 3-fold GroupKFold by `protein_id` |
| Artifact | `models/cysmutml_model.joblib` |

## MVP status

Portfolio MVP: a reproducible ML workflow for cysteine-site prioritization. Predictions are exploratory; reliable transfer to new enzymes has not been established.

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
flexibility, optional user-defined protected-site annotations, local exposed Lys
context, and native-Cys context belong to the separate engineering heuristic.

## Performance

Mean metrics across protein-grouped folds:

| Population | Model | MAE | RMSE | R² | Pearson | Spearman |
|---|---|---:|---:|---:|---:|---:|
| All mutations | Dummy mean | 0.792 | 1.040 | -0.009 | — | — |
| All mutations | Ridge | 0.667 | 0.912 | 0.223 | 0.481 | 0.467 |
| All mutations | HistGradientBoosting | 0.662 | 0.908 | 0.231 | 0.489 | 0.479 |
| X→Cys | Dummy mean | 0.731 | 0.913 | -0.189 | — | — |
| X→Cys | Ridge | 0.576 | 0.790 | 0.112 | 0.344 | 0.313 |
| X→Cys | HistGradientBoosting | 0.571 | 0.783 | 0.128 | 0.367 | 0.363 |

HGB performs slightly better. Ridge remains deployed because the gain is small
relative to the interpretability and operational simplicity of the linear model.

### Homology-aware MVP (historical, pre-audit)

A reduced, deterministic benchmark was executed on 150 proteins and 5,634 mutation
rows after MMseqs2 clustering at 30% identity and 80% coverage. Mean MAE across three
folds was 1.508 for Ridge with protein grouping and 1.523 with homology-cluster
grouping. On X→Cys rows, the corresponding values were 1.535 and 1.630. The
benchmark also includes Dummy, Random Forest, and HistGradientBoosting, fit/predict
timings, and held-out permutation importance. These results demonstrate validation
discipline and residual homology effects; they are not a claim of state-of-the-art
accuracy.


## Evaluation design

Primary evaluation groups rows by protein, so mutations from the same protein do
not appear in both train and test folds. Random mutation-level splitting is not
used as the headline result.

This design does not guarantee separation of homologous proteins. The separate
[ESM experiment](docs/ESM_CONTEXT_COMPARISON.md) evaluates corrected descriptors and
frozen embeddings on identical homology-grouped folds for 5,506 verified rows.
Its results did not justify promoting ESM. They are not directly comparable to the
full-data metrics above.

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

- package dependency requirements (not a locked training environment) and a command-line interface;
- grouped cross-validation results;
- model metadata and a serialized artifact;
- self-contained tests;
- a CI workflow for lint, tests, package build, and portfolio-notebook execution;
- a retrospective validation with explicitly frozen heuristic settings;
- an executed, reduced homology-aware MVP with versioned fold metrics and sampling audit.


## Audit status

The audit branch corrects the canonical BLOSUM62 descriptor and rejects composite mutation labels during ingestion. The committed Ridge artifact and full physicochemical benchmarks were regenerated in Actions run 34040897260 on 6 September 2026. Historical homology results were not regenerated and are not evidence for this corrected artifact. The retrospective Godoy report is retained for traceability, but BTL2 rows with unresolved residue joins are not valid validation evidence. AlphaFold pLDDT is treated as confidence rather than experimental mobility in inference.

## Data audit and scope

The [sequence recovery audit](docs/SEQUENCE_RECOVERY_AUDIT.md) found 68 conflicting WT positions across 29 protein names in records lacking sequence. Construct identity, numbering, source-value aggregation and sign provenance need reconciliation before stronger predictive claims. This does not invalidate every row, but limits interpretation of the broader benchmark. The MVP retains Ridge as a demonstration baseline; source reconstruction and further retraining are deferred. ESM remains an offline experiment, not part of app inference.

# CysMutML Interview Guide

## 30-Second Explanation

CysMutML is a hybrid protein-engineering ML project. I trained a Ridge regression model on FireProtDB mutation-stability data using interpretable physicochemical mutation features, evaluated it with protein-grouped cross-validation to avoid leakage, and then combined its X->Cys destabilization predictions with target-PDB structural heuristics like SASA and B-factor-derived rigidity to rank cysteine mutation candidates.

## 2-Minute Explanation

The project separates two problems that are often mixed together. The ML model learns mutation-associated destabilization from experimental stability data. It uses WT amino-acid properties, mutant properties, property deltas, and BLOSUM62. It does not use structural features because only a smaller, biased subset of FireProtDB can be mapped confidently to PDB structures.

For a new PDB, the pipeline generates all X->Cys candidates, predicts destabilization with the ML model, and separately calculates target-specific structural descriptors: relative SASA, a B-factor-derived rigidity proxy, optional protected-site distance, and existing-cysteine proximity. A transparent heuristic combines those components into a final ranking. The score is not a probability and not an immobilization predictor.

## Why Ridge Instead of HGB?

HistGradientBoosting was slightly better by grouped CV, but the improvement was modest. Ridge was close, faster, easier to reproduce, and easier to explain. For a v1.0 portfolio project emphasizing scientific transparency, Ridge was the better deployed model.

## Why GroupKFold?

Mutations from the same protein are correlated. If mutations from one protein appear in both train and test sets, the model can exploit protein-specific patterns and produce optimistic metrics. GroupKFold keeps all mutations from a protein in the same fold.

## Why Aggregate Repeated Measurements?

FireProtDB contains repeated measurements for the same protein/mutation. Training on every measurement can overweight heavily studied mutations. Median aggregation gives one robust central value per protein/mutation while preserving measurement count and dispersion metadata.

## Why Not Train With SASA?

Training with SASA would require reliable residue-to-structure mapping for many mutations. That substantially shrinks and biases the dataset. CysMutML v1.0 instead learns from the large mutation dataset and applies structure only to the target PDB during ranking.

## Why Use SASA in the Heuristic?

For cysteine engineering, accessibility matters. A mutation predicted to be stable but buried may be less practical for labeling or attachment. SASA is a simple, inspectable way to prioritize exposed sites.

## Why Is B-Factor Only a Proxy?

B-factors are affected by crystallographic refinement, resolution, occupancy, disorder, and local motion. They are useful as a lightweight structural signal, but they are not a direct measurement of solution dynamics.

## Why Is the Score Not a Probability?

The final score combines an ML regression output with heuristic structural components and hand-chosen weights. It was not trained or calibrated against experimental cysteine-engineering success data.

## Main Model Limitation

The model uses compact physicochemical descriptors, so it cannot fully capture detailed structural context, conformational changes, activity effects, or protein-family-specific behavior.

## Highest-Value Improvement

The best next scientific extension would be an external benchmark and overlap audit, such as S669 if it can be obtained reproducibly and kept independent from training.

## ML Concepts Demonstrated

- data cleaning and target normalization;
- duplicate handling and aggregation;
- feature engineering;
- regression;
- baseline comparison;
- leakage-safe grouped cross-validation;
- model selection;
- Cys-specific subgroup evaluation;
- interpretability;
- error analysis;
- reproducible inference pipeline.

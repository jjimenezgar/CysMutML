# CysMutML Project Plan

## 1. Scientific Objective
CysMutML predicts experimentally measured destabilization associated with single-point amino-acid substitutions, then separately ranks candidate X->Cys mutations for protein-engineering use.

## 2. Prediction Problem
Primary ML task: regression from mutation plus wild-type structural context to `destabilization_ddg_kcal_mol`.

Downstream task: heuristic ranking of candidate cysteine substitutions from a new PDB structure. This is not a supervised ML target.

## 3. Intended Dataset
Primary source: FireProtDB v2.0, accessed through its documented REST API or user-exported CSV/TSV placed in `data/raw/`.

The FireProtDB API documentation states `/api/search` supports `json`, `jsonl`, `csv`, and `tsv` formats. Its help page defines DDG sign convention as negative stabilizing and positive destabilizing.

Optional external benchmark: S669 may be supported only when a reproducible legal copy is provided locally and overlap against training entries is checked.

## 4. Feature Definitions
Feature groups:

- Model A: wild-type identity/properties, mutant identity/properties, and mutant-minus-WT delta properties.
- Model B: Model A plus structural features: SASA, relative SASA, CA contacts at 6/8/10 A, heavy-atom contacts, normalized CA distance to protein center, local density, residue B-factor and chain-normalized B-factor.
- Model C: Model B plus optional context such as secondary structure when a reliable dependency is available.

## 5. Target Definition
Internal canonical target:

`destabilization_ddg_kcal_mol`

Convention: larger positive values mean greater destabilization.

FireProtDB DDG is documented as negative = stabilizing and positive = destabilizing, so for FireProtDB:

`destabilization_ddg_kcal_mol = source_ddg_kcal_mol`

If another source uses the opposite sign, a source-specific transformer must be added and documented before use.

## 6. Data-Cleaning Strategy
Keep only:

- single substitutions;
- canonical one-letter amino acids;
- interpretable DDG values in kcal/mol;
- rows with protein identifiers;
- rows passing WT identity checks when structure mapping is attempted.

Exclude or flag insertions, deletions, multi-mutants, noncanonical residues, ambiguous target semantics, unresolved mapping, missing residues, and incompatible units.

## 7. Protein/Mutation Mapping Strategy
Never assume database residue position equals PDB residue number. Mapping aligns canonical sequence to selected PDB chain sequence, maps sequence positions to residues, and verifies dataset WT equals mapped structure WT. Failures are written to `reports/residue_mapping_report.csv`.

## 8. Train/Validation/Test Strategy
Primary evaluation uses `GroupKFold` grouped by normalized `protein_id`. A random mutation split is not a primary validation strategy.

## 9. Leakage Prevention Strategy
All preprocessing lives inside scikit-learn `Pipeline`/`ColumnTransformer` and is fit only on training folds. Proteins are not split across folds.

## 10. Model-Selection Strategy
Compare DummyRegressor, Ridge, RandomForestRegressor, and HistGradientBoostingRegressor. Hyperparameter tuning is deliberately small and performed within grouped cross-validation.

## 11. Metrics
Regression: MAE, RMSE, R2, Pearson correlation, Spearman correlation.

Derived classification: ROC-AUC, PR-AUC, balanced accuracy, precision, recall, F1, MCC, confusion matrix at configurable threshold.

## 12. Inference Strategy
For a supplied PDB and chain, generate every canonical non-Cys residue as X->C, extract the same feature schema, predict destabilization, and emit `mutation_predictions.csv`.

## 13. Cysteine Ranking Strategy
The engineering ranking is a transparent heuristic using predicted destabilization, relative SASA/exposure, optional protected-residue penalties, and existing-cysteine proximity warnings. It is not an ML prediction.

## 14. Assumptions
FireProtDB DDG values are in kcal/mol and use the documented sign convention. Structure-derived features are computed from static coordinates and representative chain models. For v1, uncertain mappings are excluded instead of guessed.

## 15. Scientific Limitations
Experimental conditions vary; databases are biased; protein families are imbalanced; PDB structures may miss residues or differ from solution ensembles; SASA/contact features are simplified descriptors; stability is not activity or immobilization performance; Cys-specific generalization may differ from overall performance.

## 16. Structural-Ablation Milestone
The next milestone tests whether explicit structural context improves mutation-associated destabilization prediction beyond physicochemical mutation descriptors alone.

Required order:

1. audit duplicate FireProtDB measurements;
2. create a mutation-level aggregated dataset using median DDG;
3. audit extreme DDG values without deleting them from the main analysis;
4. identify FireProtDB records with PDB and UniProt identifiers;
5. download and validate experimentally resolved PDB structures into `data/structures/`;
6. retrieve canonical UniProt sequences where available;
7. align canonical sequences to candidate PDB chains;
8. accept mappings only when WT identity is verified;
9. extract a small set of interpretable structural features;
10. freeze deterministic protein-grouped folds on the high-confidence structural subset;
11. compare physicochemical-only versus physicochemical-plus-structure models on exactly the same rows and folds.

The deployed physicochemical Ridge model must not be overwritten during this milestone. If a structural model performs better, it should first be saved as a candidate artifact and promoted only after documenting the evidence.

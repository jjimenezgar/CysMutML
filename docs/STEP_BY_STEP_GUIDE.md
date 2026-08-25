# Step-by-Step Guide

## Current Hybrid Architecture

CysMutML now has two deliberately separate parts:

1. A learned physicochemical ML model predicts `predicted_destabilization_ddg`.
2. A target-PDB structural heuristic calculates accessibility, rigidity, protected-site penalties, and existing-cysteine warnings.

The final `cys_suitability_score` is not a learned probability. It is a transparent weighted score.

## Executed Example: 1CSP F38C

In the real case study `examples/real_case/`, residue A:38 in PDB `1csp` is phenylalanine.

1. Virtual mutation:

```text
F38C
```

2. Physicochemical descriptors:

The model encodes WT Phe properties, mutant Cys properties, mutant-minus-WT deltas, and BLOSUM62.

3. ML prediction:

```text
predicted_destabilization_ddg = -1.040 kcal/mol
```

Larger positive values mean more predicted destabilization, so this value is favorable by the current model.

4. Stability normalization:

```text
stability_component = 1.000
```

This is a clipped score, not a calibrated probability.

5. Structure:

```text
relative_sasa = 0.562
flexibility_method = BFACTOR
flexibility_value = -0.766
```

6. Structural normalization:

```text
accessibility_component = 0.562
rigidity_component = 0.832
```

7. Penalties:

```text
protected_site_penalty = 0.000
existing_cys_penalty = 0.000
```

8. Final ranking:

```text
cys_suitability_score = 0.835
rank_engineering = 1
```

This means F38C is top-ranked by the current heuristic for this example, not experimentally validated.

## 1. What problem are we solving?
We want to predict how destabilizing a single amino-acid substitution is, then use that prediction to help rank possible cysteine mutations for protein engineering.

## 2. What does one row represent?
One row is one experimental measurement for one single substitution, such as `K42C`, in one protein.

## 3. What is X?
`X` is the feature vector: WT and mutant amino-acid properties, their deltas, and, when verified mapping is available, WT structural context.

## 4. What is y?
`y` is `destabilization_ddg_kcal_mol`.

## 5. What does DDG mean?
In this project, larger positive DDG means greater destabilization. FireProtDB already uses positive = destabilizing and negative = stabilizing, so no sign flip is applied.

## 6. What preprocessing was applied?
The implemented FireProtDB parser keeps single canonical substitutions with numeric DDG and preserves method, pH, temperature, source dataset, original DDG value, units, and sign convention.

## 7. How is a mutation mapped onto a PDB?
The implemented mapper aligns a canonical sequence against the selected PDB chain sequence, maps the dataset position to a structural residue, and rejects the mapping unless the WT amino acid matches.

## 8. What structural features are calculated?
Absolute SASA: solvent-accessible surface area from Biopython Shrake-Rupley, in square angstroms.

Relative SASA: absolute SASA divided by Tien 2013 maximum ASA for that residue type.

CA contacts: number of C-alpha atoms within 6, 8, and 10 angstroms.

Heavy-atom contacts: number of target residue heavy atoms near other heavy atoms within 4.5 angstroms.

Distance to center: CA distance to chain geometric center, normalized by the largest chain CA-center distance.

B-factor: mean residue experimental B-factor and chain-normalized B-factor. This is not a pure flexibility measure.

Secondary structure: currently `unknown` unless a future optional DSSP path is added.

## 9. What mutation features are calculated?
Hydrophobicity, volume, mass, charge, polarity, aromaticity for WT and mutant residues, plus deltas computed as mutant minus WT. BLOSUM62 substitution score is also included.

## 10. What is train/validation/test?
Training fits the model. Validation estimates generalization. External test data, such as S669, must remain untouched until final evaluation.

## 11. Why is random mutation splitting dangerous?
If mutations from the same protein are split randomly, the model can learn protein-specific behavior from training rows and appear better on test rows from that same protein.

## 12. What does GroupKFold do?
GroupKFold assigns whole proteins to folds. No protein group appears in both train and validation for a fold.

## 13. What models were tested?
Code supports DummyRegressor, Ridge, Random Forest, and histogram gradient boosting. In this run, fast full-FireProtDB grouped CV completed for Dummy and Ridge-style linear regression. Full Random Forest and gradient boosting CV did not complete interactively.

## 14. What is hyperparameter optimization?
It is selecting model settings without using the final test set. CysMutML keeps this modest; v1 uses conservative defaults.

## 15. What metrics are used?
MAE is average absolute error. RMSE emphasizes larger errors. R2 compares to a mean predictor. Pearson measures linear correlation. Spearman measures rank correlation.

## 16. What is feature importance?
Permutation importance measures how prediction quality changes when a feature is shuffled. It describes model reliance, not causality.

## 17. What is an ablation study?
An ablation compares feature groups, for example physicochemical-only versus physicochemical plus structure, to test whether structure improves generalization.

## 18. What did the model actually learn?
The first real model is a physicochemical Ridge baseline. In grouped FireProtDB CV, it improved mean MAE from 0.949 kcal/mol for the Dummy baseline to 0.852 kcal/mol. For X->Cys rows, mean MAE improved from 0.906 to 0.731 kcal/mol. These are modest baseline results, not state-of-the-art claims.

## 19. What happens when I provide a new PDB?
The CLI parses the selected chain, skips existing Cys residues, computes residue structural features, creates virtual X->Cys mutation features, and predicts destabilization with a trained artifact.

## 20. How is X->Cys generated virtually?
The WT residue stays the structural context. The mutant residue is set to C, and all delta descriptors are computed as Cys minus WT.

## 21. What does predicted destabilization mean?
It is the model-estimated DDG in kcal/mol under the internal convention. It is not experimental certainty.

## 22. How is final residue ranking generated?
The ranking combines normalized low predicted destabilization, high exposure, and optional protected-site penalty components.

## 23. Which part is ML and which part is heuristic?
The DDG prediction is ML. The final cysteine suitability score is a transparent engineering heuristic.

## 24. Main limitations
Experimental conditions vary; databases are biased; mapping is uncertain without sequence/chain metadata; structures are static; stability is not activity; accessibility does not guarantee immobilization success.

## 25. Interview explanation
30 seconds: CysMutML is an end-to-end protein stability ML project that normalizes FireProtDB mutation DDG data, builds interpretable mutation and structural features, uses protein-grouped validation to avoid leakage, and applies the resulting stability model separately from a cysteine-engineering ranking heuristic.

2 minutes: The project treats mutation stability prediction and cysteine site selection as related but separate tasks. It defines a clear DDG sign convention, preserves experimental metadata, prevents protein-level leakage with GroupKFold, implements classical interpretable models, and extracts structural descriptors from PDB coordinates. For new structures it generates X->Cys candidates, predicts destabilization, then ranks candidates using a documented exposure/stability heuristic. The current run downloaded and processed FireProtDB and verified the software pipeline on a toy fixture, but full scientific CV was not completed in this interactive session.

## Structural Model Milestone

Physicochemical descriptors alone describe what substitution happened, but not where it happened in the folded protein. A buried aromatic residue and a solvent-exposed aromatic residue can have very different mutation tolerance even if the WT->mutant descriptor deltas are identical.

Structure might help by adding context:

- SASA measures how exposed a residue is to solvent. High relative SASA usually means surface-exposed; low relative SASA suggests burial.
- Contact number counts nearby residues by C-alpha distance. Higher contact counts are a simple packing proxy.
- B-factor is an experimental crystallographic field. It may reflect disorder, refinement, resolution, occupancy, or flexibility, so CysMutML calls it B-factor, not pure flexibility.
- Sequence-to-structure mapping means aligning the canonical protein sequence to the residues actually present in a PDB chain.
- PDB residue numbering can differ from sequence numbering because of signal peptides, missing residues, engineered constructs, tags, insertion codes, or historical numbering.
- WT identity must match after mapping; otherwise the row is not used for structural training.

The structural subset is the set of FireProtDB rows with PDB/UniProt identifiers, downloaded structures/sequences, an unambiguous chain, verified WT identity, and successfully calculated structural features.

For ablation, both models must use exactly the same rows and folds. Comparing 555,932 physicochemical rows against 114 mapped structural rows would be invalid. The executed structural ablation therefore used the same 114 mapped rows for both feature sets.

Actual result in the partial structural run: adding structural features did not improve performance on the executed subset. Ridge MAE worsened from 2.760 to 3.123 kcal/mol; HistGradientBoosting MAE worsened from 2.197 to 2.385 kcal/mol. The subset is small and biased by interactive runtime constraints, so this is not a final claim that structure is unhelpful. It means the current structural candidate is not ready to replace the deployed physicochemical Ridge baseline.

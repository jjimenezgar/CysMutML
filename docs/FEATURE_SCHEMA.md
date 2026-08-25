# CysMutML v1.0 Feature Schema

The deployed v1.0 model is a Ridge regression pipeline trained on physicochemical mutation features only. Structural descriptors are not ML inputs.

## Target

| Column | Meaning |
|---|---|
| `destabilization_ddg_kcal_mol` | Canonical target. Larger positive values mean greater destabilization. Units are kcal/mol. |

## Categorical Features

| Column | Type | Meaning | Encoding |
|---|---|---|---|
| `wt_aa` | string | Wild-type one-letter amino-acid identity | imputed then one-hot encoded |
| `mut_aa` | string | Mutant one-letter amino-acid identity | imputed then one-hot encoded |

## Numeric Features

| Column | Meaning | Unit / Scale | Transformation |
|---|---|---|---|
| `wt_hydrophobicity` | WT amino-acid hydrophobicity | Kyte-Doolittle-like scale | median imputation, standard scaling |
| `wt_volume` | WT side-chain/residue volume descriptor | Angstrom^3-like tabulated scale | median imputation, standard scaling |
| `wt_mass` | WT residue molecular mass | Da | median imputation, standard scaling |
| `wt_charge` | Approximate WT charge category | -1, 0, +1 | median imputation, standard scaling |
| `wt_polarity` | WT polarity indicator | 0/1 | median imputation, standard scaling |
| `wt_aromatic` | WT aromaticity indicator | 0/1 | median imputation, standard scaling |
| `mut_hydrophobicity` | Mutant amino-acid hydrophobicity | same as WT | median imputation, standard scaling |
| `mut_volume` | Mutant volume descriptor | same as WT | median imputation, standard scaling |
| `mut_mass` | Mutant molecular mass | Da | median imputation, standard scaling |
| `mut_charge` | Approximate mutant charge category | -1, 0, +1 | median imputation, standard scaling |
| `mut_polarity` | Mutant polarity indicator | 0/1 | median imputation, standard scaling |
| `mut_aromatic` | Mutant aromaticity indicator | 0/1 | median imputation, standard scaling |
| `delta_hydrophobicity` | Change in hydrophobicity | mutant - WT | median imputation, standard scaling |
| `delta_volume` | Change in volume | mutant - WT | median imputation, standard scaling |
| `delta_mass` | Change in mass | mutant - WT | median imputation, standard scaling |
| `delta_charge` | Change in approximate charge | mutant - WT | median imputation, standard scaling |
| `delta_polarity` | Change in polarity indicator | mutant - WT | median imputation, standard scaling |
| `delta_aromatic` | Change in aromaticity indicator | mutant - WT | median imputation, standard scaling |
| `blosum62` | BLOSUM62 substitution score | substitution log-odds score | median imputation, standard scaling |

## Delta Convention

```text
delta_property = mutant_property - WT_property
```

For example, K->C has:

```text
delta_charge = charge(C) - charge(K) = 0 - 1 = -1
```

## Training/Inference Consistency

The serialized artifact `models/cysmutml_model.joblib` stores:

- `numeric_features`
- `categorical_features`
- the fitted scikit-learn preprocessing pipeline
- the fitted Ridge model
- training feature ranges for simple extrapolation warnings

The inference code generates the same physicochemical schema for every virtual X->Cys mutation before prediction.

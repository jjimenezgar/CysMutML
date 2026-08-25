# CysMutML v1.0 Scientific Audit

Date: 2026-08-25

## Data

Source:

```text
FireProtDB v2.0 API CSV export
```

Status:

- The dataset source is legitimate and reproducibly downloaded by the project pipeline.
- Raw rows: 613,208.
- Processed valid single-substitution rows: 555,932.
- Mutation-level median-aggregated rows: 352,005.
- Unique proteins in aggregated table: 542.
- Measurement-level X->Cys rows: 25,026.
- Aggregated X->Cys rows: 16,236.

Target:

```text
destabilization_ddg_kcal_mol
```

Convention:

```text
larger positive values = greater destabilization
```

Duplicates:

- Repeated protein/mutation measurements are audited.
- The v1.0 deployed model uses median aggregation to reduce overweighting of repeated measurements.
- Measurement count and dispersion metadata are preserved in the aggregated table.

Limitation:

Experimental conditions remain heterogeneous. The v1.0 model does not mechanistically correct for pH, temperature, buffer, or assay differences.

## ML

Deployed model:

```text
Ridge regression
```

Feature set:

```text
physicochemical only
```

Validation:

```text
3-fold GroupKFold by protein_id
```

Audit findings:

- Leakage-safe grouped validation is used.
- Dummy baseline is included.
- Ridge and HistGradientBoosting are compared fairly on the same aggregated table and folds.
- X->Cys subset performance is reported separately.
- HGB is slightly better, but Ridge remains deployed for simplicity and interpretability.

Leakage check:

Target-derived aggregate columns are excluded from ML features:

- `median_destabilization_ddg`
- `mean_destabilization_ddg`
- `std_destabilization_ddg`
- `min_destabilization_ddg`
- `max_destabilization_ddg`
- `n_measurements`

## Inference

The inference pipeline:

- accepts a PDB and chain;
- generates every non-Cys canonical X->Cys candidate;
- calculates the same physicochemical schema used in training;
- predicts continuous `predicted_destabilization_ddg`;
- writes `mutation_predictions.csv`.

The ML output is clearly labeled as predicted destabilization. It is not a feasibility probability.

## Structural Heuristic

Structural descriptors are calculated only on the target PDB:

- relative SASA by Biopython Shrake-Rupley;
- B-factor-derived rigidity proxy;
- optional protected-residue distance;
- existing-cysteine proximity.

The final `cys_suitability_score` is reconstructable from output columns:

```text
score =
  0.50 * stability_component
+ 0.30 * accessibility_component
+ 0.20 * rigidity_component
- 0.10 * protected_site_penalty
- 0.05 * existing_cys_penalty
```

The weights are heuristic defaults and configurable in `configs/default.yaml`.

Protected residues:

- user supplied only;
- no active sites are inferred automatically;
- distances use C-alpha coordinates across chains.

Existing cysteines:

- proximity is reported as an engineering caution;
- proximity does not imply disulfide formation.

## Claims

CysMutML v1.0 does not claim:

- experimentally validated immobilization success;
- activity retention;
- calibrated probability of success;
- disulfide prediction;
- state-of-the-art DDG prediction.

Supported claim:

```text
CysMutML v1.0 is a reproducible hybrid ML and structural-bioinformatics pipeline for transparent prioritization of candidate cysteine substitutions.
```

## Exploratory Structural ML

The previous structure-trained ML ablation is archived, not deployed:

- 114 mapped structural rows;
- 6 X->Cys observations;
- no improvement observed in that small subset;
- underpowered;
- useful as scientific iteration history.

## Remaining Scientific Limitations

- FireProtDB coverage is biased.
- Protein-family imbalance remains.
- Static PDB structures do not represent conformational ensembles.
- B-factors are not pure flexibility measurements.
- SASA/contact descriptors are simplified.
- Stability does not equal activity.
- Stability does not equal immobilization performance.
- Ranking weights require experimental calibration for application-specific deployment.

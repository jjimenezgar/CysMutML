# CysMutML v1.0 Current Project Status

Date: 2026-08-25

## Status

CysMutML v1.0 is release-prepared as a lightweight hybrid ML and structural-bioinformatics pipeline for cysteine mutation prioritization.

## Production Architecture

```text
FireProtDB physicochemical ML -> predicted_destabilization_ddg
Target PDB structural analysis -> SASA + B-factor-derived rigidity + optional penalties
Both components -> cys_suitability_score
```

## Deployed Model

- Artifact: `models/cysmutml_model.joblib`
- Model: Ridge regression
- Package version in metadata: `1.0.0`
- Dataset: FireProtDB mutation-level median-aggregated features
- Training rows: 352,005
- Unique proteins: 542
- Aggregated X->Cys rows: 16,236
- Feature set: 21 physicochemical features
- Target: `destabilization_ddg_kcal_mol`
- Convention: larger positive values mean greater destabilization

## Verified Metrics

Overall 3-fold GroupKFold by `protein_id`:

| Model | MAE | RMSE | R2 | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| Dummy mean | 0.800 | 1.049 | -0.002 | undefined | undefined |
| Ridge | 0.684 | 0.930 | 0.213 | 0.462 | 0.450 |
| HistGradientBoosting | 0.669 | 0.917 | 0.234 | 0.485 | 0.475 |

X->Cys subset:

| Model | MAE | RMSE | R2 | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| Dummy mean | 0.739 | 0.923 | -0.172 | undefined | undefined |
| Ridge | 0.587 | 0.803 | 0.115 | 0.348 | 0.345 |
| HistGradientBoosting | 0.579 | 0.795 | 0.131 | 0.364 | 0.362 |

Ridge remains deployed because HGB improves performance only modestly and Ridge is simpler, faster, and more interpretable.

## Real Case Study

- PDB: `1csp`
- Chain: `A`
- Candidates: 67
- Top candidate: `F38C`
- Outputs: `examples/real_case/`

## Release Artifacts

- `docs/FEATURE_SCHEMA.md`
- `docs/RANKING_FORMULA.md`
- `docs/RANKING_SENSITIVITY.md`
- `docs/ML_VS_HEURISTIC.md`
- `docs/SCIENTIFIC_AUDIT.md`
- `docs/RELEASE_CHECKLIST.md`
- `docs/INTERVIEW_GUIDE.md`
- `results/model_interpretability/ridge_coefficients.csv`
- `docs/figures/cysmutml_workflow.png`

## Exploratory Work

The structure-trained ML ablation is preserved as exploratory, not production:

- 114 mapped rows
- 6 X->Cys rows
- underpowered
- not promoted

## Remaining Manual Actions

- Initialize/push the GitHub repository.
- Create the GitHub release.
- Optionally add screenshots to the README.

## Final Local Checks

- `pytest -q`: 15 passed.
- `ruff check .`: all checks passed.
- `pip wheel . --no-deps`: built `cysmutml-1.0.0-py3-none-any.whl`.

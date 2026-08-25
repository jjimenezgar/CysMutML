# CysMutML v1.0 Project Status

Date: 2026-08-25

## V1.0 COMPLETE AND VERIFIED

- Hybrid architecture finalized: physicochemical ML plus target-PDB structural heuristic.
- FireProtDB v2.0 preprocessing implemented.
- Duplicate audit and mutation-level median aggregation implemented.
- Primary training table: `data/processed/fireprotdb_aggregated_features.csv`.
- Deployed model: Ridge regression in `models/cysmutml_model.joblib`.
- Target: `destabilization_ddg_kcal_mol`, larger positive values mean greater destabilization.
- Feature schema documented in `docs/FEATURE_SCHEMA.md`.
- Ranking formula documented in `docs/RANKING_FORMULA.md`.
- Model comparison completed with Dummy, Ridge, and HistGradientBoosting.
- Cys-specific metrics reported.
- Real case study completed for PDB `1csp` chain A.
- Ranking sensitivity analysis completed.
- Ridge coefficient interpretability output generated.
- PyMOL script and score-encoded PDB generation implemented.
- Release documentation created.
- Final tests pass: 15 passed.
- Final lint passes: Ruff all checks passed.
- Wheel build verified: `cysmutml-1.0.0-py3-none-any.whl`.

## CURRENT METRICS

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

## REAL CASE STUDY

- PDB: `1csp`
- Chain: `A`
- Generated X->Cys candidates: `67`
- Top candidate: `F38C`
- Outputs: `examples/real_case/`

## EXPLORATORY WORK

The structure-trained ML ablation remains archived:

- `results/structural_ablation/`
- `reports/structural_ablation_conclusion.md`
- `docs/STRUCTURAL_ERROR_ANALYSIS.md`

It used 114 mapped rows and 6 X->Cys observations. It is underpowered and not part of the v1.0 production architecture.

## KNOWN LIMITATIONS

- The model predicts mutation-associated destabilization, not immobilization success.
- Stability does not guarantee activity.
- Relative SASA does not guarantee cysteine reactivity.
- B-factor-derived rigidity is a lightweight crystallographic proxy, not a direct dynamics measurement.
- Ranking weights are heuristic and not experimentally calibrated.
- FireProtDB contains heterogeneous measurements across proteins, assays, temperatures, and pH.
- The deployed Ridge model is intentionally simple and not state-of-the-art.

## FUTURE V1.1 / V2 IDEAS

- Add an external S669 benchmark with overlap audit.
- Experimentally calibrate ranking weights if real Cys-engineering outcomes become available.
- Compare ESM sequence embeddings as an optional ML feature family.
- Add ANM/ProDy rigidity as an optional heuristic component.
- Improve uncertainty estimation.

## RELEASE NOTES

- License: MIT.
- Package version: `1.0.0`.
- Repository is not currently initialized as a git repository in this workspace.
- Large raw/processed structure/data files are ignored by `.gitignore`; selected lightweight v1.0 artifacts are allowed for portfolio review.

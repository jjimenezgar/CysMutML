# CysMutML Project Status

## Portfolio MVP closure — 9 September 2026

Ridge remains the demonstration model; ESM is an executed offline ablation. README, model card and Streamlit explicitly describe exploratory predictions and unresolved dataset identities. Source reconstruction and further retraining are deferred. This closes the portfolio implementation scope, not scientific validation of enzyme mutation outcomes.



Release baseline date: 2026-08-25

Portfolio hardening update: 2026-09-02

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
- Current automated test suite: 24 tests across Python 3.10 and 3.12 CI.
- Final lint passes: Ruff all checks passed.
- Wheel build verified: `cysmutml-1.0.0-py3-none-any.whl`.

## PORTFOLIO HARDENING

- Added GitHub Actions CI for Python 3.10 and 3.12.
- Made tests self-contained for clean repository checkouts.
- Added an executable end-to-end portfolio notebook.
- Added a model card documenting intended use, evaluation, and limitations.
- Added a concise reviewer path to the README.
- Added MMseqs2 sequence-cluster generation and homology-aware CV comparison.
- Added explicit leakage guards so cluster identifiers cannot enter model features.
- Added matched-subset coverage auditing for proteins without canonical sequences.
- Added deterministic, cluster-complete sampling for a lightweight 150-protein MVP.
- Added Dummy, Ridge, Random Forest, and HistGradientBoosting comparison with runtime measurements.
- Added X→Cys metrics, Ridge coefficients, held-out tree-model permutation importance, and portfolio figures.
- Added a four-tab Streamlit portfolio app with real PDB inference and downloadable outputs.

## CURRENT MODEL METRICS

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


## HOMOLOGY-AWARE MVP (EXECUTED)

- Run: MMseqs2 at 30% identity / 80% coverage, deterministic seed 42.
- Coverage: 543 source protein names; 171 mapped into 157 clusters; 372 excluded without usable sequence.
- Matched subset: 150 proteins and 5,634 mutation rows; four models; three folds.
- Overall MAE: protein grouped — Dummy 1.493, Ridge 1.508, RF 1.538, HGB 1.529.
- Overall MAE: homology clustered — Dummy 1.499, Ridge 1.523, RF 1.544, HGB 1.534.
- X→Cys Ridge MAE: 1.535 protein grouped versus 1.630 homology clustered.
- Full fold metrics, timings, permutation importance, and audit: results/homology_validation/.

## REAL CASE STUDY

- PDB: `1csp`
- Chain: `A`
- Generated X->Cys candidates: `67`
- Top candidate: `F38C`
- Outputs: `examples/real_case/`

## GODOY 2011 RETROSPECTIVE VALIDATION

Status: executed and documented.

- Heuristic was upgraded and frozen before extracting Godoy outcome tables.
- Frozen config: `validation/godoy2011/prevalidation_config.yaml`.
- Frozen SHA256 verified: `361b2f6cb47c1e07fc28005bcca8a4fdff4a3987d2f69b4428642120449641a6`.
- FireProtDB ML model was not retrained.
- Structures used: PGA `1K5Q`; BTL2 `2W22`.
- BTL2 cysteine-free experimental background accounted for with `2w22_c64s_c295s_context.pdb`.
- Mutants evaluated: 13 total, 6 PGA and 7 BTL2.
- Report: `validation/godoy2011/VALIDATION_REPORT.md`.

Main retrospective validation results:

| Comparison | n | Pearson | Spearman |
|---|---:|---:|---:|
| ML stability vs relative soluble activity | 13 | 0.122 | 0.087 |
| Calculated vs reported accessibility | 12 | 0.718 | 0.650 |
| Calculated vs reported exposed Lys count | 13 | 0.508 | 0.472 |
| Combined rigidification vs thermal stabilization | 13 | 0.472 | 0.522 |
| Combined rigidification vs solvent stabilization | 13 | 0.429 | 0.468 |

Top-k recovery using `final_engineering_score`:

| Enzyme | Top 5 | Top 10 | Top 20% | Top 30% |
|---|---:|---:|---:|---:|
| PGA | 0/6 | 1/6 | 4/6 | 4/6 |
| BTL2 | 0/7 | 0/7 | 3/7 | 4/7 |

Interpretation: the heuristic provides useful component-level diagnostics and partial enrichment, but this is not sufficient to claim calibrated prediction of immobilization success.

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
- Add ANM/ProDy rigidity as an optional heuristic component.
- Improve uncertainty estimation.
- Execute a 30%/40%/50% homology-threshold sensitivity study after the MVP, if compute budget justifies it.

## RELEASE NOTES

- License: MIT.
- Package version: `1.0.0`.
- Repository is initialized and has been pushed to GitHub.
- Large raw/processed structure/data files are ignored by `.gitignore`; selected lightweight v1.0 artifacts are allowed for portfolio review.


## Benchmark synchronization — 6 September 2026

The deployed Ridge model is unchanged. Full physicochemical metrics now match successful Actions run 34040897260 (training source commit 9b8613b868a333a5f6eda35b47575b19bec7c58b). Training: 351,487 aggregated rows; 543 evaluation groups; 16,208 X-to-Cys rows. Missing-name audit: 16 measurements formed 10 mutations across six UniProt identities; no cross-protein collision was found in those groups. A pre-aggregation identity guard now rejects conflicting identifiers and unresolved anonymous duplicates. No new training was required. Historical homology and invalid BTL2 results must not be interpreted as validation of the corrected model.


## Enzyme ESM feasibility — 9 September 2026

Successful run: https://github.com/jjimenezgar/CysMutML/actions/runs/34327672338

Using the frozen audited dataset, the strict reviewed-UniProt plus EC filter and exact sequence/WT checks retained 2,043 mutations from 64 unique sequences, including 46 X-to-Cys mutations. This is annotation/mapping coverage, not proof that excluded records are non-enzymes. No model was trained or promoted.

ESM-2 8M CPU pilot (two threads): 96/228/809 residues took 0.022/0.039/0.178 seconds for tokenization and inference after model load. Peak process RSS: 467.2 MiB; model load: 1.06 seconds in that runner. These timings exclude dependency installation and are not Streamlit latency guarantees. Full details and UniProt annotation snapshot are retained in the workflow artifact. The small Cys subset limits robust downstream validation; enzyme-only specialization remains experimental.


## ESM context comparison — 9 September 2026

Completed [run 34328356730](https://github.com/jjimenezgar/CysMutML/actions/runs/34328356730): 5,506 verified mutations, 170 sequences, 142 groups and three matched homology-grouped folds. Compared mean baseline, physicochemical Ridge, ESM context-only Ridge and their combination. Combined MAE: 1.466 overall and 1.075 for enzyme X→Cys (46 rows, 11 groups); respective physicochemical MAE: 1.487 and 1.501. Global combined R²: −0.206. Results are exploratory and not directly comparable to the full production benchmark. No production model or Streamlit inference change. Methods, limitations and artifacts: [ESM context comparison](https://github.com/jjimenezgar/CysMutML/blob/experiment/enzyme-esm-feasibility/docs/ESM_CONTEXT_COMPARISON.md).


## Sequence recovery audit — 9 September 2026

Completed frozen-data audit: 345,638 rows lack sequence; only four have UniProt. Found 68 conflicting WT positions across 29 protein names. Original Tsuboyama Fig. 3 source labels match 297 names (271,718 rows; 12,344 X→Cys), but these are recovery candidates, not verified full constructs. No new training-ready sequences or model changes. Reconstruct source-level sequence/assay provenance and reconcile target values before retraining. See [sequence recovery audit](docs/SEQUENCE_RECOVERY_AUDIT.md) and results/sequence_recovery_audit/.

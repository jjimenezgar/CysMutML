# Execution Log

Date/time: 2026-08-24 Europe/Madrid

Command: `python3 -m venv .venv`
Purpose: Create isolated environment.
Result: Succeeded.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/pip install -e '.[dev]'`
Purpose: Install project dependencies.
Result: Initial sandboxed run failed due DNS; escalated network run succeeded.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml prepare-data --download-fireprotdb --raw data/raw/fireprotdb.csv --output data/processed/fireprotdb_mutations.csv`
Purpose: Download and preprocess FireProtDB.
Result: Download succeeded; first parser attempt failed on live column names. Parser was fixed.
Status: PARTIALLY VERIFIED

Command: `.venv/bin/cysmutml prepare-data --raw data/raw/fireprotdb.csv --output data/processed/fireprotdb_mutations.csv`
Purpose: Normalize downloaded FireProtDB CSV.
Result: 613,208 raw rows; 555,932 valid single substitutions with DDG; 57,276 rejected.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml build-features --input data/processed/fireprotdb_mutations.csv --output data/processed/fireprotdb_features.csv`
Purpose: Build physicochemical features.
Result: 555,932 feature rows written.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml evaluate --features data/processed/fireprotdb_features.csv --results-dir results/fireprotdb`
Purpose: Full model-suite GroupKFold evaluation.
Result: Interrupted after several minutes due impractical interactive runtime.
Status: NOT EXECUTED TO COMPLETION

Command: `.venv/bin/cysmutml evaluate --features data/processed/fireprotdb_features.csv --results-dir results/fireprotdb_ridge --models dummy_mean,ridge`
Purpose: Reduced full-data GroupKFold evaluation.
Result: Interrupted after several minutes; no full scientific metrics claimed.
Status: NOT EXECUTED TO COMPLETION

Command: `.venv/bin/cysmutml evaluate-fast --features data/processed/fireprotdb_features.csv --results-dir results/fireprotdb_fast_baselines`
Purpose: Fast full FireProtDB grouped CV for Dummy and Ridge without out-of-fold prediction files.
Result: Completed. Overall mean MAE: Dummy 0.949, Ridge 0.852 kcal/mol. Cys-only mean MAE: Dummy 0.906, Ridge 0.731 kcal/mol.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml train --features data/processed/fireprotdb_features.csv --model-output models/cysmutml_model.joblib --metadata-output models/model_metadata.json --model-name ridge`
Purpose: Train first real physicochemical Ridge model on all processed FireProtDB rows.
Result: Completed; model metadata records `feature_configuration: physicochemical`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb examples/tiny_protein.pdb --chain A --model models/cysmutml_model.joblib --output examples/output_real`
Purpose: Verify X->Cys prediction using the real Ridge model.
Result: Completed; wrote `examples/output_real/mutation_predictions.csv`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml rank --predictions examples/output_real/mutation_predictions.csv --output examples/output_real/residue_ranking.csv`
Purpose: Generate separate engineering ranking from real-model predictions.
Result: Completed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml evaluate-fast --features data/processed/fireprotdb_features.csv --results-dir results/fireprotdb_fast_baselines --save-oof`
Purpose: Generate efficient out-of-fold predictions for Dummy and Ridge baselines.
Result: Completed; wrote `results/fireprotdb_fast_baselines/fast_baseline_out_of_fold_predictions.csv`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml error-analysis --predictions results/fireprotdb_fast_baselines/fast_baseline_out_of_fold_predictions.csv --output-dir results/fireprotdb_fast_baselines --model ridge_fast`
Purpose: Generate Ridge residual reports and figures from real out-of-fold predictions.
Result: Completed; wrote largest residuals, mutant-AA metrics, Cys OOF predictions, and three figures.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/pytest -q`
Purpose: Run tests.
Result: 8 passed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/ruff check .`
Purpose: Lint code.
Result: All checks passed.
Status: EXECUTED AND VERIFIED

Command: Synthetic prepare/build/evaluate/train/predict/rank/visualize commands.
Purpose: Verify software workflow without relying on scientific data.
Result: Completed; outputs in `examples/output/` and `results/synthetic/`.
Status: EXECUTED AND VERIFIED AS TOY SOFTWARE TEST

Command: `.venv/bin/cysmutml ablation --features data/processed/synthetic_features.csv --results-dir results/synthetic_ablation`
Purpose: Verify ablation machinery.
Result: Completed on toy data.
Status: EXECUTED AND VERIFIED AS TOY SOFTWARE TEST

Command: Python matplotlib one-liner to plot `results/synthetic_ablation/ablation_results.csv`.
Purpose: Verify figure output path.
Result: Wrote `results/figures/toy_ablation_mae.png`; matplotlib used a temporary cache because the home config path was not writable.
Status: EXECUTED AND VERIFIED AS TOY SOFTWARE TEST

Command: `.venv/bin/cysmutml audit-data --processed data/processed/fireprotdb_mutations.csv --raw data/raw/fireprotdb.csv`
Purpose: Audit duplicates, aggregate mutations, audit DDG distribution/outliers, and identify structural candidates.
Result: Completed; wrote duplicate, DDG, aggregation, and structure candidate reports.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml download-structures --candidates reports/structure_candidate_records.csv --structures-dir data/structures --sequences-dir data/sequences`
Purpose: Download/cache PDB structures and UniProt sequences.
Result: Initial sandboxed run failed; network-approved run downloaded 223 valid PDBs and 168/170 UniProt FASTA files.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml build-structural-features --candidates reports/structure_candidate_records_balanced_small.csv --min-identity 0.9`
Purpose: Build high-confidence structural mapping/features on a balanced-small execution subset.
Result: 298 candidates attempted; 114 mapped and feature-complete rows.
Status: PARTIALLY EXECUTED

Command: `.venv/bin/cysmutml structural-ablation --structural-features data/processed/fireprotdb_structural_features.csv --folds data/processed/structural_cv_folds.csv --results-dir results/structural_ablation`
Purpose: Compare physchem-only versus physchem-plus-structure on identical mapped rows and folds.
Result: Completed on 114 rows; structural features did not improve MAE on the executed subset.
Status: PARTIALLY EXECUTED

## 2026-08-25 Hybrid Simplification Milestone

Command: `.venv/bin/cysmutml build-features --input data/processed/fireprotdb_mutations_aggregated.csv --output data/processed/fireprotdb_aggregated_features.csv`
Purpose: Build physicochemical features for the median-aggregated mutation-level dataset.
Result: Completed; wrote 352,005 feature rows.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml compare-physchem --features data/processed/fireprotdb_aggregated_features.csv --results-dir results/physchem_model_comparison --models dummy_mean,ridge,hist_gradient_boosting`
Purpose: Compare physicochemical-only models on the aggregated dataset.
Result: Initial run produced implausibly strong metrics because target-derived aggregate columns were still eligible features. The issue was identified as leakage and fixed.
Status: EXECUTED BUT INVALIDATED BY LEAKAGE AUDIT

Command: `.venv/bin/cysmutml compare-physchem --features data/processed/fireprotdb_aggregated_features.csv --results-dir results/physchem_model_comparison --models dummy_mean,ridge,hist_gradient_boosting`
Purpose: Rerun the physicochemical-only model comparison after excluding target-derived aggregate columns.
Result: Completed. Mean MAE: Dummy 0.800, Ridge 0.684, HistGradientBoosting 0.669. Cys-only mean MAE: Dummy 0.739, Ridge 0.587, HistGradientBoosting 0.579.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml train --features data/processed/fireprotdb_aggregated_features.csv --model-output models/cysmutml_model.joblib --metadata-output models/model_metadata.json --model-name ridge`
Purpose: Retrain deployed model on the aggregated physicochemical dataset.
Result: Completed; model remains Ridge with physicochemical features only.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb examples/real_case/1csp.pdb --chain A --target CYS --model models/cysmutml_model.joblib --output examples/real_case`
Purpose: Execute the real hybrid case study.
Result: Completed; wrote ML predictions, engineering ranking, score-encoded PDB, and PyMOL script for 67 X->Cys candidates.
Status: EXECUTED AND VERIFIED

Command: Python ranking sensitivity script.
Purpose: Recalculate the 1CSP ranking under balanced, stability-heavy, and accessibility-heavy heuristic weights.
Result: Completed; wrote `results/ranking_sensitivity/`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/pytest -q`
Purpose: Verify the updated hybrid pipeline and tests.
Result: 12 passed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/ruff check .`
Purpose: Lint the updated codebase.
Result: All checks passed.
Status: EXECUTED AND VERIFIED

## 2026-08-25 v1.0 Release Audit

Command: `.venv/bin/cysmutml train --features data/processed/fireprotdb_aggregated_features.csv --model-output models/cysmutml_model.joblib --metadata-output models/model_metadata.json --model-name ridge`
Purpose: Rebuild the deployed model after setting package version to 1.0.0.
Result: Completed; metadata reports `cysmutml: 1.0.0`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb examples/real_case/1csp.pdb --chain A --target CYS --model models/cysmutml_model.joblib --output examples/real_case`
Purpose: Regenerate the v1.0 real case study outputs.
Result: Completed; generated 67 X->Cys candidates and top-ranked `F38C`.
Status: EXECUTED AND VERIFIED

Command: Python script for Ridge coefficients and workflow figure.
Purpose: Generate v1.0 interpretability and documentation figure artifacts.
Result: Wrote `results/model_interpretability/ridge_coefficients.csv`, `results/model_interpretability/ridge_coefficients_top15.png`, and `docs/figures/cysmutml_workflow.png`.
Status: EXECUTED AND VERIFIED

Command: Python script for ranking sensitivity.
Purpose: Recalculate balanced, stability-heavy, and accessibility-heavy rankings for 1CSP.
Result: Top-5 overlap versus balanced was 4 for stability-heavy and 2 for accessibility-heavy; Spearman correlations were 0.947 and 0.692.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml --help`
Purpose: Verify CLI entry point.
Result: Help displayed successfully.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --help`
Purpose: Verify primary workflow help.
Result: Help displayed successfully.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/pytest -q`
Purpose: Final v1.0 test run.
Result: 15 passed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/ruff check .`
Purpose: Final v1.0 lint run.
Result: All checks passed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/python -m build`
Purpose: Check package build.
Result: Failed because the `build` module is not installed in the environment.
Status: NOT EXECUTED TO COMPLETION

Command: `.venv/bin/pip wheel . --no-deps --wheel-dir /tmp/cysmutml_wheel`
Purpose: Verify package wheel build from `pyproject.toml`.
Result: Initial sandboxed run failed due blocked network while resolving build dependencies; approved network run succeeded and built `cysmutml-1.0.0-py3-none-any.whl`.
Status: EXECUTED AND VERIFIED

## 2026-08-25 Godoy 2011 Heuristic Validation Milestone

Command: `.venv/bin/pytest -q`
Purpose: Verify upgraded heuristic implementation before outcome-table extraction.
Result: 16 passed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/ruff check .`
Purpose: Lint upgraded heuristic implementation before outcome-table extraction.
Result: All checks passed.
Status: EXECUTED AND VERIFIED

Command: `git commit -m "Freeze Godoy prevalidation heuristic"`
Purpose: Freeze scoring formulas and code before extracting Godoy outcome tables.
Result: Commit `1715ae3` created.
Status: EXECUTED AND VERIFIED

Command: `git commit -m "Record Godoy prevalidation freeze"`
Purpose: Record frozen config hash and prevalidation metadata.
Result: Commit `3c8a972` created.
Status: EXECUTED AND VERIFIED

Command: `pdftotext -layout ...`
Purpose: Extract text from supplied Godoy main paper and supporting-information PDFs.
Result: Temporary local text extraction succeeded and was used to transcribe validation tables. The extracted article text was not retained in the repository.
Status: EXECUTED AND VERIFIED

Command: `curl -L ...`
Purpose: Download PDB structures specified by the Godoy supporting information.
Result: Downloaded `1k5q.pdb` and `2w22.pdb`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb validation/godoy2011/structures/1k5q.pdb --chain A --model models/cysmutml_model.joblib --output validation/godoy2011/predictions/PGA_alpha`
Purpose: Generate PGA alpha-chain X->Cys predictions and rankings.
Result: Completed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb validation/godoy2011/structures/1k5q.pdb --chain B --model models/cysmutml_model.joblib --output validation/godoy2011/predictions/PGA_beta`
Purpose: Generate PGA beta-chain X->Cys predictions and rankings.
Result: Completed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb validation/godoy2011/structures/2w22.pdb --chain A --model models/cysmutml_model.joblib --output validation/godoy2011/predictions/BTL2_wt --monocysteine-design`
Purpose: Generate BTL2 wild-type-context X->Cys predictions and rankings.
Result: Completed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/cysmutml predict --pdb validation/godoy2011/structures/2w22_c64s_c295s_context.pdb --chain A --model models/cysmutml_model.joblib --output validation/godoy2011/predictions/BTL2_cysfree_context --monocysteine-design`
Purpose: Generate BTL2 predictions in the cysteine-free experimental background.
Result: Completed after correcting experimental C64/C295 labels to PDB residues A:65/A:296.
Status: EXECUTED AND VERIFIED

Command: Python validation-matrix script.
Purpose: Compare CysMutML scores against Godoy accessibility, Lys counts, soluble activity, immobilization recovery, and stabilization factors.
Result: Wrote validation CSV files and figures under `validation/godoy2011/`.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/pytest -q`
Purpose: Final verification after adding Godoy validation artifacts and documentation.
Result: 16 passed.
Status: EXECUTED AND VERIFIED

Command: `.venv/bin/ruff check .`
Purpose: Final lint verification after adding Godoy validation artifacts and documentation.
Result: All checks passed.
Status: EXECUTED AND VERIFIED

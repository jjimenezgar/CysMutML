# AGENTS.md

Future agents must preserve the scientific contract of CysMutML v1.0.

## Objective
Use physicochemical ML to predict mutation-associated destabilization and separately use target-PDB structural heuristics to rank X->Cys engineering candidates.

## Never Silently Change
- Internal target: `destabilization_ddg_kcal_mol`.
- Sign convention: larger positive values mean greater destabilization.
- FireProtDB DDG is currently used directly because its docs define positive DDG as destabilizing.
- Do not use random mutation-level splits as primary validation.
- Do not report homology-clustered metrics unless the cluster mapping, thresholds, coverage, and fold strategy are recorded.
- Never allow `sequence_cluster` or representative IDs into the model feature matrix.
- Do not merge the ML stability model with the engineering ranking heuristic.
- Do not add structural descriptors to the deployed ML model unless a future task explicitly changes the project scope and validates it leakage-safely.
- Do not claim toy outputs as scientific results.
- Do not call `stability_component` or `cys_suitability_score` a calibrated probability.

## Structural Features
SASA uses Biopython Shrake-Rupley. Relative SASA uses Tien 2013 maximum ASA. Contacts are CA neighbors at configured radii and heavy-atom proximity counts.

In the current production architecture, structure is a ranking heuristic on the target PDB, not a training feature group.

## Mapping
Dataset numbering must be aligned to PDB chain sequence and WT identity must match. Failed mappings are excluded or explicitly flagged.

## Commands
```bash
.venv/bin/pytest -q
.venv/bin/ruff check .
cysmutml prepare-data --download-fireprotdb --raw data/raw/fireprotdb.csv --output data/processed/fireprotdb_mutations.csv
cysmutml build-features --input data/processed/fireprotdb_mutations.csv --output data/processed/fireprotdb_features.csv
cysmutml build-features --input data/processed/fireprotdb_mutations_aggregated.csv --output data/processed/fireprotdb_aggregated_features.csv
cysmutml compare-physchem --features data/processed/fireprotdb_aggregated_features.csv --results-dir results/physchem_model_comparison
cysmutml predict --pdb data/structures/1csp.pdb --chain A --target CYS --output examples/real_case
.venv/bin/pytest -q
.venv/bin/ruff check .
```

## Documentation
Update `docs/EXECUTION_LOG.md`, `PROJECT_STATUS.md`, `docs/SCIENTIFIC_AUDIT.md`, and `reports/run_summary.json` whenever execution status changes.

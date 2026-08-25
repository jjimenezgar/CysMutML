# CysMutML v1.0 Release Checklist

Date: 2026-08-25

- [x] Package version set to `1.0.0`.
- [x] Deployed model artifact exists: `models/cysmutml_model.joblib`.
- [x] Model metadata exists: `models/model_metadata.json`.
- [x] Deployed model verified as Ridge regression.
- [x] Training dataset verified as FireProtDB mutation-level median-aggregated features.
- [x] Aggregated training row count verified: 352,005.
- [x] Target convention documented: larger positive `destabilization_ddg_kcal_mol` means greater destabilization.
- [x] GroupKFold metrics verified from actual result files.
- [x] Cys-specific metrics verified from actual result files.
- [x] Feature schema documented.
- [x] Ranking formula documented.
- [x] ML vs heuristic distinction documented.
- [x] Real case study runs for `1csp` chain A.
- [x] Real case input PDB included at `examples/real_case/1csp.pdb`.
- [x] Real case output files generated.
- [x] Ranking sensitivity analysis generated.
- [x] Ridge coefficient interpretability output generated.
- [x] README updated for v1.0.
- [x] Scientific limitations visible.
- [x] License present: MIT.
- [x] Citation file present.
- [x] `pytest -q` passes: 15 passed.
- [x] `ruff check .` passes.
- [x] Wheel build passes: `cysmutml-1.0.0-py3-none-any.whl`.
- [x] CLI help works.
- [x] No fake metrics knowingly reported.
- [x] `.gitignore` excludes large raw/processed data and caches while allowing selected lightweight v1.0 artifacts.
- [ ] Remote GitHub repository initialized.
- [ ] GitHub release created.
- [ ] Optional README screenshots added.

Final local release status:

```text
Local v1.0 release audit passed.
```

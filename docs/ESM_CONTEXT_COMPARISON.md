# ESM sequence-context experiment

Status: completed on 9 September 2026; experimental, not deployed.

## Question and decision

Does a frozen protein language model improve mutation-associated destabilization prediction on unseen protein groups?

Adding ESM context reduced error relative to the physicochemical Ridge on this cohort, but the global improvement was small and global R² remained negative. **Keep the existing production model.** This experiment does not establish reliable position-specific mutation prediction or immobilization success.

## Data and validation

The frozen audited FireProtDB feature table contained 351,487 aggregated mutations. We retained 5,506 rows from 170 exact sequences and 142 evaluation groups. Exclusions were 345,638 rows without sequence, 266 with unsupported length/alphabet, 44 without an exact verified UniProt reference and 33 with position or wild-type mismatch. Sequences were limited to 1,000 standard amino acids.

The enzyme subset requires reviewed UniProt status and an EC annotation: 1,998 rows in 47 groups. Missing annotation does not mean a protein is not an enzyme. There are 139 X→Cys rows in 25 groups overall, including only 46 enzyme X→Cys rows in 11 groups. The enzyme pilot retained 2,043 rows before the length/alphabet restriction.

MMseqs2 clustered sequences at 30% identity and 80% bidirectional coverage (`--cov-mode 0`). Clusters sharing a protein identifier or UniProt accession were joined. Three GroupKFold splits kept each group intact, with explicit checks for sequence/group overlap. All models used identical rows and folds. Clusters reduce relatedness leakage but cannot exclude every remote homolog.

## Models

- **Mean baseline:** predicts the training-fold mean.
- **Physicochemical Ridge:** existing amino-acid descriptors, substitution deltas and BLOSUM62 features.
- **ESM context control:** Ridge using the wild-type residue embedding only. It has no mutant identity, so different substitutions at the same position receive the same prediction. This is a context-only ablation, not a general mutation-effect predictor.
- **Physicochemical + ESM:** concatenates substitution descriptors and reduced residue context.

ESM-2 8M (`facebook/esm2_t6_8M_UR50D`) was frozen at revision `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`. The final-layer wild-type residue embedding has 320 dimensions. Embeddings were standardized, reduced to 16 PCA components, then standardized again, all fitted exclusively on each training fold. Ridge alpha was fixed at 1.0. No supervised ESM fine-tuning or hyperparameter search was performed.

These embeddings supply sequence context; they are not explicit 3D contacts, solvent exposure, or a physical simulation of cysteine interactions. Training occurred in GitHub Actions, outside Streamlit.

## Results

Pooled out-of-fold MAE in kcal/mol; lower is better:

| Model | All mutations | X→Cys | Enzymes | Enzymes X→Cys |
|---|---:|---:|---:|---:|
| Mean baseline | 1.474 | 1.691 | 1.482 | 1.142 |
| Physicochemical Ridge | 1.487 | 1.811 | 1.469 | 1.501 |
| ESM context control | 1.447 | 1.626 | 1.395 | 1.011 |
| Physicochemical + ESM | 1.466 | 1.593 | 1.400 | 1.075 |

The combined model's global R² is −0.206; enzyme X→Cys R² is −0.102. Its enzyme X→Cys MAE improves substantially over physicochemical Ridge, but only about 6% over the mean baseline. With 46 observations from 11 groups, this is exploratory evidence, not a validated improvement. The context-only control also has lower MAE than the combined model, so this run does not demonstrate added value from substitution descriptors.

MAE is average absolute error. RMSE emphasizes large errors. R² below zero means worse squared-error performance than the evaluated subset's constant observed mean. Pearson measures linear association; Spearman measures ranking association. Neither correlation alone demonstrates accurate ΔΔG prediction. Pooled baseline correlations can be nonzero because each fold predicts a different training mean.

See [rounded metrics](../results/esm_context_comparison/metrics_rounded.csv) for all metrics. These are pooled held-out metrics, not averages of fold metrics. **Do not compare them directly with the production benchmark's approximately 0.667 MAE:** the cohort, grouping and aggregation of metrics differ.

## Reproducibility and provenance

- Script: [compare_esm_context.py](../scripts/compare_esm_context.py).
- Workflow: [esm-context-comparison.yml](../.github/workflows/esm-context-comparison.yml).
- Successful [run 34328356730](https://github.com/jjimenezgar/CysMutML/actions/runs/34328356730), source commit `d41bbd38b704fb4735a5414d59a7e9142d2529a5`.
- [Full artifact](https://github.com/jjimenezgar/CysMutML/actions/runs/34328356730/artifacts/10094749100): full-precision metrics, out-of-fold predictions, row manifest, group mapping, fold counts, sequence FASTA, cached embeddings, fold estimators and run metadata. Actions artifacts have limited retention; the rounded summary is versioned in Git.
- Input features: audited retraining artifact from run 34040897260. SHA256: `fa6144f729df0ebe34305306d1e9f9ea20a8692bca46542b48fbd9e2ef634d60`.
- Frozen UniProt annotations: pilot artifact from run 34327672338. SHA256: `bc5a11d4de2baff21b2984e25598eacb4758c2500f8d73fd05d2c05dd87a2894`.
- Runtime: Python 3.12.14, scikit-learn 1.9.0, torch 2.14.0+cpu, transformers 5.16.1, MMseqs2 15-6f452+ds-2.
- Embedding loop took 4.54 seconds on that CPU runner using two threads, excluding model loading, dependency installation, clustering and regression.

To rerun with the frozen input files and workflow dependencies:

```bash
python scripts/compare_esm_context.py \
  --features frozen/data/processed/fireprotdb_aggregated_features.csv \
  --annotations annotation_snapshot/uniprot_annotations.csv \
  --output experiment_output
```

Use a fresh output directory: the embedding cache filenames identify the sequence, not the model revision.

## Limits and next decision

Sequence verification cannot resolve every experimental construct or numbering ambiguity. Assay heterogeneity remains, the sequence-complete cohort may be biased, and ESM pretraining overlap is not ruled out. This is a single exploratory comparison without an independent model-selection test or confidence intervals.

Before promotion, obtain broader verified sequence coverage and repeat model selection within training groups, reserving untouched groups for final assessment. Evaluate X→Cys across proteins as the primary downstream task. The existing structural heuristic and production model are unchanged.

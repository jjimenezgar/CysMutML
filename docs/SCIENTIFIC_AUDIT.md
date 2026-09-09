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
- B-factor-derived flexibility proxy;
- local exposed Lys count;
- optional protected-residue distance;
- existing-cysteine proximity.

The final ranking is separated into `cys_site_suitability`, `rigidification_potential`, and `final_engineering_score`. All three are reconstructable from output columns.

Cys-site suitability:

```text
cys_site_suitability =
  0.60 * stability_component
+ 0.35 * accessibility_component
- 0.10 * existing_cys_penalty
- 0.15 * protected_site_penalty
```

Rigidification potential:

```text
rigidification_potential =
  0.35 * flexibility_component
+ 0.40 * lysine_environment_component
+ 0.25 * accessibility_component
- 0.05 * existing_cys_penalty
- 0.10 * protected_site_penalty
```

Final engineering score:

```text
final_engineering_score =
  0.60 * cys_site_suitability
+ 0.40 * rigidification_potential
```

The weights are heuristic defaults and configurable in `configs/default.yaml`.

Protected residues:

- user supplied only;
- no active sites are inferred automatically;
- distances use C-alpha coordinates across chains.

Existing cysteines:

- proximity is reported as an engineering caution;
- proximity does not imply disulfide formation.

## Godoy 2011 Retrospective Validation

Status:

- Prevalidation heuristic freeze completed before inspecting outcome tables.
- Frozen config hash verified after validation.
- FireProtDB model was not retrained.
- Godoy validation report created at `validation/godoy2011/VALIDATION_REPORT.md`.

Audit findings:

- ML stability did not meaningfully correlate with relative soluble activity across the 13 Godoy mutants.
- Calculated accessibility had moderate overall association with reported accessibility, but BTL2 site-level mismatches were substantial.
- Calculated local exposed Lys counts had moderate rank association with reported Lys counts, but exact counts differed.
- Combined rigidification potential had moderate association with stabilization factors across both enzymes, but per-enzyme associations were weak or inconsistent.
- The final engineering score partially enriched experimental sites in the top 20-30%, but not in the top 5-10 candidates.

Conclusion:

The validation supports CysMutML as a transparent retrospective prioritization aid, not as a calibrated predictor of immobilization success.

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


## Homology-aware MVP (2026-09-02)

A reduced benchmark was executed after enriching the FireProtDB export with UniProt
reference sequences. MMseqs2 used 30% minimum sequence identity, 80% coverage, and
coverage mode 0. The source table contained 543 protein names; 171 mapped to 157
clusters, while 372 names without usable sequences were excluded. A deterministic
complete-cluster sample (seed 42) contained 150 proteins and 5,634 mutation rows.

| Split | Ridge MAE | HGB MAE | RF MAE | Dummy MAE |
|---|---:|---:|---:|---:|
| Protein grouped | 1.508 | 1.529 | 1.538 | 1.493 |
| Homology clustered | 1.523 | 1.534 | 1.544 | 1.499 |

For X→Cys rows, Ridge MAE was 1.535 with protein grouping and 1.630 with homology
grouping. The homology-aware estimate is therefore the more conservative headline
for transfer to less-related proteins, but the sample is deliberately small and
heterogeneous. It is evidence of validation design, not a state-of-the-art claim.

The benchmark also records fit/predict runtime, Ridge interpretability via
coefficients, and held-out permutation importance for the tree models. Cluster IDs,
representative IDs, and sequence provenance IDs are excluded from model features.


## Benchmark synchronization — 6 September 2026

The deployed Ridge model is unchanged. Full physicochemical metrics now match successful Actions run 34040897260 (training source commit 9b8613b868a333a5f6eda35b47575b19bec7c58b). Training: 351,487 aggregated rows; 543 evaluation groups; 16,208 X-to-Cys rows. Missing-name audit: 16 measurements formed 10 mutations across six UniProt identities; no cross-protein collision was found in those groups. A pre-aggregation identity guard now rejects conflicting identifiers and unresolved anonymous duplicates. No new training was required. Historical homology and invalid BTL2 results must not be interpreted as validation of the corrected model.


## Enzyme ESM feasibility — 9 September 2026

Successful run: https://github.com/jjimenezgar/CysMutML/actions/runs/34327672338

Using the frozen audited dataset, the strict reviewed-UniProt plus EC filter and exact sequence/WT checks retained 2,043 mutations from 64 unique sequences, including 46 X-to-Cys mutations. This is annotation/mapping coverage, not proof that excluded records are non-enzymes. No model was trained or promoted.

ESM-2 8M CPU pilot (two threads): 96/228/809 residues took 0.022/0.039/0.178 seconds for tokenization and inference after model load. Peak process RSS: 467.2 MiB; model load: 1.06 seconds in that runner. These timings exclude dependency installation and are not Streamlit latency guarantees. Full details and UniProt annotation snapshot are retained in the workflow artifact. The small Cys subset limits robust downstream validation; enzyme-only specialization remains experimental.


## ESM context comparison — 9 September 2026

Completed [run 34328356730](https://github.com/jjimenezgar/CysMutML/actions/runs/34328356730): 5,506 verified mutations, 170 sequences, 142 groups and three matched homology-grouped folds. Compared mean baseline, physicochemical Ridge, ESM context-only Ridge and their combination. Combined MAE: 1.466 overall and 1.075 for enzyme X→Cys (46 rows, 11 groups); respective physicochemical MAE: 1.487 and 1.501. Global combined R²: −0.206. Results are exploratory and not directly comparable to the full production benchmark. No production model or Streamlit inference change. Methods, limitations and artifacts: [ESM context comparison](https://github.com/jjimenezgar/CysMutML/blob/experiment/enzyme-esm-feasibility/docs/ESM_CONTEXT_COMPARISON.md).


## Sequence recovery audit — 9 September 2026

Completed frozen-data audit: 345,638 rows lack sequence; only four have UniProt. Found 68 conflicting WT positions across 29 protein names. Original Tsuboyama Fig. 3 source labels match 297 names (271,718 rows; 12,344 X→Cys), but these are recovery candidates, not verified full constructs. No new training-ready sequences or model changes. Reconstruct source-level sequence/assay provenance and reconcile target values before retraining. See [sequence recovery audit](SEQUENCE_RECOVERY_AUDIT.md) and results/sequence_recovery_audit/.

## Portfolio MVP closure — 9 September 2026

Ridge remains the demonstration model; ESM is an executed offline ablation. README, model card and Streamlit explicitly describe exploratory predictions and unresolved dataset identities. Source reconstruction and further retraining are deferred. This closes the portfolio implementation scope, not scientific validation of enzyme mutation outcomes.


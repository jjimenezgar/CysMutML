# Homology-aware validation

## Why this exists

Grouping mutations by `protein_id` prevents mutations from the same protein from
appearing in both training and test folds. It does not prevent closely homologous
proteins from being split across folds.

CysMutML v1.2 adds a stricter evaluation in which proteins are first clustered by
sequence identity and every cluster is kept inside one fold. This estimates
generalization to less-related protein families and makes residual homology
leakage visible.

The deployed v1.0 model is not changed by this experiment.

## Requirements

Install [MMseqs2](https://github.com/soedinglab/MMseqs2) and ensure the `mmseqs`
executable is on `PATH`.

The aggregated mutation table must contain:

```text
protein_id
canonical_sequence
```

When `prepare-data` is run with `--download-fireprotdb`, CysMutML also resolves each documented `SOURCE_SEQUENCE_ID` through the FireProtDB sequence endpoint and stores `canonical_sequence`. Existing local raw exports can be enriched with `--fetch-sequences`. The median-aggregated table preserves both the sequence and its FireProtDB identifier, so previously generated tables must be rebuilt once.

## Build sequence clusters

```bash
cysmutml audit-data \
  --processed data/processed/fireprotdb_mutations.csv \
  --raw data/raw/fireprotdb.csv

cysmutml build-features \
  --input data/processed/fireprotdb_mutations_aggregated.csv \
  --output data/processed/fireprotdb_aggregated_features.csv

cysmutml build-homology-clusters \
  --input data/processed/fireprotdb_mutations_aggregated.csv \
  --output data/processed/sequence_clusters.csv \
  --min-sequence-identity 0.30 \
  --coverage 0.80
```

The command uses `MMseqs2 easy-cluster` with bidirectional coverage mode 0 and
writes:

```text
data/processed/sequence_clusters.csv
data/processed/sequence_clusters.metadata.json
```

The mapping contains one cluster per protein with a valid canonical sequence. It fails if a protein has
conflicting canonical sequences or if MMseqs2 does not return every input
sequence.

## Compare grouping strategies

```bash
cysmutml compare-grouping-strategies \
  --features data/processed/fireprotdb_aggregated_features.csv \
  --clusters data/processed/sequence_clusters.csv \
  --results-dir results/homology_validation \
  --models dummy_mean,ridge,random_forest,hist_gradient_boosting \
  --target-proteins 150 \
  --random-seed 42
```

This selects whole sequence clusters with a deterministic seed until at least 150 proteins are included, then runs the same four physicochemical models with three folds under two CV strategies:

| Strategy | Test-set isolation |
|---|---|
| `protein_grouped` | No protein occurs in train and test |
| `homology_clustered` | No sequence cluster occurs in train and test |

Generated artifacts:

```text
results/homology_validation/protein_grouped/
results/homology_validation/homology_clustered/
results/homology_validation/split_comparison_fold_metrics.csv
results/homology_validation/split_comparison_summary.csv
results/homology_validation/split_comparison_cys_metrics.csv
results/homology_validation/split_comparison_cys_summary.csv
results/homology_validation/mvp_protein_manifest.csv
results/homology_validation/tree_permutation_importance.csv
results/homology_validation/cluster_audit.json
results/homology_validation/figures/*.png
```

Proteins without a cluster are excluded explicitly. Sampling never cuts a homology cluster: every protein from a selected cluster is included, so the realized total can slightly exceed 150. Both strategies then use exactly the same rows and precomputed fold manifest. `cluster_audit.json` records source, mapped, included, and excluded counts plus the sampling parameters. Grouping columns are excluded from model features. The comparison therefore changes only the split, not the dataset or model inputs.

Runtime is recorded separately for fitting and prediction. Ridge coefficients provide global linear interpretability. Random Forest and HistGradientBoosting receive held-out permutation importance on the first homology fold; this avoids presenting impurity importance as if it were model-agnostic.

## Interpretation

Report both results. Do not select whichever split gives the most favorable
metric.

Expected possibilities:

- Similar performance suggests the physicochemical signal transfers reasonably
  across the chosen sequence-identity boundary.
- Lower homology-clustered performance indicates that protein-grouped CV was
  benefiting from relationships between proteins in different folds.
- High fold variance indicates that the number or composition of clusters may be
  insufficient for a stable estimate.

A 30% identity threshold is intentionally strict but is not uniquely correct.
Sensitivity at 30%, 40%, and 50% identity is a useful follow-up. Thresholds must
be selected before inspecting downstream performance.

## Leakage guarantees

Automated tests verify that:

- each protein maps to exactly one cluster;
- incomplete and conflicting mappings fail loudly;
- a sequence cluster appears in exactly one test fold;
- cluster identifiers and representative identifiers never enter the ML feature matrix;
- reduced sampling is deterministic and preserves complete clusters;
- evaluation records fit and prediction runtime.

## Current status

The infrastructure, leakage guards, figure generation, and Streamlit integration are implemented. The benchmark workflow is manual-only because it downloads and processes FireProtDB; selected compact results are versioned after a successful run.

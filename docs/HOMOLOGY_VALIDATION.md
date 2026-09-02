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

CysMutML now preserves `canonical_sequence` during median aggregation. Existing
locally generated aggregated tables must be rebuilt once.

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
  --results-dir results/homology_validation
```

This runs identical physicochemical models under two CV strategies:

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
results/homology_validation/cluster_audit.json
```

Proteins without a cluster are excluded explicitly. Both strategies are then run
on exactly the same mapped rows, and `cluster_audit.json` records source, included,
and excluded counts. Grouping columns are excluded from model features. The
comparison therefore changes only the split, not the dataset or model inputs.

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
- cluster identifiers and representative identifiers never enter the ML feature matrix.

## Current status

The infrastructure and leakage guards are implemented and CI-tested. Numerical
homology-split metrics are intentionally not claimed in the repository until the
local FireProtDB feature table is regenerated and the MMseqs2 experiment is run.

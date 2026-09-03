# CysMutML

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

CysMutML is an interpretable ML and structural-ranking pipeline for designing cysteine substitutions as site-specific handles for enzyme immobilization.

It combines two separate signals:

1. a physicochemical ML model trained on FireProtDB mutation data;
2. structural information from a target PDB structure.

CysMutML does not model the immobilization chemistry itself; it prioritizes mutation sites that are predicted to preserve stability while remaining structurally accessible. The result is a design hypothesis, not a probability of immobilization success.

![CysMutML](docs/figures/cysmutml_github_cover.jpg)

[![Tests](https://img.shields.io/badge/tests-24%20passed-2ea44f?style=for-the-badge&logo=pytest&logoColor=white)](https://github.com/jjimenezgar/CysMutML/actions)
[![CI](https://img.shields.io/badge/CI-passing-2ea44f?style=for-the-badge&logo=githubactions&logoColor=white)](https://github.com/jjimenezgar/CysMutML/actions)
[![Launch Streamlit](https://img.shields.io/badge/Launch%20Streamlit-2ea44f?style=for-the-badge&logo=streamlit&logoColor=white)](https://cysmutml.streamlit.app)

## What the project demonstrates

- Data cleaning and aggregation for a heterogeneous protein dataset.
- Protein-aware and homology-aware cross-validation.
- Interpretable regression with a reproducible benchmark.
- Leakage checks and explicit feature contracts.
- End-to-end inference on a PDB structure.
- A lightweight Streamlit interface and downloadable analysis files.

## Current model

The deployed model is a Ridge regressor using amino-acid physicochemical descriptors, mutation deltas and BLOSUM62 features.

It was trained on 352,005 median-aggregated FireProtDB rows from 542 proteins. The target is:

```text
destabilization_ddg_kcal_mol
```

Larger positive values indicate greater predicted destabilisation. Structural descriptors are not used by this model.

CysMutML captures average physicochemical tendencies of substitutions, but it is not a position-specific mutational effect predictor.

## Homology-aware MVP

A reduced benchmark was run on 150 proteins and 5,634 mutation rows. Sequences were clustered with MMseqs2 at 30% identity and 80% coverage. The experiment used seed 42, complete sequence clusters and three folds.

| Split | Dummy | Ridge | Random Forest | HistGradientBoosting |
|---|---:|---:|---:|---:|
| Protein grouped | 1.493 | 1.508 | 1.538 | 1.529 |
| Homology clustered | 1.499 | 1.523 | 1.544 | 1.534 |

Mean MAE, lower is better. On the X→Cys subset, Ridge scored 1.535 with protein grouping and 1.630 with homology clustering.

**How to read the splits.** `Protein grouped` keeps every mutation from a protein in one fold, so test proteins are unseen during training. `Homology clustered` first groups similar sequences with MMseqs2 and keeps each whole cluster in one fold; it is a stricter test against residual sequence relatedness.

This is a small portfolio benchmark, not a state-of-the-art claim. The full fold metrics, timing measurements, sampling audit and permutation importance are in [docs/HOMOLOGY_VALIDATION.md](docs/HOMOLOGY_VALIDATION.md).

**What the results say.** The learned models improve on the mean baseline, but the explained variance is modest. Gradient boosting is slightly better than Ridge in this benchmark; Ridge remains deployed for its simplicity and interpretability. The increase in error under homology clustering is an explicit warning about generalisation to related proteins.

## Streamlit app

Run the app locally:

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[app]'
.venv/bin/streamlit run streamlit_app.py
```

The app supports:

- the bundled 1CSP example, PDB/mmCIF upload, an RCSB/PDB identifier or a UniProt accession;
- automatic structure retrieval from RCSB and UniProt/AlphaFold fallback;
- X→Cys predictions for a selected chain;
- a ranked table with readable component names and short explanations;
- an interactive 3D structure viewer with top-percentile and amino-acid filters;
- CSV, score-encoded PDB and PyMOL downloads;
- benchmark, methods and limitations summaries.

## Command-line example

```bash
cysmutml predict \
  --pdb examples/real_case/1csp.pdb \
  --chain A \
  --output examples/real_case
```

To reproduce the homology-aware benchmark, install MMseqs2 and run:

```bash
cysmutml build-homology-clusters \
  --input data/processed/fireprotdb_mutations_aggregated.csv \
  --output data/processed/sequence_clusters.csv \
  --min-sequence-identity 0.30 \
  --coverage 0.80

cysmutml compare-grouping-strategies \
  --features data/processed/fireprotdb_aggregated_features.csv \
  --clusters data/processed/sequence_clusters.csv \
  --models dummy_mean,ridge,random_forest,hist_gradient_boosting \
  --target-proteins 150 \
  --random-seed 42
```

## How the ranking is constructed

The ML model supplies a stability term for each possible X→Cys substitution. The structure supplies five simple, interpretable signals:

- **Relative exposure (SASA):** solvent accessibility of the residue, normalised within the chain.
- **Flexibility:** chain-normalised B-factor, used as a proxy for local mobility.
- **Secondary-structure penalty:** soft penalty for residues in α-helices or β-sheets, assigned with MDTraj/DSSP.
- **Nearby Lys boost:** rewards exposed lysines within the configured distance threshold.
- **Existing Cys penalty:** reduces the score when nearby cysteines are already present.

The final engineering score is:

```text
0.30 × ML stability
+ 0.20 × relative SASA
+ 0.20 × flexibility
+ 0.10 × nearby Lys boost
− 0.10 × existing Cys penalty
− 0.10 × secondary-structure penalty
```

Weights are configurable in `configs/default.yaml`. Secondary structure is assigned with MDTraj/DSSP; unknown assignments receive no penalty. Protected residues can be supplied as an optional exclusion annotation, but are not part of the default MVP score. Implementation details are in [docs/RANKING_FORMULA.md](docs/RANKING_FORMULA.md). Scores are intended for ranking candidates, not as calibrated probabilities.

## Reproducibility

```bash
.venv/bin/pip install -e '.[dev]'
.venv/bin/pytest -q
.venv/bin/ruff check .
```

GitHub Actions also runs the Python 3.10/3.12 test matrix, package build, portfolio notebook and Streamlit health check.

## Repository guide

- [Streamlit entry point](streamlit_app.py)
- [Portfolio notebook](notebooks/CysMutML_Portfolio_Demo.ipynb)
- [Model card](MODEL_CARD.md)
- [Homology validation](docs/HOMOLOGY_VALIDATION.md)
- [Feature schema](docs/FEATURE_SCHEMA.md)
- [Ranking formula](docs/RANKING_FORMULA.md)
- [Scientific audit](docs/SCIENTIFIC_AUDIT.md)
- [Retrospective validation](validation/godoy2011/VALIDATION_REPORT.md)
- [Versioned benchmark results](results/homology_validation/)

## Limitations

CysMutML does not predict immobilisation yield, retained activity, cysteine reactivity, disulfide formation or experimental success. FireProtDB measurements are heterogeneous, B-factors are only a rigidity proxy, and the ranking weights have not been experimentally calibrated.

The repository also contains an exploratory structure-trained ablation under `results/structural_ablation/`. It is retained for transparency but is not part of the deployed pipeline.

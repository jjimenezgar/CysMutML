# Decision Log

## Decision

Use a hybrid production architecture: physicochemical ML for mutation destabilization plus target-PDB structural heuristic for cysteine engineering ranking.

### Alternatives considered

- Train a structure-based DDG model requiring residue-to-PDB mapping for training.
- Use only ML without structural ranking.
- Use only structural heuristics without FireProtDB learning.

### Reason

FireProtDB provides hundreds of thousands of mutation observations, but only a much smaller subset can be mapped conservatively to structures. Training on physicochemical descriptors keeps the ML dataset large and reproducible. Structure remains important for cysteine engineering, so it is applied deterministically to the target PDB through SASA, B-factor flexibility, protected-site distance, and existing-cysteine proximity.

### Risk

The ML model cannot directly learn structural context. The final score is a heuristic, not an experimentally calibrated immobilization predictor.

## Decision

Use mutation-level median aggregation as the primary training dataset for the deployed model.

### Alternatives considered

- Train on all measurement-level rows.
- Drop duplicate mutation measurements.
- Model experimental conditions explicitly.

### Reason

Aggregation reduces overweighting of repeatedly measured protein/mutation pairs while preserving one interpretable target per mutation. Median aggregation is robust to outliers.

### Risk

Condition-specific variation is compressed into one value. The project preserves this as a limitation rather than pretending all experimental conditions are equivalent.

## Decision

Keep Ridge as the deployed model after the aggregated physicochemical comparison.

### Alternatives considered

- Promote HistGradientBoostingRegressor.
- Keep only the earlier measurement-level Ridge model.

### Reason

HGB had slightly better grouped CV metrics, but Ridge was close, simpler, faster, and easier to interpret. The project philosophy prioritizes scientific clarity over small performance gains.

### Risk

Ridge may miss nonlinear mutation-property interactions captured by HGB.

## FireProtDB Selection
Decision: Use FireProtDB as the primary dataset.
Alternatives considered: ProTherm, ProtaBank, ThermoMutDB, S669-only.
Reason: FireProtDB has public API documentation and curated mutation stability data.
Risk: Export fields and database composition can change over time.

## Target Definition
Decision: Use regression target `destabilization_ddg_kcal_mol`.
Alternatives considered: binary stabilizing/destabilizing label.
Reason: DDG is continuous and preserves more information.
Risk: Experimental DDG values are heterogeneous across conditions.

## Sign Convention
Decision: FireProtDB DDG is used directly because FireProtDB documents positive as destabilizing.
Alternatives considered: sign flipping.
Reason: Silent sign flips are scientifically dangerous.
Risk: Other datasets need source-specific normalization.

## Feature Selection
Decision: Start with interpretable physicochemical and structural descriptors.
Alternatives considered: ESM or graph neural networks.
Reason: v1 prioritizes interpretability and reproducibility.
Risk: Simple features may underperform richer representations.

## Validation
Decision: Use GroupKFold by `protein_id`.
Alternatives considered: random mutation split.
Reason: Prevents leakage across mutations from the same protein.
Risk: Protein names are not perfect homology clusters.

## Structure Mapping
Decision: Require sequence alignment and WT identity verification.
Alternatives considered: assume PDB numbering equals database numbering.
Reason: Numbering mismatches are common.
Risk: Many records are excluded or remain unmapped unless metadata is available.

## Models
Decision: Implement Dummy, Ridge, Random Forest, and histogram gradient boosting.
Alternatives considered: XGBoost, neural networks.
Reason: Classical models are adequate for a defensible v1.
Risk: Full random forest CV can be slow on large exports.

## Engineering Ranking
Decision: Keep ranking as a separate transparent heuristic.
Alternatives considered: train suitability directly.
Reason: No large supervised immobilization dataset is available.
Risk: Ranking weights are not experimentally optimized.

## Duplicate Handling
Decision: Preserve measurement-level rows and create a separate median-aggregated mutation dataset.
Alternatives considered: drop duplicates or average all values in-place.
Reason: repeated measurements can reflect real condition/assay differences.
Risk: median aggregation hides condition-specific effects.

## Outlier Handling
Decision: audit extreme DDG values but retain them in the main dataset.
Alternatives considered: percentile trimming for main training.
Reason: extreme values may be scientifically valid.
Risk: simple models are strongly challenged by extreme measurements.

## Structural Mapping Threshold
Decision: accept only exact PDB-chain subsequence mappings to UniProt during the executed partial run, plus WT identity verification.
Alternatives considered: global pairwise alignment for all chains.
Reason: full pairwise alignment was too slow interactively; exact subsequence mapping is conservative and reproducible.
Risk: many valid but non-exact engineered constructs are excluded.

## Structural Ablation Design
Decision: compare physicochemical-only and physicochemical-plus-structure on identical mapped rows and frozen GroupKFold folds.
Alternatives considered: compare full physicochemical data against structural subset.
Reason: different rows would make the ablation scientifically invalid.
Risk: the executed subset is small, so conclusions are preliminary.

## Candidate-Model Promotion
Decision: do not overwrite `models/cysmutml_model.joblib`.
Alternatives considered: promote the structural model.
Reason: partial ablation did not show improvement.
Risk: deployed model remains a simple physicochemical Ridge baseline.

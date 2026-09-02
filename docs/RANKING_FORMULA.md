# CysMutML ranking score

The ranking is a deterministic heuristic built on top of the deployed ML prediction. It is not trained on immobilization outcomes and it is not a probability.

## Inputs

For each possible X→Cys substitution, the pipeline computes:

- **ML stability score**: the predicted destabilization DDG mapped to [0, 1]. Higher means a more favourable stability estimate.
- **Relative exposure**: relative SASA of the mutated residue in the input structure, clipped to [0, 1].
- **Flexibility**: a min–max scaling of the chain-normalized mean B-factor. Higher values identify locally more flexible positions.
- **Nearby Lys boost**: a saturating function of exposed lysines within 20 Å. It represents a simple proxy for a possible multipoint attachment environment.
- **Nearby Cys penalty**: a penalty for proximity to existing cysteines. It is an engineering caution, not a prediction of disulfide formation.
- **Protected-site penalty**: an optional penalty for user-supplied protected residues within 8 Å.

## Final score

The default weights are stored in `configs/default.yaml`:

```text
final_priority =
    0.50 * ML stability
  + 0.20 * relative exposure
  + 0.15 * flexibility
  + 0.15 * nearby Lys boost
  - 0.10 * nearby Cys penalty
  - 0.10 * protected-site penalty
```

The score is used only to order candidates. A higher value means that the candidate is more attractive under these explicit assumptions.

## Interpretation

The model contribution comes from FireProtDB. SASA and B-factor are calculated from the target structure. The lysine boost and cysteine penalty are simple structural heuristics. None of these terms proves that a mutation will improve immobilization, activity, reactivity or experimental yield.

The CSV export retains the component columns so the final score can be reconstructed. The Streamlit app shows the main terms using user-facing labels and keeps implementation fields out of the default table.

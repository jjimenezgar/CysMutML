# CysMutML ranking score

The ranking is a deterministic heuristic built on top of the deployed ML prediction. It is not trained on immobilization outcomes and it is not a probability.

## Inputs

For each possible X→Cys substitution, the pipeline computes:

- **ML stability score**: predicted destabilization DDG mapped to [0, 1]. Higher means a more favourable stability estimate.
- **Relative exposure**: relative SASA of the mutated residue in the input structure, clipped to [0, 1].
- **Flexibility**: min–max scaling of the chain-normalized mean B-factor. Higher values identify locally more flexible positions.
- **Nearby Lys boost**: a saturating function of exposed lysines within 20 Å.
- **Nearby Cys penalty**: a penalty for proximity to existing cysteines. It is an engineering caution, not a prediction of disulfide formation.

Protected residues can still be supplied through the CLI as an optional exclusion annotation. They are retained in exported files for traceability, but are not part of the default MVP score.

## Final score

The default weights are stored in `configs/default.yaml`:

```text
final_priority =
    0.30 * ML stability
  + 0.25 * relative exposure
  + 0.25 * flexibility
  + 0.10 * nearby Lys boost
  - 0.10 * nearby Cys penalty
```

Positive terms sum to 0.90; the remaining 0.10 is the maximum penalty contribution. The score is used only to order candidates.

## Interpretation

The model contribution comes from FireProtDB. SASA and B-factor are calculated from the target structure. The lysine boost and cysteine penalty are simple structural heuristics. None of these terms proves that a mutation will improve immobilization, activity, reactivity or experimental yield.

The CSV export retains the component columns so the final score can be reconstructed. The Streamlit app shows the main terms using user-facing labels and keeps implementation fields out of the default table.

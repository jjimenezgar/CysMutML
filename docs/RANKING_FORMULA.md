# CysMutML ranking score

The ranking is a deterministic heuristic built on top of the deployed ML prediction. It is not trained on immobilization outcomes and it is not a probability.

## Inputs

For each possible X→Cys substitution, the pipeline computes:

- **ML stability score**: predicted destabilization DDG mapped to [0, 1]. Higher means a more favourable stability estimate.
- **Relative exposure**: relative SASA of the mutated residue in the input structure, clipped to [0, 1].
- **Flexibility**: min–max scaling of the chain-normalized mean B-factor. Higher values identify locally more flexible positions.
- **Nearby Lys boost**: a saturating function of exposed lysines within 20 Å.
- **Nearby Cys penalty**: a penalty for proximity to existing cysteines.
- **Secondary-structure penalty**: a soft penalty for residues assigned by MDTraj/DSSP to an α-helix or β-sheet. Loops and unknown assignments are not penalized.

## Final score

The default weights are stored in `configs/default.yaml`:

```text
final_priority =
    0.30 * ML stability
  + 0.20 * relative exposure
  + 0.20 * flexibility
  + 0.10 * nearby Lys boost
  - 0.10 * nearby Cys penalty
  - 0.10 * secondary-structure penalty
```

The score is used only to order candidates. A higher value means that the candidate is more attractive under these explicit assumptions.

## Interpretation

The model contribution comes from FireProtDB. SASA, B-factor and secondary structure are calculated from the target structure. The lysine boost and cysteine penalties are simple structural heuristics. The secondary-structure term is a soft prior, not evidence that every helix or sheet position is unsuitable. None of these terms proves that a mutation will improve immobilization, activity, reactivity or experimental yield.

The CSV export retains the component columns so the final score can be reconstructed. The Streamlit app shows the main terms using user-facing labels and keeps implementation fields out of the default table.

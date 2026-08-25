# ML vs Heuristic in CysMutML

## What Is Learned From Data?

CysMutML learns a mutation-associated destabilization model from FireProtDB experimental stability measurements.

The production ML model uses physicochemical mutation descriptors only:

- WT residue identity and properties;
- mutant residue identity and properties;
- mutant-minus-WT property changes;
- BLOSUM62 substitution score.

The model output is:

```text
predicted_destabilization_ddg
```

This is a continuous stability-risk estimate. It is not a feasibility probability.

## What Is Calculated From Structure?

For a new target PDB, CysMutML calculates structural descriptors for each possible X->Cys candidate:

- relative SASA;
- B-factor based flexibility proxy;
- distance to user-supplied protected residues, if provided;
- distance to existing cysteines, if present.

These values are not used to train the deployed ML model.

## What Is Heuristic?

The final `cys_suitability_score` is a transparent engineering heuristic:

```text
score =
  stability_weight * stability_component
+ accessibility_weight * accessibility_component
+ rigidity_weight * rigidity_component
- protected_penalty_weight * protected_site_penalty
- existing_cys_penalty_weight * existing_cys_penalty
```

Default weights are in `configs/default.yaml`.

The normalized components are:

- `stability_component`: clipped linear transform of predicted destabilization, high is favorable;
- `accessibility_component`: clipped relative SASA, high is more exposed;
- `rigidity_component`: within-structure min-max transform of the flexibility proxy, high is less flexible by the selected proxy.

These weights are not experimentally optimized.

## What Is Not Claimed?

CysMutML does not claim to predict:

- immobilization yield;
- catalytic activity after mutation;
- calibrated probability of experimental success;
- disulfide formation;
- a state-of-the-art structure-based DDG.

The project deliberately separates learned mutation tolerance from structure-based engineering prioritization.

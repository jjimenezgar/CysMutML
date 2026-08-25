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
- B-factor-based flexibility proxy;
- local exposed Lys count within a configurable radius;
- distance to user-supplied protected residues, if provided;
- distance/counts to existing cysteines, if present.

These values are not used to train the deployed ML model.

## What Is Heuristic?

The ranking now separates three transparent heuristic outputs:

```text
cys_site_suitability
rigidification_potential
final_engineering_score
```

Default weights are in `configs/default.yaml`.

The normalized components are:

- `stability_component`: clipped linear transform of predicted destabilization, high is favorable;
- `accessibility_component`: clipped relative SASA, high is more exposed;
- `flexibility_component`: within-structure min-max transform of the B-factor proxy, high is more flexible by the selected proxy;
- `lysine_environment_component`: saturating transform of nearby exposed Lys count;
- `existing_cys_penalty`: native Cys proximity caution;
- `protected_site_penalty`: optional user-defined protected-residue caution.

These weights are not experimentally optimized.

`cys_site_suitability` is intended to represent mutation tolerance and practical Cys exposure. `rigidification_potential` is intended to represent local opportunity for multipoint immobilization/rigidification. `final_engineering_score` combines both. None of these is an ML prediction.

## What Is Not Claimed?

CysMutML does not claim to predict:

- immobilization yield;
- catalytic activity after mutation;
- calibrated probability of experimental success;
- disulfide formation;
- a state-of-the-art structure-based DDG.

The project deliberately separates learned mutation tolerance from structure-based engineering prioritization.

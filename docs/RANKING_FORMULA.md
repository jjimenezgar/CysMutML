# CysMutML v1.0 Ranking Formula

The final `cys_suitability_score` is a deterministic engineering heuristic. It is not learned from immobilization data and is not a calibrated probability.

## Predicted DDG

The ML model predicts:

```text
predicted_destabilization_ddg
```

Convention:

```text
larger positive values = greater predicted destabilization
```

## DDG to Stability Component

Default references from `configs/default.yaml`:

```text
favorable_ddg = -1.0
unfavorable_ddg = 2.0
```

Formula:

```text
stability_component =
  clip(1 - (predicted_destabilization_ddg - favorable_ddg)
           / (unfavorable_ddg - favorable_ddg), 0, 1)
```

Interpretation:

- `1.0` means favorable by this heuristic transform.
- `0.0` means unfavorable by this heuristic transform.
- Values are clipped to `[0, 1]`.
- This is not a probability of success.

Worked examples:

| Predicted DDG | Calculation | Stability component |
|---:|---|---:|
| -1.0 | `1 - (-1 - -1) / 3` | 1.000 |
| 0.5 | `1 - (0.5 - -1) / 3` | 0.500 |
| 2.0 | `1 - (2 - -1) / 3` | 0.000 |
| 3.0 | clipped below 0 | 0.000 |

## SASA and Accessibility

SASA is calculated with Biopython Shrake-Rupley on the input structure.

Relative SASA:

```text
relative_sasa = absolute_residue_sasa / max_asa_for_residue_type
```

The maximum ASA reference table is Tien et al. 2013, stored in `cysmutml.amino_acids`.

Accessibility component:

```text
accessibility_component = clip(relative_sasa, 0, 1)
```

High relative SASA means the residue is more solvent exposed in the static input structure. It does not guarantee cysteine chemistry or immobilization success.

## B-Factor-Derived Rigidity Proxy

The v1.0 rigidity component uses the chain-normalized residue B-factor:

```text
normalized_b_factor = (residue_mean_b_factor - chain_mean_b_factor) / chain_b_factor_std
```

The ranking then maps this flexibility proxy to:

```text
rigidity_component = minmax_low_good(flexibility_value)
```

Lower B-factor values receive higher `rigidity_component` within the target protein. If fewer than two valid values are available, the fallback is neutral:

```text
rigidity_component = 0.5
```

B-factors are not pure dynamics measurements. They depend on crystallographic refinement, disorder, occupancy, resolution, and model quality.

## Protected-Site Penalty

Protected residues are optional user-supplied residues such as:

```text
A:45,A:48,B:120
```

Distances use C-alpha coordinates. Protected residues can be on any chain in the input structure.

Default radius:

```text
protected_site_radius_angstrom = 8.0
```

Formula:

```text
protected_site_penalty =
  clip((radius - distance_to_nearest_protected) / radius, 0, 1)
```

If no protected residues are supplied, the penalty is `0.0`.

## Existing-Cys Warning and Penalty

Distances to existing cysteines use C-alpha coordinates across all chains in the input structure.

Default warning radius:

```text
existing_cys_warning_radius_angstrom = 6.0
```

Default penalty radius:

```text
existing_cys_penalty_radius_angstrom = 6.0
```

Formula:

```text
existing_cys_penalty =
  clip((radius - distance_to_existing_cys) / radius, 0, 1)
```

Proximity to an existing cysteine is an engineering caution. It does not imply disulfide formation.

## Final Score

Default v1.0 weights:

```text
stability_weight = 0.50
accessibility_weight = 0.30
rigidity_weight = 0.20
protected_penalty_weight = 0.10
existing_cys_penalty_weight = 0.05
```

Exact implemented formula:

```text
cys_suitability_score =
  0.50 * stability_component
+ 0.30 * accessibility_component
+ 0.20 * rigidity_component
- 0.10 * protected_site_penalty
- 0.05 * existing_cys_penalty
```

All component columns are preserved in `residue_ranking.csv`, so the score can be reconstructed.

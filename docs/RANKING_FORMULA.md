# CysMutML Ranking Formula

The final engineering ranking is a deterministic heuristic. It is not learned from immobilization data and is not a calibrated probability.

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

## B-Factor-Derived Flexibility Proxy

The current rigidification heuristic uses the chain-normalized residue B-factor:

```text
normalized_b_factor = (residue_mean_b_factor - chain_mean_b_factor) / chain_b_factor_std
```

For rigidification, higher local flexibility is treated as potentially more improvable by multipoint attachment:

```text
flexibility_component = minmax_high_good(normalized_b_factor)
```

If fewer than two valid values are available, the fallback is neutral:

```text
flexibility_component = 0.5
```

B-factors are not pure dynamics measurements. They depend on crystallographic refinement, disorder, occupancy, resolution, and model quality.

The legacy `rigidity_component` column is preserved for backward compatibility, but it is not the main rigidification term.

## Local Exposed Lys Environment

The rigidification heuristic includes the number of exposed Lys residues near the candidate Cys.

Defaults:

```text
radius = 20 Angstrom
exposed Lys threshold = relative_sasa >= 0.25
saturation_k = 3.0
```

Formula:

```text
lysine_environment_component =
  local_exposed_lys_count / (local_exposed_lys_count + saturation_k)
```

This is a saturating boost, not a fitted probability. It reflects that glyoxyl-style multipoint attachment often benefits from nearby accessible amino groups, but the exact chemistry remains system-specific.

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
existing_cys_warning_radius_angstrom = 10.0
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

The current ranking separates two intermediate scores.

### Cys Site Suitability

```text
cys_site_suitability =
  0.60 * stability_component
+ 0.35 * accessibility_component
- 0.10 * existing_cys_penalty
- 0.15 * protected_site_penalty
```

This score asks whether the site is a reasonable Cys mutation candidate.

### Rigidification Potential

```text
rigidification_potential =
  0.35 * flexibility_component
+ 0.40 * lysine_environment_component
+ 0.25 * accessibility_component
- 0.05 * existing_cys_penalty
- 0.10 * protected_site_penalty
```

This score asks whether the local structural environment may be practically useful for rigidifying immobilization chemistry.

### Final Engineering Score

```text
final_engineering_score =
  0.60 * cys_site_suitability
+ 0.40 * rigidification_potential
```

`cys_suitability_score` remains as an alias for `final_engineering_score` for backward compatibility. All component columns are preserved in `residue_ranking.csv`, so every score can be reconstructed.

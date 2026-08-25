# Duplicate Measurement Summary

Repeated protein/mutation measurements are preserved in the measurement-level dataset.
The aggregated dataset uses the median `destabilization_ddg_kcal_mol`
as the primary target.

- Measurement-level rows: 555,932
- Aggregated mutation rows: 352,005
- Repeated protein/mutation groups: 47,223
- Repeated measurement rows: 251,150

## Duplicate Categories

- replicate_or_unclear_repeat: 46,279
- condition_or_assay_variation: 814
- exact_or_near_exact_duplicate: 130

Categories are heuristic audit labels, not automatic exclusion rules.
No measurements are removed solely because repeated observations differ.

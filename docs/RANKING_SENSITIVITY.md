# Ranking Sensitivity

The v1.0 ranking weights are heuristic defaults, not experimentally optimized parameters.

For the real 1CSP chain A case study, the ranking was recalculated under three scenarios:

| Scenario | Stability | Accessibility | Rigidity |
|---|---:|---:|---:|
| balanced | 0.50 | 0.30 | 0.20 |
| stability-heavy | 0.70 | 0.20 | 0.10 |
| accessibility-heavy | 0.35 | 0.50 | 0.15 |

Penalties were kept fixed:

```text
protected penalty weight = 0.10
existing-cys penalty weight = 0.05
```

## Results

Compared with the balanced ranking:

| Scenario | Top-5 overlap | Top-10 overlap | Spearman rank correlation |
|---|---:|---:|---:|
| balanced | 5 | 10 | 1.000 |
| stability-heavy | 4 | 7 | 0.947 |
| accessibility-heavy | 2 | 4 | 0.692 |

## Interpretation

The ranking is fairly stable when the ML stability component is emphasized. It is more sensitive when accessibility receives the largest weight. This is expected because exposed residues can differ substantially from residues favored by mutation physicochemistry.

This analysis does not optimize weights. It only shows how the heuristic ranking changes under reasonable alternative assumptions.

Files:

- `results/ranking_sensitivity/sensitivity_summary.csv`
- `results/ranking_sensitivity/top10_by_scenario.csv`
- `results/ranking_sensitivity/balanced_ranking.csv`
- `results/ranking_sensitivity/stability_heavy_ranking.csv`
- `results/ranking_sensitivity/accessibility_heavy_ranking.csv`

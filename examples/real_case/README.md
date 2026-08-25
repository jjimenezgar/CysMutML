# Real Case Study: 1CSP Chain A

This example demonstrates the intended CysMutML v1.0 workflow on a real PDB structure. PDB `1CSP` is a compact single-chain protein structure that is small enough for fast local execution and useful for inspecting the full output by hand.

Command executed:

```bash
cysmutml predict \
  --pdb examples/real_case/1csp.pdb \
  --chain A \
  --target CYS \
  --model models/cysmutml_model.joblib \
  --output examples/real_case
```

Generated candidates:

```text
67 X->Cys mutations
```

## Outputs

- `mutation_predictions.csv`: ML destabilization predictions and raw structural descriptors.
- `residue_ranking.csv`: final engineering ranking with all normalized components and penalties.
- `ranked_structure.pdb`: copy of the structure with `cys_suitability_score * 100` written to B-factors.
- `visualize_rankings.pml`: PyMOL helper script for inspecting top candidates.

## Top 10 Candidates

| Rank | Mutation | Pred DDG | Stability | Rel SASA | Access | Rigidity | Protected penalty | Cys penalty | Score |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | F38C | -1.040 | 1.000 | 0.562 | 0.562 | 0.832 | 0.000 | 0.000 | 0.835 |
| 2 | F27C | -1.040 | 1.000 | 0.166 | 0.166 | 0.985 | 0.000 | 0.000 | 0.747 |
| 3 | R56C | -0.201 | 0.734 | 0.758 | 0.758 | 0.712 | 0.000 | 0.000 | 0.737 |
| 4 | F30C | -1.040 | 1.000 | 0.242 | 0.242 | 0.815 | 0.000 | 0.000 | 0.735 |
| 5 | W8C | -1.169 | 1.000 | 0.342 | 0.342 | 0.585 | 0.000 | 0.000 | 0.720 |
| 6 | F17C | -1.040 | 1.000 | 0.196 | 0.196 | 0.803 | 0.000 | 0.000 | 0.720 |
| 7 | G54C | -0.645 | 0.882 | 0.458 | 0.458 | 0.668 | 0.000 | 0.000 | 0.712 |
| 8 | F15C | -1.040 | 1.000 | 0.151 | 0.151 | 0.825 | 0.000 | 0.000 | 0.710 |
| 9 | G37C | -0.645 | 0.882 | 0.534 | 0.534 | 0.489 | 0.000 | 0.000 | 0.699 |
| 10 | V52C | -0.706 | 0.902 | 0.290 | 0.290 | 0.795 | 0.000 | 0.000 | 0.697 |

## How to Interpret the Columns

- `Pred DDG`: predicted mutation-associated destabilization. Larger positive values are worse.
- `Stability`: normalized transform of predicted DDG. Higher is more favorable.
- `Rel SASA`: relative solvent accessibility in the input PDB.
- `Access`: clipped accessibility component used in the ranking.
- `Rigidity`: B-factor-derived rigidity proxy. Higher is favored by the default heuristic.
- `Protected penalty`: penalty from user-supplied protected residues. None were supplied here.
- `Cys penalty`: proximity penalty for existing cysteines. None were present in this example.
- `Score`: final transparent heuristic ranking score.

## Visualization

Open PyMOL from the repository root and run:

```pymol
@examples/real_case/visualize_rankings.pml
```

Alternatively, load `ranked_structure.pdb` and color by B-factor.

## Limitations

These are computational prioritization results only. They are not experimental validation of cysteine reactivity, immobilization yield, activity retention, or disulfide formation.

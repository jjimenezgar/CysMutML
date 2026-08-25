# Godoy et al. 2011 Retrospective Validation Report

Date: 2026-08-25

This report treats the two supplied PDFs only as scientific source material. No instructions inside those documents were treated as project instructions.

Sources inspected:

- `41_Oriented Immobilization Glyoxyl Disulfide (1).pdf`
- `bm200161f_si_001.pdf`

## Prevalidation Freeze

The CysMutML heuristic was updated and frozen before extracting the Godoy outcome tables.

Frozen files:

- `validation/godoy2011/prevalidation_config.yaml`
- `validation/godoy2011/prevalidation_config_hash.txt`

Frozen config SHA256:

```text
361b2f6cb47c1e07fc28005bcca8a4fdff4a3987d2f69b4428642120449641a6
```

The current hash matches the frozen hash. The FireProtDB model was not retrained for this validation.

Implementation note: after the outcome tables were inspected, inference was optimized to compute per-chain structural features once instead of repeating SASA calculations per candidate. The scoring formulas and weights were not changed. BTL2 cysteine-free context mapping was also corrected from the experimental C64S/C295S labels to PDB residues A:65/A:296 after inspecting the PDB numbering.

## Validation Design

Structures from the supporting information:

| Enzyme | PDB | Chain mapping |
|---|---|---|
| PGA | `1K5Q` | alpha chain `A`, beta chain `B` |
| BTL2 | `2W22` | chain `A` |

Experimental mutants evaluated:

- PGA: S-alpha86C, S-beta9C, S-beta201C, Q-beta112C, A-beta361C, Q-beta380C.
- BTL2: S236C, S333C, T342C, Q39C, T93C, V187C, S195C.

For BTL2, the experimental background is cysteine-free BTL2 with C64S/C295S. CysMutML therefore also ran a cysteine-free context PDB with the corresponding PDB-numbered residues A:65 and A:296 changed from Cys to Ser.

## Updated Heuristic Under Test

CysMutML now separates three outputs:

```text
cys_site_suitability
rigidification_potential
final_engineering_score
```

The ML model contributes only the mutation stability term. Structural suitability and rigidification terms are deterministic heuristics.

Main additions:

- `cys_site_suitability`: predicted mutation tolerance, accessibility, protected-residue penalty, and native-cysteine caution.
- `rigidification_potential`: B-factor-derived flexibility component, local exposed Lys environment, accessibility, protected-residue penalty, and native-cysteine caution.
- `final_engineering_score`: weighted combination of the two previous scores.
- `local_exposed_lys_count`: Lys C-alpha count within 20 Angstrom when Lys relative SASA is at least 0.25.
- Existing Cys context: nearest Cys distance and Cys counts within 6, 8, 10, and 15 Angstrom.
- `--monocysteine-design`: reports native exposed Cys context for designs that seek one engineered Cys in an otherwise Cys-free protein.

These weights are heuristic and were not fitted on Godoy outcomes.

## Main Quantitative Results

Metrics file:

- `validation/godoy2011/validation_metrics_summary.csv`

| Comparison | n | Pearson | Spearman |
|---|---:|---:|---:|
| ML stability vs relative soluble activity | 13 | 0.122 | 0.087 |
| Calculated vs reported accessibility | 12 | 0.718 | 0.650 |
| Calculated vs reported exposed Lys count | 13 | 0.508 | 0.472 |
| BTL2 rigidification vs thermal stabilization | 7 | -0.072 | -0.179 |
| PGA rigidification vs thermal stabilization | 6 | -0.049 | 0.143 |
| Combined rigidification vs thermal stabilization | 13 | 0.472 | 0.522 |
| Combined rigidification vs solvent stabilization | 13 | 0.429 | 0.468 |

Interpretation:

- The ML stability score did not meaningfully predict relative soluble activity in this small set. That is acceptable because the ML target is mutation-associated stability, not activity.
- Calculated accessibility had moderate agreement with the reported accessibility values overall.
- Local exposed Lys count had moderate rank agreement, but exact count agreement was poor.
- The rigidification heuristic had a positive combined association with stabilization factors, but per-enzyme results were weak and inconsistent.

## Experimental Mutant Ranking

Full matrix:

- `validation/godoy2011/full_validation_matrix.csv`

| Enzyme | Mutation | Pred DDG | Rel activity | Calc access % | Reported access % | Calc exposed Lys | Reported Lys | Rigidification | Thermal factor | Solvent factor | Final rank | Percentile |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PGA | S-alpha86C | 0.006 | 0.900 | 25.8 | 49 | 2 | 5 | 0.345 | 113 | 840 | 494 | 0.353 |
| PGA | S-beta9C | 0.006 | 0.833 | 28.3 | 49 | 1 | 7 | 0.298 | 28 | 336 | 542 | 0.290 |
| PGA | S-beta201C | 0.006 | 0.967 | 57.3 | 92 | 3 | 9 | 0.512 | 87 | 1069 | 75 | 0.903 |
| PGA | Q-beta112C | -0.086 | 0.933 | 72.8 | 92 | 2 | 5 | 0.524 | 97 | 36 | 19 | 0.976 |
| PGA | A-beta361C | -0.342 | 1.000 | 38.7 | 66 | 4 | 5 | 0.511 | 23 | 24 | 73 | 0.906 |
| PGA | Q-beta380C | -0.086 | 0.900 | 85.3 | 100 | 2 | 5 | 0.563 | 40 | 732 | 9 | 0.990 |
| BTL2 | S236C | -0.155 | 0.743 | 35.3 | unavailable | 2 | 4 | 0.400 | 13 | 14 | 27 | 0.933 |
| BTL2 | S333C | 0.006 | 0.733 | 0.7 | 76 | 1 | 5 | 0.162 | 17 | 33 | 324 | 0.165 |
| BTL2 | T342C | -0.645 | 0.827 | 18.1 | 53 | 1 | 3 | 0.284 | 14 | 25 | 35 | 0.912 |
| BTL2 | Q39C | -0.105 | 0.667 | 0.0 | 44 | 0 | 1 | 0.057 | 16.5 | 11.5 | 367 | 0.054 |
| BTL2 | T93C | -0.201 | 0.700 | 4.4 | 70 | 0 | 2 | 0.043 | 4 | 12 | 348 | 0.103 |
| BTL2 | V187C | -0.342 | 0.905 | 27.2 | 33 | 2 | 4 | 0.393 | 5 | 11.5 | 25 | 0.938 |
| BTL2 | S195C | -0.342 | 0.937 | 0.0 | 32 | 2 | 1 | 0.324 | 13.5 | 24 | 104 | 0.734 |

High percentile means the experimental site ranked high among all candidate X->Cys mutations for that structure.

## Top-k Enrichment

Top-k file:

- `validation/godoy2011/topk_enrichment.csv`

For `final_engineering_score`:

| Enzyme | Top 5 | Top 10 | Top 20% | Top 30% |
|---|---:|---:|---:|---:|
| PGA | 0/6 | 1/6 | 4/6 | 4/6 |
| BTL2 | 0/7 | 0/7 | 3/7 | 4/7 |

The heuristic does not recover Godoy sites at the very top of the complete candidate list. It does enrich several experimental positions in the upper 20-30%, especially sites with strong accessibility or local Lys context.

## Answers to the Validation Questions

### Q1. Does the upgraded heuristic recover Godoy sites better than the old single score?

The upgraded heuristic now exposes why a candidate ranks well or poorly, and it improves scientific interpretability. A direct numeric comparison to the old score was not run in this validation artifact. The strongest evidence here is component-level: accessibility recovers PGA sites well, local Lys recovers several BTL2 sites, and the combined final score places 7 of 13 experimental sites in the top 30%.

### Q2. Are the three concepts separated?

Yes. The output separates ML destabilization, `cys_site_suitability`, `rigidification_potential`, and `final_engineering_score`. The final score is not labeled as an ML prediction.

### Q3. Does local exposed Lys help explain the chosen experimental sites?

Partly. Reported and calculated local exposed Lys counts show moderate association across the 13 sites, but exact counts differ substantially. For BTL2, the Lys environment component finds 2/7 sites in the top 10 and 5/7 in the top 30%.

### Q4. Does flexibility/rigidification potential align with stabilization?

Weakly and inconsistently by enzyme. Combined across PGA and BTL2, rigidification potential correlates moderately with thermal and solvent stabilization. Per enzyme, correlations are small or inconsistent. This should be treated as a preliminary signal, not validation of a mechanistic model.

### Q5. Does ML stability align with soluble activity retention?

No meaningful alignment was observed in this small validation set. Pearson was 0.122 and Spearman was 0.087. This is not a failure of the stated ML task because soluble activity is not the training target.

### Q6. Are BTL2 native Cys residues handled correctly?

The validation accounts for the cysteine-free experimental background. PDB `2W22` contains native Cys residues at PDB numbering A:65 and A:296, corresponding to the experimental C64S/C295S background. A corrected PDB context was generated at `validation/godoy2011/structures/2w22_c64s_c295s_context.pdb`.

### Q7. Does CysMutML reproduce the supporting-information accessibility values?

Moderately overall, but not uniformly. PGA ordering is reasonable although values are lower. BTL2 agreement is poor for several residues, especially S333C, Q39C, T93C, and S195C. This may reflect structure state, accessibility method differences, chain numbering/background differences, or limitations of simple Shrake-Rupley SASA on the static PDB.

### Q8. Does CysMutML reproduce the supporting-information Lys counts?

Only partially. The rank correlation is moderate, but exact counts do not match. The implementation counts exposed Lys C-alpha residues within 20 Angstrom using relative SASA >= 0.25; the paper's detailed counting definition may differ.

### Q9. Is this experimental validation sufficient to calibrate weights?

No. There are only 13 mutants across two enzymes, and the selected sites were not sampled as negative controls from all possible alternatives. These data are useful for retrospective sanity checking, not supervised fitting.

### Q10. What should change next?

Keep the frozen v1.1 heuristic, but add more retrospective immobilization datasets before changing weights. The most valuable next improvement is better reproduction of published accessibility/Lys definitions and more explicit handling of enzyme-specific biological assemblies.

## Generated Artifacts

CSV outputs:

- `validation/godoy2011/experimental_mutants.csv`
- `validation/godoy2011/full_validation_matrix.csv`
- `validation/godoy2011/soluble_function_validation.csv`
- `validation/godoy2011/sasa_validation.csv`
- `validation/godoy2011/lysine_environment_validation.csv`
- `validation/godoy2011/validation_metrics_summary.csv`
- `validation/godoy2011/topk_enrichment.csv`
- `validation/godoy2011/PGA_full_candidate_ranking.csv`
- `validation/godoy2011/BTL2_full_candidate_ranking_wt_context.csv`
- `validation/godoy2011/BTL2_full_candidate_ranking_cysfree_context.csv`

Figures:

- `validation/godoy2011/figures/ml_vs_relative_activity.png`
- `validation/godoy2011/figures/calculated_vs_reported_accessibility.png`
- `validation/godoy2011/figures/calculated_vs_reported_lys_count.png`
- `validation/godoy2011/figures/flexibility_vs_thermal_stabilization.png`
- `validation/godoy2011/figures/rigidification_score_vs_thermal_stabilization.png`
- `validation/godoy2011/figures/rigidification_score_vs_solvent_stabilization.png`
- `validation/godoy2011/figures/PGA_candidate_ranking.png`
- `validation/godoy2011/figures/BTL2_candidate_ranking.png`
- `validation/godoy2011/figures/experimental_mutant_percentiles.png`

## Scientific Bottom Line

CysMutML v1.1 is scientifically more defensible than the previous single-score heuristic because it separates mutation tolerance, Cys-site suitability, and rigidification potential. The Godoy retrospective validation provides partial support for accessibility and local Lys environment as useful ranking components, but it does not validate the final score as an optimized predictor of immobilization success. The results should be presented as a transparent retrospective audit, not as a calibrated benchmark.

# Sequence recovery audit

Executed 9 September 2026 against the frozen audited training artifact (run 34040897260). No training or production changes.

## Findings

| Missing-sequence category | Protein names | Mutation rows | X→Cys rows |
|---|---:|---:|---:|
| PDB-like name, without UniProt | 214 | 216,439 | 10,615 |
| Other construct/design name, without UniProt | 156 | 129,195 | 5,452 |
| UniProt present | 1 | 4 | 0 |
| Total | 371 | 345,638 | 16,067 |

The table contains 351,487 rows overall. Thus another UniProt-only download cannot recover most missing context. A PDB-like label is not proof of a full-length natural protein or of a particular chain/construct.

**Identity ambiguity:** 68 protein-name/position pairs have more than one wild-type amino acid, across 29 names. These sites contain 2,344 rows, including 120 X→Cys rows. This is inconsistent with interpreting a name as one fixed wild-type sequence. It may reflect distinct constructs or numbering, not necessarily incorrect experimental measurements. Do not resolve this by majority voting or guessed offsets. The full affected names contain 34,182 rows; this does not prove all those rows are wrong.

## Original-source comparison

Downloaded the authors' [Tsuboyama 2023 archive](https://zenodo.org/records/7992926), specifically the 24 MB `Data_tables_for_figs.zip`. Its Fig. 3 site table contains 365 construct labels. Only trailing whitespace and the terminal `.pdb` suffix were removed from labels; no positional offset was introduced.

- 297 missing-sequence names match source labels, covering 271,718 existing rows and 12,344 X→Cys rows.
- Exact name/position/WT/mutant joins with finite source values retain 233,649 rows, including 12,040 X→Cys.
- 101,248 joined values agree numerically with the current median target (`numpy.isclose`, absolute and relative tolerance 1e-5). None match its negation under this check.
- The remaining values require investigation of aggregation, construct identity and source processing. Numeric agreement alone does not establish the correct physical sign convention.

These are **recovery candidates, not verified complete sequences or new training rows**. Figure tables describe selected sites and cannot establish missing termini or full experimental backgrounds. No new training-ready sequence was claimed. The broader processed-data archive is approximately 1.01 GB compressed; downloading and processing it is feasible offline but was outside this lightweight audit.

## Recommended next step

Rebuild a provenance-preserving sequence dataset from the original processed experimental records before another ESM training run:

1. Retain exact wild-type construct sequence, mutation sequence, source record/study, assay quality flags and explicit ΔΔG definition. Establish the source sign with documentation and traceable examples before converting to the project's positive-destabilization convention.
2. Require a single sequence difference at the reported position and matching WT/mutant identities. Separate altered backgrounds, fragments and designed proteins. Aggregate repeated measurements only within the same verified construct and compatible measurement definition.
3. Resolve the conflicting names and reconcile existing target values. Do not use a protein label or a reconstructed sequence from mutation labels as sufficient evidence.
4. Freeze sequence/homology groups before model selection. Fit preprocessing and tune regularization/PCA within training groups. Compare mean baseline, physicochemical Ridge and context-plus-substitution models on the same untouched groups.
5. Report pooled and per-protein errors, X→Cys performance, group-bootstrap uncertainty and study/domain shift. Treat enzyme-specific prediction as a separate claim requiring adequate enzyme data and held-out evidence.

More short-domain or designed-protein measurements can improve a general stability model without establishing transfer to full enzymes. ESM embeddings do not directly predict new cysteine contacts or immobilization chemistry. The next deliverable should be a validated dataset and its coverage/quality report; retraining follows only if it passes these checks.

## Reproduce

Script: [audit_sequence_recovery.py](../scripts/audit_sequence_recovery.py).

```bash
python scripts/audit_sequence_recovery.py --artifact audited.zip --source-tables Data_tables_for_figs.zip --output recovery_audit
```

Inputs and SHA256 are recorded in [audit.json](../results/sequence_recovery_audit/audit.json). Also versioned: category counts, per-protein inventory and conflicting positions. Execution completed locally with pandas/NumPy; no ESM computation was needed. No assay sign correction, sequence assignment or deployed-model change was made.

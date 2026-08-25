# Data Acquisition

## FireProtDB
Primary data source: FireProtDB v2.0 public REST API.

Source pages inspected:

- https://loschmidt.chemi.muni.cz/fireprotdb/api-docs/
- https://loschmidt.chemi.muni.cz/fireprotdb/help/
- FireProtDB and FireProtDB 2.0 Nucleic Acids Research/PMC articles.

The implemented command uses the documented `/api/search` endpoint with CSV output and filters for mutant entries with non-empty DDG:

```bash
cysmutml prepare-data --download-fireprotdb \
  --raw data/raw/fireprotdb.csv \
  --output data/processed/fireprotdb_mutations.csv
```

Executed on 2026-08-24. The API export completed and wrote `data/raw/fireprotdb.csv`.

Observed live CSV fields included `SUBSTITUTION`, `DDG`, `PROTEIN`, `WWPDB`, `METHOD`, `MEASURE`, `PH`, `EXP_TEMPERATURE`, and `SOURCE_DATASET`.

FireProtDB documents DDG such that negative values denote stabilizing mutations and positive values denote destabilizing mutations. CysMutML therefore preserves the source value as:

```text
destabilization_ddg_kcal_mol = FireProtDB DDG
```

## Local Resume
If the API is unavailable, place a FireProtDB CSV export in `data/raw/fireprotdb.csv` and run:

```bash
cysmutml prepare-data --raw data/raw/fireprotdb.csv \
  --output data/processed/fireprotdb_mutations.csv
```

## S669
S669 support is not bundled because no legal, single official downloadable source was integrated during this run. A future external evaluation must check overlap against FireProtDB by protein/PDB, mutation, and preferably sequence identity before claiming independence.


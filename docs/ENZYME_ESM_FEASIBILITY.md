# Enzyme ESM feasibility


## Enzyme ESM feasibility — 9 September 2026

Successful run: https://github.com/jjimenezgar/CysMutML/actions/runs/34327672338

Using the frozen audited dataset, the strict reviewed-UniProt plus EC filter and exact sequence/WT checks retained 2,043 mutations from 64 unique sequences, including 46 X-to-Cys mutations. This is annotation/mapping coverage, not proof that excluded records are non-enzymes. No model was trained or promoted.

ESM-2 8M CPU pilot (two threads): 96/228/809 residues took 0.022/0.039/0.178 seconds for tokenization and inference after model load. Peak process RSS: 467.2 MiB; model load: 1.06 seconds in that runner. These timings exclude dependency installation and are not Streamlit latency guarantees. Full details and UniProt annotation snapshot are retained in the workflow artifact. The small Cys subset limits robust downstream validation; enzyme-only specialization remains experimental.

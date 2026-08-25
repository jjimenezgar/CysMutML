# Future ESM Extension

ESM is intentionally not mandatory for v1.

A future v2 could compare:

A. physicochemical descriptors;
B. structural descriptors;
C. ESM residue-level embeddings;
D. structural plus ESM features.

For a mutation WT->MUT, residue embeddings could be extracted from the WT sequence at the mutated position and concatenated with mutation delta features. More advanced variants could include local sequence windows or mutant-sequence embeddings.

Evaluation must keep the same leakage controls: protein-grouped or homology-clustered folds, with no use of an external benchmark during model selection.


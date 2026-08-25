# Data Quality Summary

Status: EXECUTED AND VERIFIED for FireProtDB preprocessing on 2026-08-24.

- Raw records: 613,208
- Single mutations with valid numeric DDG target: 555,932
- Rejected during parsing: 57,276
- Successful structure mappings: unavailable; bulk mapping not executed
- Failed mappings: unavailable; bulk mapping not executed
- Unique proteins: 542
- Unique structures/PDB fields: 175
- Unique protein/mutation pairs: 352,005
- X->Cys mutation rows: 25,026
- Missing PDB field rows: 547,955
- Duplicate protein/mutation measurements: 203,927
- Final physicochemical ML-ready rows: 555,932

Important caveat: the downloaded FireProtDB CSV did not include canonical sequences or reliable chain IDs for most rows, so verified residue-to-structure mapping was not bulk executed.


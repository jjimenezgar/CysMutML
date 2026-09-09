"""Audit frozen mutation identities and candidate source matches; never train."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--artifact', required=True)
    parser.add_argument('--source-tables', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.artifact) as archive:
        df = pd.read_csv(archive.open(
            'data/processed/fireprotdb_mutations_aggregated.csv'), low_memory=False)
    missing = df[df.canonical_sequence.isna()].copy()
    missing['category'] = 'construct_or_design_name'
    mask = missing.protein_id.str.fullmatch(r'[0-9][A-Za-z0-9]{3}', na=False)
    missing.loc[mask, 'category'] = 'pdb_like_name'
    missing.loc[missing.uniprot_id.notna(), 'category'] = 'uniprot_present'
    counts = missing.groupby('category').agg(
        rows=('mutation', 'size'), proteins=('protein_id', 'nunique'),
        cys=('mut_aa', lambda x: x.eq('C').sum()))
    counts.to_csv(out / 'coverage.csv')
    sites = missing.groupby(['protein_id', 'position']).wt_aa.nunique()
    conflicts = sites[sites > 1].rename('wt_identities').reset_index()
    conflicts.to_csv(out / 'conflicting_sites.csv', index=False)
    affected = missing.merge(conflicts, on=['protein_id', 'position'])
    with zipfile.ZipFile(args.source_tables) as archive:
        source = pd.read_csv(archive.open(
            'Data_tables_for_figs/dG_site_feature_Fig3.csv'))
    source['protein_id'] = source.pdb_name.str.strip().str.replace(r'\.pdb$', '', regex=True)
    overlap = missing[missing.protein_id.isin(source.protein_id)]
    values = source.melt(id_vars=['protein_id', 'pos', 'wt_aa'],
                         value_vars=[c for c in source if c.startswith('ddg_')],
                         var_name='mut_aa', value_name='source_ddg')
    values.mut_aa = values.mut_aa.str[-1]
    values = values.rename(columns={'pos': 'position'})
    keys = ['protein_id', 'position', 'wt_aa', 'mut_aa']
    if values.duplicated(keys).any():
        raise ValueError('Ambiguous source mutation keys')
    joined = missing.merge(values, on=keys, validate='many_to_one')
    joined = joined[np.isfinite(joined.source_ddg)]
    same = np.isclose(joined.median_destabilization_ddg, joined.source_ddg,
                      atol=1e-5, rtol=1e-5)
    opposite = np.isclose(joined.median_destabilization_ddg, -joined.source_ddg,
                          atol=1e-5, rtol=1e-5)
    missing.groupby('protein_id').agg(
        rows=('mutation', 'size'), cys=('mut_aa', lambda x: x.eq('C').sum())
    ).assign(source_name_match=lambda x: x.index.isin(source.protein_id)).to_csv(
        out / 'protein_inventory.csv')
    report = dict(
        rows=len(df), missing_sequence_rows=len(missing),
        missing_sequence_cys_rows=int(missing.mut_aa.eq('C').sum()),
        missing_sequence_protein_names=missing.protein_id.nunique(),
        missing_sequence_without_uniprot=int(missing.uniprot_id.isna().sum()),
        categories=counts.reset_index().to_dict('records'),
        conflicting_sites=len(conflicts), conflicting_protein_names=conflicts.protein_id.nunique(),
        conflicting_site_rows=len(affected),
        conflicting_site_cys_rows=int(affected.mut_aa.eq('C').sum()),
        source_name_matches=overlap.protein_id.nunique(),
        source_name_matched_rows=len(overlap),
        source_name_matched_cys_rows=int(overlap.mut_aa.eq('C').sum()),
        exact_source_key_rows=len(joined),
        exact_source_key_cys_rows=int(joined.mut_aa.eq('C').sum()),
        same_value_rows=int(same.sum()), opposite_value_rows=int(opposite.sum()),
        new_training_ready_sequences=0, model_changed=False,
        source='https://zenodo.org/records/7992926',
        source_member='Data_tables_for_figs/dG_site_feature_Fig3.csv',
        caveat='Name/site matches are candidates, not verified full construct sequences. '
               'Source keys use supplied positions without offsets. Reconcile construct '
               'identity, aggregation, assay quality and sign before training.',
        input_sha256={Path(p).name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                      for p in [args.artifact, args.source_tables]},
    )
    (out / 'audit.json').write_text(json.dumps(report, indent=2, default=int) + '\n')
    print(json.dumps(report, indent=2, default=int))


if __name__ == '__main__':
    main()

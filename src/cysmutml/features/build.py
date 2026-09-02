"""Feature matrix construction."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pandas.api.types import is_numeric_dtype

from cysmutml.amino_acids import physicochemical_features
from cysmutml.data.fireprotdb import normalize_fireprotdb_table

STRUCTURAL_COLUMNS = [
    "abs_sasa",
    "relative_sasa",
    "residue_sasa_abs",
    "residue_sasa_rel",
    "ca_neighbors_6a",
    "ca_neighbors_8a",
    "ca_neighbors_10a",
    "ca_contacts_6A",
    "ca_contacts_8A",
    "ca_contacts_10A",
    "heavy_atom_contact_count",
    "normalized_ca_distance_to_center",
    "normalized_ca_distance_to_structure_center",
    "local_density_10a",
    "mean_b_factor",
    "normalized_b_factor",
    "mean_residue_b_factor",
    "chain_normalized_b_factor",
    "secondary_structure",
]


def add_physicochemical_features(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        features = physicochemical_features(str(row["wt_aa"]), str(row["mut_aa"]))
        rows.append({**row.to_dict(), **features})
    out = pd.DataFrame(rows)
    for column in out.columns:
        if column not in {
            "wt_aa",
            "mut_aa",
            "protein_id",
            "pdb_id",
            "chain",
            "canonical_sequence",
            "mutation",
            "source_ddg_units",
            "source_ddg_sign_convention",
            "secondary_structure",
        }:
            try:
                out[column] = pd.to_numeric(out[column], errors="raise")
            except (TypeError, ValueError):
                pass
    return out


def build_feature_table(input_csv: str | Path, output_csv: str | Path) -> pd.DataFrame:
    df = pd.read_csv(input_csv, low_memory=False)
    if (
        "destabilization_ddg_kcal_mol" not in df.columns
        and "median_destabilization_ddg" in df.columns
    ):
        df["destabilization_ddg_kcal_mol"] = df["median_destabilization_ddg"]
    if not {"wt_aa", "mut_aa", "destabilization_ddg_kcal_mol"}.issubset(df.columns):
        df, _ = normalize_fireprotdb_table(df)
    featured = add_physicochemical_features(df)
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    featured.to_csv(output_csv, index=False)
    return featured


def feature_columns(
    df: pd.DataFrame, include_structural: bool = True
) -> tuple[list[str], list[str]]:
    excluded = {
        "destabilization_ddg_kcal_mol",
        "source_ddg_value",
        "source_ddg_units",
        "source_ddg_sign_convention",
        "canonical_sequence",
        "source_row_index",
        "uniprot_id",
        "fireprotdb_sequence_id",
        "pdb_id",
        "chain",
        "position",
        "protein_id",
        "sequence_cluster",
        "representative_protein_id",
        "mutation",
        "method",
        "measure",
        "ph",
        "exp_temperature",
        "source_dataset",
        "publication_pmid",
        "publication_doi",
        "publication_year",
        "row_id",
        "execution_subset_max_rows",
        "dataset_chain",
        "selected_chain",
        "original_mutation",
        "dataset_position",
        "dataset_wt",
        "dataset_mut",
        "mapped_pdb_resseq",
        "insertion_code",
        "mapped_wt",
        "sequence_identity",
        "alignment_coverage",
        "mapping_status",
        "failure_reason",
        "median_destabilization_ddg",
        "mean_destabilization_ddg",
        "std_destabilization_ddg",
        "min_destabilization_ddg",
        "max_destabilization_ddg",
        "n_measurements",
        "pdb_id_values",
    }
    structural = set(STRUCTURAL_COLUMNS)
    cols = [col for col in df.columns if col not in excluded]
    if not include_structural:
        cols = [col for col in cols if col not in structural]
    categorical = [col for col in cols if not is_numeric_dtype(df[col])]
    numeric = [col for col in cols if col not in categorical]
    return numeric, categorical

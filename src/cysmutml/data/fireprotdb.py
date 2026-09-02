"""FireProtDB acquisition and normalization utilities."""

from __future__ import annotations

import json
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import urlopen

import pandas as pd

from cysmutml.mutations import parse_mutation

FIREPROTDB_SEARCH_URL = "https://loschmidt.chemi.muni.cz/fireprotdb/api/search"
FIREPROTDB_SEQUENCE_URL = "https://loschmidt.chemi.muni.cz/fireprotdb/api/sequences"


def download_fireprotdb_csv(output_path: str | Path) -> Path:
    """Download mutant entries with non-empty DDG through the documented API."""
    query = {
        "tree": {
            "operator": "AND",
            "children": [
                {"variable": "MUTATED_POSITION", "operator": "IS_EMPTY", "value": False},
                {"variable": "DDG", "operator": "IS_EMPTY", "value": False},
            ],
        }
    }
    params = urlencode({"query": json.dumps(query), "format": "csv", "inline": "true"})
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(f"{FIREPROTDB_SEARCH_URL}?{params}") as response:
        output_path.write_bytes(response.read())
    return output_path


def download_fireprotdb_sequences(
    sequence_ids: list[str],
) -> tuple[dict[str, str], list[str]]:
    """Fetch canonical sequences from FireProtDB's documented sequence endpoint."""
    sequences: dict[str, str] = {}
    failed: list[str] = []
    for sequence_id in sorted(set(sequence_ids)):
        try:
            url = f"{FIREPROTDB_SEQUENCE_URL}/{quote(sequence_id, safe='')}/sequence"
            with urlopen(url, timeout=30) as response:
                payload = json.loads(response.read().decode("utf-8"))
            sequence = payload if isinstance(payload, str) else payload.get("sequence")
            normalized = "".join(str(sequence).split()).upper()
            if not normalized or normalized in {"NAN", "NONE"}:
                raise ValueError("empty sequence")
            sequences[sequence_id] = normalized
        except (OSError, TypeError, ValueError):
            failed.append(sequence_id)
    if not sequences:
        raise RuntimeError("FireProtDB returned no canonical sequences")
    return sequences, failed


def _identifier(value: object) -> str | None:
    if pd.isna(value):
        return None
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    return text or None


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    normalized = {column.lower().strip().replace(" ", "_"): column for column in df.columns}
    for candidate in candidates:
        key = candidate.lower().strip().replace(" ", "_")
        if key in normalized:
            return normalized[key]
    return None


def normalize_fireprotdb_table(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int | None]]:
    """Normalize a FireProtDB-like CSV/export into the internal schema."""
    raw_rows = len(df)
    mutation_col = find_column(
        df, ["mutation", "mutations", "mutant", "mutation_code", "substitution"]
    )
    ddg_col = find_column(df, ["ddg", "ΔΔg", "delta_delta_g", "ddg_(kcal/mol)", "ddg_kcal_mol"])
    protein_col = find_column(
        df, ["protein", "protein_name", "source_sequence_id", "uniprotkb", "uniprot"]
    )
    pdb_col = find_column(df, ["pdb", "pdb_id", "structure", "wwpdb"])
    chain_col = find_column(df, ["chain", "pdb_chain"])
    sequence_col = find_column(df, ["sequence", "protein_sequence", "source_sequence"])
    source_sequence_id_col = find_column(df, ["source_sequence_id"])
    method_col = find_column(df, ["method"])
    measure_col = find_column(df, ["measure"])
    ph_col = find_column(df, ["ph"])
    temp_col = find_column(df, ["exp_temperature"])
    source_dataset_col = find_column(df, ["source_dataset"])
    uniprot_col = find_column(df, ["uniprotkb", "uniprot"])
    pmid_col = find_column(df, ["publication_pmid", "pmid"])
    doi_col = find_column(df, ["publication_doi", "doi"])
    year_col = find_column(df, ["publication_year", "year"])

    if mutation_col is None or ddg_col is None:
        raise ValueError("Could not identify mutation and DDG columns in FireProtDB export")

    rows = []
    rejected = 0
    for source_row_index, row in df.iterrows():
        try:
            mutation = parse_mutation(str(row[mutation_col]).split(",")[0])
            ddg = pd.to_numeric(row[ddg_col], errors="raise")
        except Exception:
            rejected += 1
            continue
        rows.append(
            {
                "protein_id": str(row[protein_col]) if protein_col else "unknown",
                "source_row_index": int(source_row_index),
                "pdb_id": str(row[pdb_col]).lower() if pdb_col and pd.notna(row[pdb_col]) else None,
                "uniprot_id": str(row[uniprot_col])
                if uniprot_col and pd.notna(row[uniprot_col])
                else None,
                "chain": str(row[chain_col]) if chain_col and pd.notna(row[chain_col]) else None,
                "canonical_sequence": str(row[sequence_col])
                if sequence_col and pd.notna(row[sequence_col])
                else None,
                "fireprotdb_sequence_id": _identifier(row[source_sequence_id_col])
                if source_sequence_id_col
                else None,
                "mutation": mutation.label,
                "wt_aa": mutation.wt,
                "position": mutation.position,
                "mut_aa": mutation.mut,
                "source_ddg_value": float(ddg),
                "source_ddg_units": "kcal/mol",
                "source_ddg_sign_convention": (
                    "FireProtDB: negative stabilizing, positive destabilizing"
                ),
                "destabilization_ddg_kcal_mol": float(ddg),
                "method": str(row[method_col])
                if method_col and pd.notna(row[method_col])
                else None,
                "measure": str(row[measure_col])
                if measure_col and pd.notna(row[measure_col])
                else None,
                "ph": row[ph_col] if ph_col and pd.notna(row[ph_col]) else None,
                "exp_temperature": row[temp_col] if temp_col and pd.notna(row[temp_col]) else None,
                "source_dataset": str(row[source_dataset_col])
                if source_dataset_col and pd.notna(row[source_dataset_col])
                else None,
                "publication_pmid": str(row[pmid_col])
                if pmid_col and pd.notna(row[pmid_col])
                else None,
                "publication_doi": str(row[doi_col])
                if doi_col and pd.notna(row[doi_col])
                else None,
                "publication_year": row[year_col] if year_col and pd.notna(row[year_col]) else None,
            }
        )
    out = pd.DataFrame(rows)
    summary = {
        "raw_records": raw_rows,
        "single_mutations_with_valid_target": len(out),
        "rejected_records": rejected,
    }
    return out, summary


def prepare_data(
    raw_path: str | Path,
    output_path: str | Path,
    fetch_sequences: bool = False,
) -> dict[str, int | None]:
    df = pd.read_csv(raw_path, low_memory=False)
    normalized, summary = normalize_fireprotdb_table(df)
    if fetch_sequences:
        sequence_ids = normalized["fireprotdb_sequence_id"].dropna().astype(str).unique().tolist()
        if not sequence_ids:
            raise ValueError("FireProtDB export has no SOURCE_SEQUENCE_ID values")
        sequences, failed = download_fireprotdb_sequences(sequence_ids)
        fetched = normalized["fireprotdb_sequence_id"].map(sequences)
        normalized["canonical_sequence"] = normalized["canonical_sequence"].fillna(fetched)
        summary["canonical_sequences_downloaded"] = len(sequences)
        summary["canonical_sequence_download_failures"] = len(failed)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(output_path, index=False)
    return summary

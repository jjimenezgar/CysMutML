"""Structure and sequence acquisition for structural mapping."""

from __future__ import annotations

from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

import pandas as pd

from cysmutml.structures.io import parse_pdb


def _valid_pdb(path: Path) -> tuple[bool, str | None]:
    try:
        structure = parse_pdb(path)
        has_atom = any(True for _ in structure.get_atoms())
    except Exception as exc:
        return False, str(exc)
    return (True, None) if has_atom else (False, "no_atoms")


def download_pdbs_for_candidates(
    candidate_csv: str | Path = "reports/structure_candidate_records.csv",
    structures_dir: str | Path = "data/structures",
    report_csv: str | Path = "reports/structure_download_report.csv",
) -> pd.DataFrame:
    candidates = pd.read_csv(candidate_csv, low_memory=False)
    pdb_ids = sorted(candidates["pdb_id"].dropna().astype(str).str.lower().unique())
    structures_dir = Path(structures_dir)
    structures_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb_id in pdb_ids:
        path = structures_dir / f"{pdb_id}.pdb"
        if path.exists():
            valid, reason = _valid_pdb(path)
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "download_status": "cached_valid" if valid else "cached_invalid",
                    "local_path": str(path) if valid else None,
                    "failure_reason": reason,
                }
            )
            continue
        try:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            with urlopen(url, timeout=30) as resp:
                path.write_bytes(resp.read())
            valid, reason = _valid_pdb(path)
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "download_status": "downloaded_valid" if valid else "downloaded_invalid",
                    "local_path": str(path) if valid else None,
                    "failure_reason": reason,
                }
            )
        except (OSError, URLError) as exc:
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "download_status": "failed",
                    "local_path": None,
                    "failure_reason": str(exc),
                }
            )
    report = pd.DataFrame(rows)
    Path(report_csv).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(report_csv, index=False)
    return report


def download_uniprot_sequences_for_candidates(
    candidate_csv: str | Path = "reports/structure_candidate_records.csv",
    sequence_dir: str | Path = "data/sequences",
    report_csv: str | Path = "reports/sequence_download_report.csv",
) -> pd.DataFrame:
    candidates = pd.read_csv(candidate_csv, low_memory=False)
    accessions = sorted(candidates["uniprot_id"].dropna().astype(str).unique())
    sequence_dir = Path(sequence_dir)
    sequence_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for accession in accessions:
        path = sequence_dir / f"{accession}.fasta"
        if path.exists() and path.stat().st_size > 0:
            rows.append(
                {
                    "uniprot_id": accession,
                    "download_status": "cached",
                    "local_path": str(path),
                    "failure_reason": None,
                }
            )
            continue
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
            with urlopen(url, timeout=30) as resp:
                payload = resp.read()
            if not payload.startswith(b">"):
                raise ValueError("response_not_fasta")
            path.write_bytes(payload)
            rows.append(
                {
                    "uniprot_id": accession,
                    "download_status": "downloaded",
                    "local_path": str(path),
                    "failure_reason": None,
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "uniprot_id": accession,
                    "download_status": "failed",
                    "local_path": None,
                    "failure_reason": str(exc),
                }
            )
    report = pd.DataFrame(rows)
    Path(report_csv).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(report_csv, index=False)
    return report


def read_fasta_sequence(path: str | Path) -> str:
    lines = [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ]
    return "".join(lines)

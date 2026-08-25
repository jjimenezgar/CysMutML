"""Simple ranking visualization outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio.PDB import PDBIO

from cysmutml.structures.features import residue_key
from cysmutml.structures.io import parse_pdb


def write_ranked_bfactor_pdb(
    ranking_csv: str | Path,
    pdb_path: str | Path,
    output_pdb: str | Path,
    score_column: str = "cys_suitability_score",
) -> Path:
    """Write a PDB copy with candidate ranking score encoded in B-factor."""
    df = pd.read_csv(ranking_csv)
    score_lookup = {
        (str(row["chain"]), str(row["residue_number"])): float(row[score_column]) * 100.0
        for _, row in df.iterrows()
    }
    structure = parse_pdb(pdb_path)
    for model in structure:
        for chain in model:
            for residue in chain:
                score = score_lookup.get((chain.id, residue_key(residue)), 0.0)
                for atom in residue:
                    atom.bfactor = score
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb))
    return output_pdb


def write_pymol_script(
    ranking_csv: str | Path, pdb_path: str | Path, output_pml: str | Path
) -> Path:
    df = pd.read_csv(ranking_csv)
    output_pml = Path(output_pml)
    output_pml.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"load {Path(pdb_path).resolve()}, ranked_structure",
        "hide everything",
        "show cartoon",
        "spectrum b, blue_white_red, ranked_structure",
    ]
    top = df.head(10)
    for _, row in top.iterrows():
        rank = int(row["rank_engineering"])
        lines.append(
            f"select candidate_{rank}, chain {row['chain']} "
            f"and resi {row['residue_number']}"
        )
        lines.append(f"show sticks, candidate_{rank}")
    output_pml.write_text("\n".join(lines) + "\n")
    return output_pml

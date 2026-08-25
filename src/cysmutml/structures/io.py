"""PDB parsing helpers."""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlretrieve

from Bio.PDB import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


def parse_pdb(path: str | Path) -> Structure:
    parser = PDBParser(QUIET=True)
    return parser.get_structure(Path(path).stem, str(path))


def first_model(structure: Structure) -> Model:
    return next(structure.get_models())


def get_chain(structure: Structure, chain_id: str) -> Chain:
    model = first_model(structure)
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id!r} not found in structure")
    return model[chain_id]


def download_pdb(pdb_id: str, output_dir: str | Path = "data/structures") -> Path:
    pdb_id = pdb_id.lower()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{pdb_id}.pdb"
    if path.exists():
        return path
    urlretrieve(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", path)
    return path

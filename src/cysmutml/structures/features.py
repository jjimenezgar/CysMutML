"""Interpretable structural features around a residue."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from Bio.PDB import ShrakeRupley
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from cysmutml.amino_acids import MAX_ASA_TIEN_2013, THREE_TO_ONE
from cysmutml.structures.io import get_chain, parse_pdb


def is_standard_residue(residue: Residue) -> bool:
    return residue.id[0] == " " and residue.resname in THREE_TO_ONE


def residue_one_letter(residue: Residue) -> str:
    return THREE_TO_ONE[residue.resname]


def chain_residues(chain: Chain) -> list[Residue]:
    return [res for res in chain if is_standard_residue(res)]


def ca_coord(residue: Residue) -> np.ndarray | None:
    if "CA" not in residue:
        return None
    return np.asarray(residue["CA"].coord, dtype=float)


def residue_key(residue: Residue) -> str:
    het, number, insertion = residue.id
    return f"{number}{insertion.strip()}" if het == " " else f"{het}:{number}{insertion.strip()}"


def compute_sasa_by_residue(pdb_path: str | Path) -> dict[tuple[str, str], float]:
    structure = parse_pdb(pdb_path)
    ShrakeRupley().compute(structure, level="R")
    out: dict[tuple[str, str], float] = {}
    for chain in next(structure.get_models()):
        for residue in chain_residues(chain):
            out[(chain.id, residue_key(residue))] = float(getattr(residue, "sasa", math.nan))
    return out


def residue_structural_features(
    pdb_path: str | Path,
    chain_id: str,
    residue: Residue,
    contact_radii: tuple[float, ...] = (6.0, 8.0, 10.0),
    heavy_atom_contact_radius: float = 4.5,
) -> dict[str, float | int | str]:
    structure = parse_pdb(pdb_path)
    chain = get_chain(structure, chain_id)
    residues = chain_residues(chain)
    residue_id = residue.id
    target = next((res for res in residues if res.id == residue_id), None)
    if target is None:
        raise ValueError("Residue not found in reparsed chain")

    ShrakeRupley().compute(structure, level="R")
    aa = residue_one_letter(target)
    target_ca = ca_coord(target)
    if target_ca is None:
        raise ValueError(f"Residue {residue_key(target)} has no CA atom")

    ca_points = [
        ca_coord(res) for res in residues if res.id != target.id and ca_coord(res) is not None
    ]
    ca_points_arr = np.asarray(ca_points, dtype=float)
    distances = (
        np.linalg.norm(ca_points_arr - target_ca, axis=1) if len(ca_points_arr) else np.array([])
    )

    all_ca = np.asarray(
        [ca_coord(res) for res in residues if ca_coord(res) is not None], dtype=float
    )
    center = all_ca.mean(axis=0)
    max_center_distance = max(float(np.linalg.norm(point - center)) for point in all_ca) or 1.0

    target_atoms = [atom for atom in target.get_atoms() if atom.element != "H"]
    other_atoms = [
        atom
        for res in residues
        if res.id != target.id
        for atom in res.get_atoms()
        if atom.element != "H"
    ]
    heavy_contacts = 0
    for atom in target_atoms:
        coord = np.asarray(atom.coord, dtype=float)
        if any(
            float(np.linalg.norm(coord - np.asarray(other.coord))) <= heavy_atom_contact_radius
            for other in other_atoms
        ):
            heavy_contacts += 1

    b_factors = [float(atom.bfactor) for atom in target.get_atoms()]
    chain_b = [float(atom.bfactor) for res in residues for atom in res.get_atoms()]
    chain_b_mean = float(np.mean(chain_b)) if chain_b else math.nan
    chain_b_std = float(np.std(chain_b)) if chain_b else math.nan
    residue_b = float(np.mean(b_factors)) if b_factors else math.nan

    sasa = float(getattr(target, "sasa", math.nan))
    max_asa = MAX_ASA_TIEN_2013[aa]
    features: dict[str, float | int | str] = {
        "abs_sasa": sasa,
        "relative_sasa": min(sasa / max_asa, 1.5) if not math.isnan(sasa) else math.nan,
        "normalized_ca_distance_to_center": float(
            np.linalg.norm(target_ca - center) / max_center_distance
        ),
        "local_density_10a": int(np.sum(distances <= 10.0)) if len(distances) else 0,
        "heavy_atom_contact_count": heavy_contacts,
        "mean_b_factor": residue_b,
        "normalized_b_factor": (residue_b - chain_b_mean) / chain_b_std if chain_b_std > 0 else 0.0,
        "secondary_structure": "unknown",
    }
    for radius in contact_radii:
        features[f"ca_neighbors_{int(radius)}a"] = (
            int(np.sum(distances <= radius)) if len(distances) else 0
        )
    return features


def chain_feature_rows(pdb_path: str | Path, chain_id: str) -> list[dict[str, float | int | str]]:
    structure = parse_pdb(pdb_path)
    chain = get_chain(structure, chain_id)
    rows = []
    for residue in chain_residues(chain):
        features = residue_structural_features(pdb_path, chain_id, residue)
        features.update(
            {
                "chain": chain_id,
                "residue_number": residue_key(residue),
                "wt_aa": residue_one_letter(residue),
            }
        )
        rows.append(features)
    return rows

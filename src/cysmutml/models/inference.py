"""Prediction on new structures."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from cysmutml.amino_acids import MAX_ASA_TIEN_2013, physicochemical_features
from cysmutml.config import load_config
from cysmutml.ranking.engineering import stability_component_from_ddg
from cysmutml.structures.features import (
    ca_coord,
    chain_residues,
    compute_sasa_by_residue,
    residue_key,
    residue_one_letter,
)
from cysmutml.structures.io import get_chain, parse_pdb


def _secondary_structure_map(pdb_path: str | Path) -> dict[tuple[str, str], str]:
    """Return DSSP secondary-structure codes keyed by (chain ID, residue number)."""
    try:
        import mdtraj as md

        trajectory = md.load(str(pdb_path))
        dssp = trajectory.compute_dssp()[0]
        structure = parse_pdb(pdb_path)
        chain_ids = [chain.id for chain in next(structure.get_models())]
        output: dict[tuple[str, str], str] = {}
        for residue, code in zip(trajectory.topology.residues, dssp):
            chain_index = residue.chain.index
            chain_id = chain_ids[chain_index] if chain_index < len(chain_ids) else str(chain_index)
            output[(chain_id, str(residue.resSeq))] = str(code)
        return output
    except (ImportError, OSError, ValueError, IndexError, RuntimeError):
        # DSSP is an informative structural signal; an unavailable assignment
        # should not prevent the rest of the ranking from running.
        return {}


def _secondary_structure_penalty(code: str) -> float:
    """Penalize residues assigned to helices or beta sheets, not loops."""
    return 1.0 if code in {"H", "B", "E", "G", "I"} else 0.0


def _parse_protected_residues(protected_residues: str | None) -> set[tuple[str, str]]:
    if not protected_residues:
        return set()
    parsed = set()
    for item in protected_residues.split(","):
        if not item.strip():
            continue
        try:
            chain, residue = item.strip().split(":", maxsplit=1)
        except ValueError as exc:
            raise ValueError(
                "Protected residues must use comma-separated CHAIN:RESIDUE entries, "
                "for example A:45,A:48,B:120."
            ) from exc
        if not chain or not residue:
            raise ValueError(
                "Protected residues must use comma-separated CHAIN:RESIDUE entries, "
                "for example A:45,A:48,B:120."
            )
        parsed.add((chain, residue))
    return parsed


def _relative_sasa_for_residue(
    sasa_by_residue: dict[tuple[str, str], float], chain_id: str, residue
) -> float:
    aa = residue_one_letter(residue)
    sasa = sasa_by_residue.get((chain_id, residue_key(residue)), np.nan)
    if np.isnan(sasa):
        return np.nan
    return min(float(sasa) / MAX_ASA_TIEN_2013[aa], 1.5)


def _chain_structural_feature_map(
    residues,
    chain_id: str,
    sasa_by_residue: dict[tuple[str, str], float],
    secondary_structure_by_residue: dict[tuple[str, str], str] | None = None,
) -> dict[str, dict[str, float | int | str]]:
    ca_coords = {residue_key(res): ca_coord(res) for res in residues if ca_coord(res) is not None}
    all_ca = np.asarray(list(ca_coords.values()), dtype=float)
    center = all_ca.mean(axis=0)
    max_center_distance = max(float(np.linalg.norm(point - center)) for point in all_ca) or 1.0
    chain_b = [float(atom.bfactor) for res in residues for atom in res.get_atoms()]
    chain_b_mean = float(np.mean(chain_b)) if chain_b else np.nan
    chain_b_std = float(np.std(chain_b)) if chain_b else np.nan
    out = {}
    for residue in residues:
        key = residue_key(residue)
        aa = residue_one_letter(residue)
        target_ca = ca_coords.get(key)
        if target_ca is None:
            continue
        other_ca = np.asarray(
            [coord for other_key, coord in ca_coords.items() if other_key != key], dtype=float
        )
        distances = (
            np.linalg.norm(other_ca - target_ca, axis=1) if len(other_ca) else np.array([])
        )
        b_factors = [float(atom.bfactor) for atom in residue.get_atoms()]
        residue_b = float(np.mean(b_factors)) if b_factors else np.nan
        sasa = float(sasa_by_residue.get((chain_id, key), np.nan))
        out[key] = {
            "abs_sasa": sasa,
            "relative_sasa": min(sasa / MAX_ASA_TIEN_2013[aa], 1.5)
            if not np.isnan(sasa)
            else np.nan,
            "normalized_ca_distance_to_center": float(
                np.linalg.norm(target_ca - center) / max_center_distance
            ),
            "local_density_10a": int(np.sum(distances <= 10.0)) if len(distances) else 0,
            "heavy_atom_contact_count": np.nan,
            "mean_b_factor": residue_b,
            "normalized_b_factor": (residue_b - chain_b_mean) / chain_b_std
            if chain_b_std > 0
            else 0.0,
            "secondary_structure": (secondary_structure_by_residue or {}).get(
                (chain_id, key), "unknown"
            ),
            "ca_neighbors_6a": int(np.sum(distances <= 6.0)) if len(distances) else 0,
            "ca_neighbors_8a": int(np.sum(distances <= 8.0)) if len(distances) else 0,
            "ca_neighbors_10a": int(np.sum(distances <= 10.0)) if len(distances) else 0,
        }
    return out


def generate_cys_feature_rows(
    pdb_path: str | Path,
    chain_id: str,
    protected_residues: str | None = None,
    config_path: str | Path = "configs/default.yaml",
    monocysteine_design: bool = False,
) -> pd.DataFrame:
    config = load_config(config_path)
    lys_config = config.get("lysine_environment", {})
    existing_cys_config = config.get("existing_cys", {})
    monocys_config = config.get("monocysteine_design", {})
    lys_radius = float(lys_config.get("radius_angstrom", 20.0))
    lys_exposure_threshold = float(lys_config.get("exposure_relative_sasa_threshold", 0.25))
    cys_count_radii = [
        float(radius) for radius in existing_cys_config.get("count_radii_angstrom", [6, 8, 10, 15])
    ]
    exposed_native_cys_threshold = float(
        monocys_config.get("exposed_cys_relative_sasa_threshold", 0.25)
    )
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    structure = parse_pdb(pdb_path)
    chain = get_chain(structure, chain_id)
    residues = chain_residues(chain)
    all_residues = [
        (other_chain.id, residue)
        for model in structure
        for other_chain in model
        for residue in chain_residues(other_chain)
    ]
    sasa_by_residue = compute_sasa_by_residue(pdb_path)
    secondary_structure_by_residue = _secondary_structure_map(pdb_path)
    structural_by_residue = _chain_structural_feature_map(
        residues, chain_id, sasa_by_residue, secondary_structure_by_residue
    )
    existing_cys_entries = [
        (other_chain_id, residue, ca_coord(residue))
        for other_chain_id, residue in all_residues
        if residue_one_letter(residue) == "C"
    ]
    existing_cys_entries = [
        (other_chain_id, residue, coord)
        for other_chain_id, residue, coord in existing_cys_entries
        if coord is not None
    ]
    exposed_native_cys_count = sum(
        _relative_sasa_for_residue(sasa_by_residue, other_chain_id, residue)
        >= exposed_native_cys_threshold
        for other_chain_id, residue, _ in existing_cys_entries
    )
    protected = _parse_protected_residues(protected_residues)
    protected_coords = [
        ca_coord(residue)
        for other_chain_id, residue in all_residues
        if (other_chain_id, residue_key(residue)) in protected and ca_coord(residue) is not None
    ]
    rows = []
    for residue in residues:
        wt = residue_one_letter(residue)
        if wt == "C":
            continue
        physchem = physicochemical_features(wt, "C")
        number = residue_key(residue)
        structural = structural_by_residue[number]
        target_ca = ca_coord(residue)
        lys_distances = []
        exposed_lys_distances = []
        cys_distances = []
        if target_ca is not None:
            for other_chain_id, other_residue in all_residues:
                if other_chain_id == chain_id and other_residue.id == residue.id:
                    continue
                other_ca = ca_coord(other_residue)
                if other_ca is None:
                    continue
                distance = float(np.linalg.norm(target_ca - other_ca))
                if residue_one_letter(other_residue) == "K" and distance <= lys_radius:
                    lys_distances.append(distance)
                    lys_rel_sasa = _relative_sasa_for_residue(
                        sasa_by_residue, other_chain_id, other_residue
                    )
                    if lys_rel_sasa >= lys_exposure_threshold:
                        exposed_lys_distances.append(distance)
            for _, _cys_residue, cys_coord in existing_cys_entries:
                cys_distances.append(float(np.linalg.norm(target_ca - cys_coord)))
        if cys_distances:
            distance_to_existing_cys = min(cys_distances)
        else:
            distance_to_existing_cys = np.nan
        if target_ca is not None and protected_coords:
            distance_to_nearest_protected = min(
                float(np.linalg.norm(target_ca - coord)) for coord in protected_coords
            )
        else:
            distance_to_nearest_protected = np.nan
        rows.append(
            {
                "chain": chain_id,
                "residue_number": number,
                "position": int(residue.id[1]),
                "wt_aa": wt,
                "mutation": f"{wt}{number}C",
                "distance_to_existing_cys": distance_to_existing_cys,
                "nearest_existing_cys_distance": distance_to_existing_cys,
                "local_lys_count": len(lys_distances),
                "local_exposed_lys_count": len(exposed_lys_distances),
                "lysine_radius_angstrom": lys_radius,
                "lys_exposure_relative_sasa_threshold": lys_exposure_threshold,
                "native_cys_count_total": len(existing_cys_entries),
                "native_exposed_cys_count": exposed_native_cys_count,
                "native_cys_context": "monocysteine_caution"
                if monocysteine_design and exposed_native_cys_count
                else "native_cys_present"
                if existing_cys_entries
                else "no_native_cys_detected",
                "distance_to_nearest_protected": distance_to_nearest_protected,
                **physchem,
                **structural,
                "secondary_structure_penalty": _secondary_structure_penalty(
                    str(structural.get("secondary_structure", "unknown"))
                ),
                **{
                    f"existing_cys_count_{int(radius)}A": int(
                        sum(distance <= radius for distance in cys_distances)
                    )
                    for radius in cys_count_radii
                },
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(
            f"No eligible non-cysteine standard residues found in chain {chain_id}; "
            "cannot generate X->Cys candidates."
        )
    return out


def out_of_domain_warnings(
    features: pd.DataFrame, artifact: dict, margin_fraction: float = 0.05
) -> list[str]:
    warnings = []
    for column, limits in artifact.get("training_feature_ranges", {}).items():
        if column not in features:
            continue
        low, high = limits["min"], limits["max"]
        margin = (high - low) * margin_fraction
        outside = int(
            ((features[column] < low - margin) | (features[column] > high + margin)).sum()
        )
        if outside:
            warnings.append(
                f"{column}: {outside} values outside training range [{low:.3g}, {high:.3g}]"
            )
    return warnings


def predict_cys_mutations(
    pdb_path: str | Path,
    chain_id: str,
    model_path: str | Path,
    output_dir: str | Path,
    protected_residues: str | None = None,
    config_path: str | Path = "configs/default.yaml",
    monocysteine_design: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    config = load_config(config_path)
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    artifact = joblib.load(model_path)
    features = generate_cys_feature_rows(
        pdb_path,
        chain_id,
        protected_residues,
        config_path=config_path,
        monocysteine_design=monocysteine_design,
    )
    model_features = artifact["numeric_features"] + artifact["categorical_features"]
    for col in artifact["numeric_features"]:
        if col not in features:
            features[col] = np.nan
        else:
            features[col] = pd.to_numeric(features[col], errors="coerce")
    for col in artifact["categorical_features"]:
        if col not in features:
            features[col] = "missing"
    predictions = artifact["model"].predict(features[model_features])
    out = features.copy()
    out["predicted_destabilization_ddg"] = predictions
    ranking_config = config.get("ranking", {})
    out["stability_component"] = stability_component_from_ddg(
        out["predicted_destabilization_ddg"],
        favorable_ddg=float(ranking_config.get("stability_reference_ddg_low", -1.0)),
        unfavorable_ddg=float(ranking_config.get("stability_reference_ddg_high", 2.0)),
    )
    if "normalized_b_factor" in out:
        out["local_flexibility_proxy"] = pd.to_numeric(out["normalized_b_factor"], errors="coerce")
        out["flexibility_value"] = out["local_flexibility_proxy"]
        out["flexibility_method"] = "BFACTOR"
    else:
        out["local_flexibility_proxy"] = np.nan
        out["flexibility_value"] = np.nan
        out["flexibility_method"] = "UNAVAILABLE"
    out = out.sort_values("predicted_destabilization_ddg", ascending=True).reset_index(drop=True)
    out.insert(0, "rank_ml", range(1, len(out) + 1))
    warnings = out_of_domain_warnings(
        features,
        artifact,
        margin_fraction=float(
            config.get("out_of_domain", {}).get("numeric_margin_fraction", 0.05)
        ),
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_dir / "mutation_predictions.csv", index=False)
    if warnings:
        (output_dir / "out_of_domain_warnings.txt").write_text("\n".join(warnings) + "\n")
    return out, warnings

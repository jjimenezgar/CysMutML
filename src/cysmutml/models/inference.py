"""Prediction on new structures."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from cysmutml.amino_acids import physicochemical_features
from cysmutml.config import load_config
from cysmutml.ranking.engineering import stability_component_from_ddg
from cysmutml.structures.features import (
    ca_coord,
    chain_residues,
    residue_key,
    residue_one_letter,
    residue_structural_features,
)
from cysmutml.structures.io import get_chain, parse_pdb


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


def generate_cys_feature_rows(
    pdb_path: str | Path,
    chain_id: str,
    protected_residues: str | None = None,
) -> pd.DataFrame:
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
    existing_cys_coords = [
        ca_coord(residue) for _, residue in all_residues if residue_one_letter(residue) == "C"
    ]
    existing_cys_coords = [coord for coord in existing_cys_coords if coord is not None]
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
        structural = residue_structural_features(pdb_path, chain_id, residue)
        physchem = physicochemical_features(wt, "C")
        number = residue_key(residue)
        target_ca = ca_coord(residue)
        if target_ca is not None and existing_cys_coords:
            distance_to_existing_cys = min(
                float(np.linalg.norm(target_ca - coord)) for coord in existing_cys_coords
            )
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
                "distance_to_nearest_protected": distance_to_nearest_protected,
                **physchem,
                **structural,
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
) -> tuple[pd.DataFrame, list[str]]:
    config = load_config(config_path)
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    artifact = joblib.load(model_path)
    features = generate_cys_feature_rows(pdb_path, chain_id, protected_residues)
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
        out["flexibility_value"] = pd.to_numeric(out["normalized_b_factor"], errors="coerce")
        out["flexibility_method"] = "BFACTOR"
    else:
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

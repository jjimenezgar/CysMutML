"""Transparent cysteine engineering ranking heuristic."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cysmutml.config import load_config


def minmax_low_good(values: pd.Series) -> pd.Series:
    if values.max() == values.min():
        return pd.Series(1.0, index=values.index)
    return 1.0 - (values - values.min()) / (values.max() - values.min())


def minmax_high_good(values: pd.Series) -> pd.Series:
    if values.max() == values.min():
        return pd.Series(1.0, index=values.index)
    return (values - values.min()) / (values.max() - values.min())


def stability_component_from_ddg(
    values: pd.Series,
    favorable_ddg: float = -1.0,
    unfavorable_ddg: float = 2.0,
) -> pd.Series:
    """Map predicted DDG to [0, 1], where high means lower predicted destabilization."""
    if unfavorable_ddg <= favorable_ddg:
        raise ValueError("unfavorable_ddg must be greater than favorable_ddg")
    scaled = 1.0 - (values.astype(float) - favorable_ddg) / (unfavorable_ddg - favorable_ddg)
    return pd.Series(np.clip(scaled, 0.0, 1.0), index=values.index)


def accessibility_component_from_sasa(values: pd.Series) -> pd.Series:
    """Map relative SASA to [0, 1]; values >=1 are treated as fully exposed."""
    return pd.Series(np.clip(values.astype(float).fillna(0.0), 0.0, 1.0), index=values.index)


def rigidity_component_from_flexibility(values: pd.Series) -> pd.Series:
    """Map a target-protein flexibility proxy to [0, 1], with lower flexibility favored."""
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() < 2:
        return pd.Series(0.5, index=values.index)
    filled = numeric.fillna(numeric.median())
    return minmax_low_good(filled)


def rank_predictions(
    predictions_csv: str | Path,
    output_csv: str | Path,
    config_path: str | Path = "configs/default.yaml",
) -> pd.DataFrame:
    config = load_config(config_path)
    ranking_config = config.get("ranking", {})
    weights = ranking_config.get("weights", {})
    stability_weight = float(weights.get("stability", 0.50))
    accessibility_weight = float(weights.get("accessibility", weights.get("exposure", 0.30)))
    rigidity_weight = float(weights.get("rigidity", 0.20))
    protected_penalty_weight = float(weights.get("protected_penalty", 0.10))
    existing_cys_penalty_weight = float(weights.get("existing_cys_penalty", 0.05))
    existing_cys_warning_radius = float(
        ranking_config.get("existing_cys_warning_radius_angstrom", 6.0)
    )
    existing_cys_penalty_radius = float(
        ranking_config.get("existing_cys_penalty_radius_angstrom", existing_cys_warning_radius)
    )
    protected_radius = float(ranking_config.get("protected_site_radius_angstrom", 8.0))

    df = pd.read_csv(predictions_csv)
    if "stability_component" not in df:
        df["stability_component"] = stability_component_from_ddg(
            df["predicted_destabilization_ddg"].astype(float),
            favorable_ddg=float(ranking_config.get("stability_reference_ddg_low", -1.0)),
            unfavorable_ddg=float(ranking_config.get("stability_reference_ddg_high", 2.0)),
        )
    df["accessibility_component"] = accessibility_component_from_sasa(df["relative_sasa"])
    if "flexibility_value" not in df:
        if "normalized_b_factor" in df:
            df["flexibility_value"] = pd.to_numeric(df["normalized_b_factor"], errors="coerce")
            df["flexibility_method"] = "BFACTOR"
        else:
            df["flexibility_value"] = np.nan
            df["flexibility_method"] = "UNAVAILABLE"
    df["rigidity_component"] = rigidity_component_from_flexibility(df["flexibility_value"])

    if "distance_to_nearest_protected" in df:
        protected_distance = pd.to_numeric(df["distance_to_nearest_protected"], errors="coerce")
        df["protected_site_penalty"] = np.clip(
            (protected_radius - protected_distance.fillna(protected_radius)) / protected_radius,
            0.0,
            1.0,
        )
    else:
        df["protected_site_penalty"] = 0.0

    if "distance_to_existing_cys" in df:
        cys_distance = pd.to_numeric(df["distance_to_existing_cys"], errors="coerce")
        df["existing_cys_proximity_warning"] = (
            cys_distance <= existing_cys_warning_radius
        )
        df["existing_cys_penalty"] = np.clip(
            (existing_cys_penalty_radius - cys_distance.fillna(existing_cys_penalty_radius))
            / existing_cys_penalty_radius,
            0.0,
            1.0,
        )
    else:
        df["existing_cys_proximity_warning"] = False
        df["existing_cys_penalty"] = 0.0

    df["ranking_formula"] = (
        "stability_weight*stability_component + "
        "accessibility_weight*accessibility_component + "
        "rigidity_weight*rigidity_component - "
        "protected_penalty_weight*protected_site_penalty - "
        "existing_cys_penalty_weight*existing_cys_penalty"
    )
    df["cys_suitability_score"] = (
        stability_weight * df["stability_component"]
        + accessibility_weight * df["accessibility_component"]
        + rigidity_weight * df["rigidity_component"]
        - protected_penalty_weight * df["protected_site_penalty"]
        - existing_cys_penalty_weight * df["existing_cys_penalty"]
    )
    df = df.sort_values(
        ["cys_suitability_score", "predicted_destabilization_ddg"],
        ascending=[False, True],
    ).reset_index(drop=True)
    df.insert(0, "rank_engineering", range(1, len(df) + 1))
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    return df

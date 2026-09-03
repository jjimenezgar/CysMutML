"""Transparent cysteine engineering and rigidification ranking heuristics."""

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
    """Legacy transform: lower flexibility proxy receives higher rigidity score."""
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() < 2:
        return pd.Series(0.5, index=values.index)
    filled = numeric.fillna(numeric.median())
    return minmax_low_good(filled)


def flexibility_component_from_proxy(values: pd.Series) -> pd.Series:
    """Map a local flexibility proxy to [0, 1], with higher flexibility favored."""
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() < 2:
        return pd.Series(0.5, index=values.index)
    filled = numeric.fillna(numeric.median())
    return minmax_high_good(filled)


def lysine_environment_component_from_count(
    counts: pd.Series, saturation_k: float = 3.0
) -> pd.Series:
    """Saturating transform for exposed nearby Lys counts."""
    if saturation_k <= 0:
        raise ValueError("saturation_k must be positive")
    numeric = pd.to_numeric(counts, errors="coerce").fillna(0.0).clip(lower=0.0)
    return pd.Series(1.0 - np.exp(-numeric / saturation_k), index=counts.index)


def _series_or_default(df: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in df:
        return pd.Series(default, index=df.index)
    return pd.to_numeric(df[column], errors="coerce").fillna(default)


def rank_predictions(
    predictions_csv: str | Path,
    output_csv: str | Path,
    config_path: str | Path = "configs/default.yaml",
) -> pd.DataFrame:
    config = load_config(config_path)
    ranking_config = config.get("ranking", {})
    lys_config = config.get("lysine_environment", {})
    existing_cys_config = config.get("existing_cys", {})
    # Used only to report the optional protected-residue annotation.
    protected_radius = float(ranking_config.get("protected_site_radius_angstrom", 8.0))
    existing_cys_warning_radius = float(
        existing_cys_config.get("warning_radius_angstrom", 10.0)
    )
    existing_cys_penalty_radius = float(
        existing_cys_config.get("strong_penalty_radius_angstrom", 6.0)
    )

    df = pd.read_csv(predictions_csv)
    if "stability_component" not in df:
        df["stability_component"] = stability_component_from_ddg(
            df["predicted_destabilization_ddg"].astype(float),
            favorable_ddg=float(ranking_config.get("stability_reference_ddg_low", -1.0)),
            unfavorable_ddg=float(ranking_config.get("stability_reference_ddg_high", 2.0)),
        )

    df["accessibility_component"] = accessibility_component_from_sasa(df["relative_sasa"])

    if "local_flexibility_proxy" not in df:
        if "flexibility_value" in df:
            df["local_flexibility_proxy"] = pd.to_numeric(df["flexibility_value"], errors="coerce")
        elif "normalized_b_factor" in df:
            df["local_flexibility_proxy"] = pd.to_numeric(
                df["normalized_b_factor"], errors="coerce"
            )
        else:
            df["local_flexibility_proxy"] = np.nan
    if "flexibility_method" not in df:
        df["flexibility_method"] = "BFACTOR" if "normalized_b_factor" in df else "UNAVAILABLE"
    df["flexibility_value"] = df["local_flexibility_proxy"]
    df["flexibility_component"] = flexibility_component_from_proxy(df["local_flexibility_proxy"])
    df["rigidity_component"] = rigidity_component_from_flexibility(df["local_flexibility_proxy"])

    if "local_exposed_lys_count" not in df:
        df["local_exposed_lys_count"] = 0
    df["lysine_environment_component"] = lysine_environment_component_from_count(
        df["local_exposed_lys_count"],
        saturation_k=float(lys_config.get("saturation_k", 3.0)),
    )

    if "distance_to_nearest_protected" in df:
        protected_distance = pd.to_numeric(df["distance_to_nearest_protected"], errors="coerce")
        df["protected_site_penalty"] = np.clip(
            (protected_radius - protected_distance.fillna(protected_radius)) / protected_radius,
            0.0,
            1.0,
        )
    else:
        df["protected_site_penalty"] = 0.0

    if "nearest_existing_cys_distance" not in df and "distance_to_existing_cys" in df:
        df["nearest_existing_cys_distance"] = df["distance_to_existing_cys"]
    if "nearest_existing_cys_distance" in df:
        cys_distance = pd.to_numeric(df["nearest_existing_cys_distance"], errors="coerce")
        df["existing_cys_proximity_warning"] = cys_distance <= existing_cys_warning_radius
        distance_penalty = np.clip(
            (existing_cys_penalty_radius - cys_distance.fillna(existing_cys_penalty_radius))
            / existing_cys_penalty_radius,
            0.0,
            1.0,
        )
    else:
        df["existing_cys_proximity_warning"] = False
        distance_penalty = pd.Series(0.0, index=df.index)

    count_col = "existing_cys_count_10A" if "existing_cys_count_10A" in df else None
    count_penalty = (
        pd.Series(0.0, index=df.index)
        if count_col is None
        else np.clip(_series_or_default(df, count_col) / 3.0, 0.0, 1.0)
    )
    df["existing_cys_penalty"] = np.clip(0.7 * distance_penalty + 0.3 * count_penalty, 0.0, 1.0)

    # Single transparent engineering score. Every positive signal is normalized
    # to [0, 1]; nearby cysteines act as a penalty. Protected residues, when supplied,
    # remain available as an informational column but are not part of the MVP score.
    df["stability_score"] = df["stability_component"]
    df["sasa_score"] = df["accessibility_component"]
    df["flexibility_score"] = df["flexibility_component"]
    df["lysine_boost"] = df["lysine_environment_component"]
    df["secondary_structure_penalty"] = _series_or_default(
        df, "secondary_structure_penalty", default=0.0
    ).clip(lower=0.0, upper=1.0)

    ranking_weights = {
        "stability": float(ranking_config.get("score_weights", {}).get("stability", 0.30)),
        "sasa": float(ranking_config.get("score_weights", {}).get("sasa", 0.20)),
        "flexibility": float(ranking_config.get("score_weights", {}).get("flexibility", 0.20)),
        "lysine_boost": float(ranking_config.get("score_weights", {}).get("lysine_boost", 0.10)),
        "existing_cys_penalty": float(
            ranking_config.get("score_weights", {}).get("existing_cys_penalty", 0.10)
        ),
        "secondary_structure_penalty": float(
            ranking_config.get("score_weights", {}).get("secondary_structure_penalty", 0.10)
        ),
    }
    df["final_engineering_score"] = (
        ranking_weights["stability"] * df["stability_score"]
        + ranking_weights["sasa"] * df["sasa_score"]
        + ranking_weights["flexibility"] * df["flexibility_score"]
        + ranking_weights["lysine_boost"] * df["lysine_boost"]
        - ranking_weights["existing_cys_penalty"] * df["existing_cys_penalty"]
        - ranking_weights["secondary_structure_penalty"] * df["secondary_structure_penalty"]
    )

    # Backward-compatible aliases for exported CSV consumers.
    df["cys_site_suitability"] = (
        0.60 * df["stability_score"]
        + 0.35 * df["sasa_score"]
        - 0.10 * df["existing_cys_penalty"]
        - 0.15 * df["protected_site_penalty"]
    )
    df["rigidification_potential"] = (
        0.35 * df["flexibility_score"]
        + 0.40 * df["lysine_boost"]
        + 0.25 * df["sasa_score"]
        - 0.05 * df["existing_cys_penalty"]
        - 0.10 * df["protected_site_penalty"]
    )
    df["cys_suitability_score"] = df["final_engineering_score"]
    df["ranking_formula"] = (
        "final_priority = 0.30*ml_stability + 0.20*relative_exposure + "
        "0.20*flexibility + 0.10*nearby_lys_boost - 0.10*nearby_cys_penalty "
        "- 0.10*secondary_structure_penalty"
    )

    df = df.sort_values(
        ["final_engineering_score", "cys_site_suitability", "predicted_destabilization_ddg"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    df.insert(0, "rank_engineering", range(1, len(df) + 1))
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    return df

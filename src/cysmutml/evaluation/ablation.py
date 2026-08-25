"""Feature-set ablation runner."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from cysmutml.models.train import evaluate_models


def run_ablation(feature_csv: str | Path, results_dir: str | Path = "results") -> pd.DataFrame:
    rows = []
    for feature_set, include_structural in (("physicochemical", False), ("structural", True)):
        metrics, _, _ = evaluate_models(
            feature_csv, results_dir, include_structural=include_structural
        )
        grouped = metrics.groupby("model").agg(
            {"mae": ["mean", "std"], "rmse": ["mean", "std"], "r2": ["mean", "std"]}
        )
        grouped.columns = ["_".join(col).strip("_") for col in grouped.columns]
        grouped = grouped.reset_index()
        grouped.insert(0, "feature_set", feature_set)
        rows.append(grouped)
    out = pd.concat(rows, ignore_index=True)
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    out.to_csv(Path(results_dir) / "ablation_results.csv", index=False)
    return out

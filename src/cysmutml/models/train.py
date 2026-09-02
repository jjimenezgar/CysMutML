"""Model training and group-aware evaluation."""

from __future__ import annotations

import json
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from cysmutml import __version__
from cysmutml.config import load_config
from cysmutml.evaluation.metrics import regression_metrics
from cysmutml.features.build import STRUCTURAL_COLUMNS, feature_columns
from cysmutml.models.pipeline import make_regressors


def git_commit() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


def evaluate_models(
    feature_csv: str | Path,
    results_dir: str | Path = "results",
    config_path: str | Path = "configs/default.yaml",
    include_structural: bool = True,
    model_names: list[str] | None = None,
    group_column: str = "protein_id",
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    config = load_config(config_path)
    df = pd.read_csv(feature_csv, low_memory=False)
    target = "destabilization_ddg_kcal_mol"
    if group_column not in df.columns:
        raise ValueError(f"Grouping column {group_column!r} is missing from the feature table")
    if df[group_column].isna().any():
        raise ValueError(f"Grouping column {group_column!r} contains missing values")
    groups = df[group_column].astype(str)
    n_groups = groups.nunique()
    n_splits = min(int(config["cv_folds"]), n_groups)
    if n_splits < 2:
        raise ValueError("At least two protein groups are required for GroupKFold evaluation")

    numeric, categorical = feature_columns(df, include_structural=include_structural)
    X = df[numeric + categorical]
    y = df[target].astype(float)
    models = make_regressors(numeric, categorical, int(config["random_seed"]))
    if model_names is not None:
        models = {name: models[name] for name in model_names}
    splitter = GroupKFold(n_splits=n_splits)
    metrics_rows = []
    prediction_rows = []

    for model_name, estimator in models.items():
        for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups), start=1):
            fit_start = time.perf_counter()
            estimator.fit(X.iloc[train_idx], y.iloc[train_idx])
            fit_seconds = time.perf_counter() - fit_start
            predict_start = time.perf_counter()
            pred = estimator.predict(X.iloc[test_idx])
            predict_seconds = time.perf_counter() - predict_start
            fold_metrics = regression_metrics(y.iloc[test_idx], pred)
            metrics_rows.append(
                {
                    "model": model_name,
                    "fold": fold,
                    "group_column": group_column,
                    "n_train_groups": int(groups.iloc[train_idx].nunique()),
                    "n_test_groups": int(groups.iloc[test_idx].nunique()),
                    "fit_seconds": fit_seconds,
                    "predict_seconds": predict_seconds,
                    **fold_metrics,
                }
            )
            for row_i, y_hat in zip(test_idx, pred, strict=True):
                prediction_rows.append(
                    {
                        "model": model_name,
                        "fold": fold,
                        "protein_id": df.iloc[row_i]["protein_id"],
                        "group_column": group_column,
                        "group_id": groups.iloc[row_i],
                        "mutation": df.iloc[row_i]["mutation"],
                        "observed": y.iloc[row_i],
                        "predicted": float(y_hat),
                        "residual": float(y.iloc[row_i] - y_hat),
                        "mut_aa": df.iloc[row_i]["mut_aa"],
                    }
                )

    metrics_df = pd.DataFrame(metrics_rows)
    preds_df = pd.DataFrame(prediction_rows)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(results_dir / "regression_cv_metrics.csv", index=False)
    preds_df.to_csv(results_dir / "regression_out_of_fold_predictions.csv", index=False)
    summary = metrics_df.groupby("model")["mae"].mean().sort_values()
    best_model = str(summary.index[0])
    return metrics_df, preds_df, best_model


def evaluate_fast_baselines(
    feature_csv: str | Path,
    results_dir: str | Path = "results/fireprotdb_fast_baselines",
    config_path: str | Path = "configs/default.yaml",
    save_oof: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fast full-data GroupKFold for Dummy and Ridge without OOF prediction files.

    This path is intentionally narrow: it gives an executable scientific baseline
    on large FireProtDB exports while the full model-suite evaluation remains
    available for longer runs.
    """
    config = load_config(config_path)
    df = pd.read_csv(feature_csv, low_memory=False)
    groups = df["protein_id"].fillna("unknown_protein").astype(str)
    n_splits = min(int(config["cv_folds"]), groups.nunique())
    if n_splits < 2:
        raise ValueError("At least two protein groups are required for GroupKFold evaluation")

    numeric, categorical = feature_columns(df, include_structural=False)
    X_num = df[numeric].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    if categorical:
        X_cat = pd.get_dummies(df[categorical].fillna("missing"), dtype=np.float32).to_numpy()
        X = np.hstack([X_num, X_cat])
    else:
        X = X_num
    y = df["destabilization_ddg_kcal_mol"].astype(float).to_numpy()
    is_cys = df["mut_aa"].astype(str).eq("C").to_numpy()

    from sklearn.dummy import DummyRegressor
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline

    models = {
        "dummy_mean": DummyRegressor(strategy="mean"),
        "ridge_fast": make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            Ridge(alpha=1.0, solver="lsqr"),
        ),
    }

    metrics_rows = []
    cys_rows = []
    oof_frames = []
    splitter = GroupKFold(n_splits=n_splits)
    for model_name, estimator in models.items():
        for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups), start=1):
            estimator.fit(X[train_idx], y[train_idx])
            pred = estimator.predict(X[test_idx])
            fold_metrics = regression_metrics(y[test_idx], pred)
            metrics_rows.append(
                {
                    "model": model_name,
                    "fold": fold,
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                    **fold_metrics,
                }
            )
            cys_mask = is_cys[test_idx]
            if int(cys_mask.sum()) >= 2:
                cys_metrics = regression_metrics(y[test_idx][cys_mask], pred[cys_mask])
                cys_rows.append(
                    {
                        "model": model_name,
                        "fold": fold,
                        "n_cys_test": int(cys_mask.sum()),
                        **cys_metrics,
                    }
                )
            else:
                cys_rows.append(
                    {
                        "model": model_name,
                        "fold": fold,
                        "n_cys_test": int(cys_mask.sum()),
                        "mae": np.nan,
                        "rmse": np.nan,
                        "r2": np.nan,
                        "pearson": np.nan,
                        "spearman": np.nan,
                    }
                )
            if save_oof:
                oof_frames.append(
                    pd.DataFrame(
                        {
                            "model": model_name,
                            "fold": fold,
                            "protein_id": df.iloc[test_idx]["protein_id"].to_numpy(),
                            "mutation": df.iloc[test_idx]["mutation"].to_numpy(),
                            "wt_aa": df.iloc[test_idx]["wt_aa"].to_numpy(),
                            "mut_aa": df.iloc[test_idx]["mut_aa"].to_numpy(),
                            "observed": y[test_idx],
                            "predicted": pred,
                            "residual": y[test_idx] - pred,
                            "absolute_residual": np.abs(y[test_idx] - pred),
                        }
                    )
                )

    metrics_df = pd.DataFrame(metrics_rows)
    cys_df = pd.DataFrame(cys_rows)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(results_dir / "fast_baseline_cv_metrics.csv", index=False)
    cys_df.to_csv(results_dir / "fast_baseline_cys_metrics.csv", index=False)
    if save_oof:
        pd.concat(oof_frames, ignore_index=True).to_csv(
            results_dir / "fast_baseline_out_of_fold_predictions.csv", index=False
        )
    return metrics_df, cys_df


def evaluate_physchem_model_comparison(
    feature_csv: str | Path,
    results_dir: str | Path = "results/physchem_model_comparison",
    config_path: str | Path = "configs/default.yaml",
    model_names: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare compact physicochemical-only regressors with grouped CV."""
    selected_models = model_names or ["dummy_mean", "ridge", "hist_gradient_boosting"]
    metrics, preds, _ = evaluate_models(
        feature_csv,
        results_dir,
        config_path=config_path,
        include_structural=False,
        model_names=selected_models,
    )
    cys = preds[preds["mut_aa"].astype(str).eq("C")].copy()
    cys_rows = []
    for (model, fold), fold_df in cys.groupby(["model", "fold"]):
        if len(fold_df) >= 2:
            values = regression_metrics(fold_df["observed"], fold_df["predicted"])
        else:
            values = {
                "mae": np.nan,
                "rmse": np.nan,
                "r2": np.nan,
                "pearson": np.nan,
                "spearman": np.nan,
            }
        cys_rows.append({"model": model, "fold": int(fold), "n_cys_test": len(fold_df), **values})
    cys_metrics = pd.DataFrame(cys_rows)
    results_dir = Path(results_dir)
    cys_metrics.to_csv(results_dir / "cys_specific_metrics.csv", index=False)

    summary = (
        metrics.groupby("model")[["mae", "rmse", "r2", "pearson", "spearman"]]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary.to_csv(results_dir / "physchem_model_comparison_summary.csv", index=False)
    return metrics, cys_metrics, preds


def train_final_model(
    feature_csv: str | Path,
    model_path: str | Path = "models/cysmutml_model.joblib",
    metadata_path: str | Path = "models/model_metadata.json",
    config_path: str | Path = "configs/default.yaml",
    model_name: str | None = None,
    include_structural: bool = False,
) -> dict:
    config = load_config(config_path)
    df = pd.read_csv(feature_csv, low_memory=False)
    numeric, categorical = feature_columns(df, include_structural=include_structural)
    models = make_regressors(numeric, categorical, int(config["random_seed"]))
    selected = model_name or "ridge"
    estimator = models[selected]
    X = df[numeric + categorical]
    y = df["destabilization_ddg_kcal_mol"].astype(float)
    estimator.fit(X, y)
    selected_features = numeric + categorical
    actual_feature_configuration = (
        "structural"
        if any(column in selected_features for column in STRUCTURAL_COLUMNS)
        else "physicochemical"
    )
    artifact = {
        "model": estimator,
        "numeric_features": numeric,
        "categorical_features": categorical,
        "training_feature_ranges": {
            col: {"min": float(np.nanmin(df[col])), "max": float(np.nanmax(df[col]))}
            for col in numeric
            if np.issubdtype(df[col].dtype, np.number)
        },
        "target_training_distribution": {
            "min": float(np.nanmin(y)),
            "max": float(np.nanmax(y)),
            "p05": float(np.nanpercentile(y, 5)),
            "p50": float(np.nanpercentile(y, 50)),
            "p95": float(np.nanpercentile(y, 95)),
        },
    }
    model_path = Path(model_path)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(artifact, model_path)
    metadata = {
        "training_date": datetime.now(timezone.utc).isoformat(),
        "dataset_source": "synthetic test data"
        if "synthetic" in str(feature_csv)
        else "FireProtDB v2.0 API CSV export",
        "feature_configuration": actual_feature_configuration,
        "target_definition": "destabilization_ddg_kcal_mol; positive means destabilizing",
        "model_class": selected,
        "hyperparameters": estimator.named_steps["model"].get_params(),
        "software_versions": {"python": platform.python_version(), "cysmutml": __version__},
        "git_commit": git_commit(),
    }
    Path(metadata_path).write_text(json.dumps(metadata, indent=2, default=str))
    return metadata

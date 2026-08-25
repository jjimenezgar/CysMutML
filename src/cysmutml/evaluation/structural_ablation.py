"""Controlled structural ablation on identical mapped rows and folds."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupKFold

from cysmutml.config import load_config
from cysmutml.evaluation.metrics import regression_metrics
from cysmutml.features.build import feature_columns
from cysmutml.models.pipeline import make_regressors


def make_structural_cv_folds(
    structural_csv: str | Path = "data/processed/fireprotdb_structural_features.csv",
    output_csv: str | Path = "data/processed/structural_cv_folds.csv",
    config_path: str | Path = "configs/default.yaml",
) -> pd.DataFrame:
    config = load_config(config_path)
    df = pd.read_csv(structural_csv, low_memory=False).reset_index(drop=True)
    groups = df["protein_id"].fillna("unknown_protein").astype(str)
    n_splits = min(int(config["cv_folds"]), groups.nunique())
    if n_splits < 2:
        raise ValueError("At least two protein groups are required")
    fold_df = pd.DataFrame({"row_id": df.index, "protein_id": groups, "fold": -1})
    splitter = GroupKFold(n_splits=n_splits)
    for fold, (_, test_idx) in enumerate(
        splitter.split(df, df["destabilization_ddg_kcal_mol"], groups), start=1
    ):
        fold_df.loc[test_idx, "fold"] = fold
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    fold_df.to_csv(output_csv, index=False)
    return fold_df


def _assert_no_leakage(folds: pd.DataFrame) -> None:
    for fold in sorted(folds["fold"].unique()):
        train = set(folds.loc[folds["fold"] != fold, "protein_id"])
        test = set(folds.loc[folds["fold"] == fold, "protein_id"])
        overlap = train.intersection(test)
        if overlap:
            raise ValueError(f"Protein leakage in fold {fold}: {sorted(overlap)[:5]}")


def run_structural_ablation(
    structural_csv: str | Path = "data/processed/fireprotdb_structural_features.csv",
    folds_csv: str | Path = "data/processed/structural_cv_folds.csv",
    results_dir: str | Path = "results/structural_ablation",
    config_path: str | Path = "configs/default.yaml",
) -> dict[str, pd.DataFrame]:
    config = load_config(config_path)
    df = pd.read_csv(structural_csv, low_memory=False).reset_index(drop=True)
    if len(df) < 10:
        raise ValueError("Structural dataset is too small for ablation")
    folds = make_structural_cv_folds(structural_csv, folds_csv, config_path)
    _assert_no_leakage(folds)
    y = df["destabilization_ddg_kcal_mol"].astype(float)
    feature_sets = {
        "physchem_only_structural_subset": False,
        "physchem_plus_structure": True,
    }
    metrics_rows = []
    prediction_frames = []
    models_to_run = ["dummy_mean", "ridge", "hist_gradient_boosting"]
    for feature_set, include_structure in feature_sets.items():
        numeric, categorical = feature_columns(df, include_structural=include_structure)
        X = df[numeric + categorical]
        models = make_regressors(numeric, categorical, int(config["random_seed"]))
        for model_name in models_to_run:
            estimator = models[model_name]
            for fold in sorted(folds["fold"].unique()):
                train_idx = folds.index[folds["fold"] != fold].to_numpy()
                test_idx = folds.index[folds["fold"] == fold].to_numpy()
                estimator.fit(X.iloc[train_idx], y.iloc[train_idx])
                pred = estimator.predict(X.iloc[test_idx])
                metrics_rows.append(
                    {
                        "feature_set": feature_set,
                        "model": model_name,
                        "fold": int(fold),
                        "n_train": int(len(train_idx)),
                        "n_test": int(len(test_idx)),
                        **regression_metrics(y.iloc[test_idx], pred),
                    }
                )
                prediction_frames.append(
                    pd.DataFrame(
                        {
                            "feature_set": feature_set,
                            "model": model_name,
                            "fold": int(fold),
                            "row_id": test_idx,
                            "protein_id": df.iloc[test_idx]["protein_id"].to_numpy(),
                            "mutation": df.iloc[test_idx]["original_mutation"].to_numpy(),
                            "mut_aa": df.iloc[test_idx]["mut_aa"].to_numpy(),
                            "observed": y.iloc[test_idx].to_numpy(),
                            "predicted": pred,
                        }
                    )
                )
    results_dir = Path(results_dir)
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    metrics = pd.DataFrame(metrics_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    predictions["residual"] = predictions["observed"] - predictions["predicted"]
    predictions["absolute_residual"] = predictions["residual"].abs()
    summary = (
        metrics.groupby(["feature_set", "model"])
        .agg({metric: ["mean", "std"] for metric in ["mae", "rmse", "r2", "pearson", "spearman"]})
        .reset_index()
    )
    summary.columns = ["_".join(col).strip("_") for col in summary.columns]
    metrics.to_csv(results_dir / "fold_metrics.csv", index=False)
    summary.to_csv(results_dir / "ablation_metrics.csv", index=False)
    predictions.to_csv(results_dir / "out_of_fold_predictions.csv", index=False)

    paired_rows = []
    for model_name in models_to_run:
        phys = metrics[
            (metrics["feature_set"] == "physchem_only_structural_subset")
            & (metrics["model"] == model_name)
        ]
        struct = metrics[
            (metrics["feature_set"] == "physchem_plus_structure")
            & (metrics["model"] == model_name)
        ]
        merged = phys.merge(struct, on=["model", "fold"], suffixes=("_physchem", "_structure"))
        for _, row in merged.iterrows():
            paired_rows.append(
                {
                    "model": model_name,
                    "fold": int(row["fold"]),
                    "delta_mae_structure_minus_physchem": row["mae_structure"]
                    - row["mae_physchem"],
                    "delta_rmse_structure_minus_physchem": row["rmse_structure"]
                    - row["rmse_physchem"],
                    "delta_r2_structure_minus_physchem": row["r2_structure"] - row["r2_physchem"],
                }
            )
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(results_dir / "paired_ablation_comparison.csv", index=False)

    cys_rows = []
    for (feature_set, model_name), sub in predictions.groupby(["feature_set", "model"]):
        cys = sub[sub["mut_aa"] == "C"]
        metrics_dict = (
            regression_metrics(cys["observed"], cys["predicted"]) if len(cys) >= 2 else {}
        )
        cys_rows.append(
            {
                "feature_set": feature_set,
                "model": model_name,
                "n_cys_observations": int(len(cys)),
                "n_cys_proteins": int(cys["protein_id"].nunique()) if len(cys) else 0,
                **metrics_dict,
            }
        )
    cys_metrics = pd.DataFrame(cys_rows)
    cys_metrics.to_csv(results_dir / "cys_specific_metrics.csv", index=False)

    _plot_ablation(summary, figures_dir)
    _write_permutation_importance(df, folds, results_dir, config_path)
    return {"fold_metrics": metrics, "ablation_metrics": summary, "cys_metrics": cys_metrics}


def _plot_ablation(summary: pd.DataFrame, figures_dir: Path) -> None:
    for metric in ["mae", "r2"]:
        column = f"{metric}_mean"
        plt.figure(figsize=(8, 4))
        labels = summary["feature_set"] + "\n" + summary["model"]
        plt.bar(range(len(summary)), summary[column], color="#4c78a8")
        plt.xticks(range(len(summary)), labels, rotation=45, ha="right")
        plt.ylabel(metric.upper())
        plt.title(f"Structural Ablation {metric.upper()}")
        plt.tight_layout()
        plt.savefig(figures_dir / f"{metric}_comparison_by_feature_set.png", dpi=160)
        plt.close()


def _write_permutation_importance(
    df: pd.DataFrame, folds: pd.DataFrame, results_dir: Path, config_path: str | Path
) -> None:
    config = load_config(config_path)
    numeric, categorical = feature_columns(df, include_structural=True)
    X = df[numeric + categorical]
    y = df["destabilization_ddg_kcal_mol"].astype(float)
    models = make_regressors(numeric, categorical, int(config["random_seed"]))
    estimator = models["hist_gradient_boosting"]
    fold = sorted(folds["fold"].unique())[0]
    train_idx = folds.index[folds["fold"] != fold].to_numpy()
    test_idx = folds.index[folds["fold"] == fold].to_numpy()
    estimator.fit(X.iloc[train_idx], y.iloc[train_idx])
    importance = permutation_importance(
        estimator,
        X.iloc[test_idx],
        y.iloc[test_idx],
        n_repeats=5,
        random_state=int(config["random_seed"]),
        scoring="neg_mean_absolute_error",
    )
    out = pd.DataFrame(
        {
            "feature": numeric + categorical,
            "importance_mean": importance.importances_mean,
            "importance_std": importance.importances_std,
        }
    ).sort_values("importance_mean", ascending=False)
    out.to_csv(results_dir / "permutation_importance.csv", index=False)
    top = out.head(15)
    plt.figure(figsize=(7, 4))
    plt.barh(top["feature"][::-1], top["importance_mean"][::-1], color="#59a14f")
    plt.xlabel("Permutation importance, neg-MAE decrease")
    plt.tight_layout()
    plt.savefig(results_dir / "figures" / "permutation_importance.png", dpi=160)
    plt.close()

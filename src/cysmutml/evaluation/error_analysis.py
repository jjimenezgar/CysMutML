"""Error-analysis reports for out-of-fold predictions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from cysmutml.evaluation.metrics import regression_metrics


def write_error_analysis(
    predictions_csv: str | Path,
    output_dir: str | Path = "results/fireprotdb_fast_baselines",
    model_name: str = "ridge_fast",
) -> dict[str, Path]:
    output_dir = Path(output_dir)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(predictions_csv, low_memory=False)
    df = df[df["model"] == model_name].copy()
    if df.empty:
        raise ValueError(f"No rows found for model {model_name!r}")

    df["absolute_residual"] = df["absolute_residual"].astype(float)
    worst = df.sort_values("absolute_residual", ascending=False).head(100)
    worst_path = output_dir / f"{model_name}_largest_residuals.csv"
    worst.to_csv(worst_path, index=False)

    by_mut = []
    for mut_aa, sub in df.groupby("mut_aa"):
        if len(sub) < 2:
            continue
        metrics = regression_metrics(sub["observed"], sub["predicted"])
        by_mut.append({"mut_aa": mut_aa, "n": len(sub), **metrics})
    by_mut_df = pd.DataFrame(by_mut).sort_values("mae")
    by_mut_path = output_dir / f"{model_name}_metrics_by_mutant_aa.csv"
    by_mut_df.to_csv(by_mut_path, index=False)

    cys = df[df["mut_aa"] == "C"].copy()
    cys_path = output_dir / f"{model_name}_cys_out_of_fold_predictions.csv"
    cys.to_csv(cys_path, index=False)

    pred_plot = figures_dir / f"{model_name}_predicted_vs_observed.png"
    sample = df.sample(min(len(df), 100_000), random_state=42)
    plt.figure(figsize=(5, 5))
    plt.scatter(sample["observed"], sample["predicted"], s=3, alpha=0.15)
    low = min(sample["observed"].min(), sample["predicted"].min())
    high = max(sample["observed"].max(), sample["predicted"].max())
    plt.plot([low, high], [low, high], color="black", linewidth=1)
    plt.xlabel("Observed destabilization DDG (kcal/mol)")
    plt.ylabel("Predicted destabilization DDG (kcal/mol)")
    plt.title("Ridge OOF Predictions")
    plt.tight_layout()
    plt.savefig(pred_plot, dpi=160)
    plt.close()

    residual_plot = figures_dir / f"{model_name}_residual_distribution.png"
    plt.figure(figsize=(6, 4))
    plt.hist(df["residual"], bins=100, color="#4c78a8")
    plt.xlabel("Residual: observed - predicted (kcal/mol)")
    plt.ylabel("Count")
    plt.title("Ridge OOF Residuals")
    plt.tight_layout()
    plt.savefig(residual_plot, dpi=160)
    plt.close()

    mut_plot = figures_dir / f"{model_name}_mae_by_mutant_aa.png"
    plt.figure(figsize=(7, 4))
    plt.bar(by_mut_df["mut_aa"], by_mut_df["mae"], color="#59a14f")
    plt.xlabel("Mutant amino acid")
    plt.ylabel("MAE (kcal/mol)")
    plt.title("Ridge OOF MAE by Mutant Amino Acid")
    plt.tight_layout()
    plt.savefig(mut_plot, dpi=160)
    plt.close()

    return {
        "largest_residuals": worst_path,
        "metrics_by_mutant_aa": by_mut_path,
        "cys_predictions": cys_path,
        "predicted_vs_observed": pred_plot,
        "residual_distribution": residual_plot,
        "mae_by_mutant_aa": mut_plot,
    }


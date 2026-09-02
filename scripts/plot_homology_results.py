"""Generate compact portfolio figures from the homology MVP benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLORS = {
    "protein_grouped": "#2563eb",
    "homology_clustered": "#f59e0b",
}


def _grouped_mae_plot(table: pd.DataFrame, output: Path, title: str) -> None:
    summary = table.groupby(["model", "split_strategy"])["mae"].agg(["mean", "std"])
    models = list(summary.index.get_level_values("model").unique())
    strategies = ["protein_grouped", "homology_clustered"]
    positions = np.arange(len(models))
    width = 0.36

    figure, axis = plt.subplots(figsize=(9, 4.8))
    for index, strategy in enumerate(strategies):
        values = summary.xs(strategy, level="split_strategy").reindex(models)
        axis.bar(
            positions + (index - 0.5) * width,
            values["mean"],
            width,
            yerr=values["std"].fillna(0),
            capsize=4,
            label=strategy.replace("_", " ").title(),
            color=COLORS[strategy],
        )
    axis.set_xticks(positions, [name.replace("_", " ").title() for name in models], rotation=15)
    axis.set_ylabel("MAE (kcal/mol)")
    axis.set_title(title)
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def _importance_plot(table: pd.DataFrame, output: Path) -> None:
    models = list(table["model"].unique())
    figure, axes = plt.subplots(1, len(models), figsize=(6 * len(models), 5), squeeze=False)
    for axis, model in zip(axes[0], models, strict=True):
        subset = table[table["model"] == model].nlargest(12, "importance_mean_mae")
        subset = subset.sort_values("importance_mean_mae")
        axis.barh(subset["feature"], subset["importance_mean_mae"], color="#14b8a6")
        axis.set_title(model.replace("_", " ").title())
        axis.set_xlabel("Increase in MAE after permutation")
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle("Held-out permutation importance")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def generate_figures(results_dir: str | Path) -> list[Path]:
    results_dir = Path(results_dir)
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    overall = pd.read_csv(results_dir / "split_comparison_fold_metrics.csv")
    cys = pd.read_csv(results_dir / "split_comparison_cys_metrics.csv")
    importance = pd.read_csv(results_dir / "tree_permutation_importance.csv")

    outputs = [
        figures_dir / "mae_by_split_and_model.png",
        figures_dir / "cys_mae_by_split_and_model.png",
        figures_dir / "tree_permutation_importance.png",
    ]
    _grouped_mae_plot(overall, outputs[0], "Protein vs homology-grouped validation")
    _grouped_mae_plot(cys, outputs[1], "X→Cys performance by validation strategy")
    _importance_plot(importance, outputs[2])
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results/homology_validation")
    args = parser.parse_args()
    for path in generate_figures(args.results_dir):
        print(path)


if __name__ == "__main__":
    main()

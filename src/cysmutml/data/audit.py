"""FireProtDB data-quality audits for the structural milestone."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from cysmutml.mutations import MUTATION_RE

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def split_pdb_ids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    parts = re.split(r"[,;|\s]+", str(value).strip())
    return [part.lower() for part in parts if PDB_ID_RE.match(part)]


def _joined_unique(values: pd.Series, limit: int = 12) -> str:
    unique = sorted({str(v) for v in values.dropna().unique() if str(v).strip()})
    if len(unique) > limit:
        return "; ".join(unique[:limit]) + f"; ... ({len(unique)} unique)"
    return "; ".join(unique)


def audit_duplicates_and_aggregate(
    processed_csv: str | Path = "data/processed/fireprotdb_mutations.csv",
    duplicate_audit_csv: str | Path = "reports/duplicate_measurement_audit.csv",
    duplicate_summary_md: str | Path = "reports/duplicate_measurement_summary.md",
    aggregated_csv: str | Path = "data/processed/fireprotdb_mutations_aggregated.csv",
) -> dict[str, int]:
    df = pd.read_csv(processed_csv, low_memory=False)
    key = ["protein_id", "wt_aa", "position", "mut_aa"]
    df["destabilization_ddg_kcal_mol"] = df["destabilization_ddg_kcal_mol"].astype(float)
    grouped = df.groupby(key, dropna=False)
    aggregate_df = grouped.agg(
        n_measurements=("destabilization_ddg_kcal_mol", "size"),
        median_destabilization_ddg=("destabilization_ddg_kcal_mol", "median"),
        mean_destabilization_ddg=("destabilization_ddg_kcal_mol", "mean"),
        std_destabilization_ddg=("destabilization_ddg_kcal_mol", "std"),
        min_destabilization_ddg=("destabilization_ddg_kcal_mol", "min"),
        max_destabilization_ddg=("destabilization_ddg_kcal_mol", "max"),
        pdb_id_values=("pdb_id", "first"),
    ).reset_index()
    aggregate_df["std_destabilization_ddg"] = aggregate_df["std_destabilization_ddg"].fillna(0.0)
    aggregate_df["mutation"] = (
        aggregate_df["wt_aa"].astype(str)
        + aggregate_df["position"].astype(int).astype(str)
        + aggregate_df["mut_aa"].astype(str)
    )

    duplicate_df = aggregate_df[aggregate_df["n_measurements"] > 1].copy()
    repeated_keys = duplicate_df[key]
    repeated = df.merge(repeated_keys, on=key, how="inner")
    raw_values = (
        repeated.groupby(key, dropna=False)["destabilization_ddg_kcal_mol"]
        .apply(lambda values: ";".join(f"{float(v):.6g}" for v in values))
        .rename("raw_destabilization_ddg_values")
        .reset_index()
    )
    duplicate_df = duplicate_df.merge(raw_values, on=key, how="left")
    for column, output in (
        ("exp_temperature", "temperature_values"),
        ("ph", "ph_values"),
        ("method", "method_values"),
        ("measure", "measure_values"),
        ("source_dataset", "source_dataset_values"),
        ("pdb_id", "pdb_id_values"),
    ):
        values = (
            repeated.groupby(key, dropna=False)[column]
            .apply(_joined_unique)
            .rename(output)
            .reset_index()
        )
        duplicate_df = duplicate_df.drop(columns=[output], errors="ignore").merge(
            values, on=key, how="left"
        )
    exact = duplicate_df["min_destabilization_ddg"].eq(duplicate_df["max_destabilization_ddg"])
    condition_variation = (
        duplicate_df["temperature_values"].str.contains(";", regex=False, na=False)
        | duplicate_df["ph_values"].str.contains(";", regex=False, na=False)
        | duplicate_df["method_values"].str.contains(";", regex=False, na=False)
        | duplicate_df["measure_values"].str.contains(";", regex=False, na=False)
    )
    duplicate_df["duplicate_category"] = "replicate_or_unclear_repeat"
    duplicate_df.loc[condition_variation, "duplicate_category"] = "condition_or_assay_variation"
    duplicate_df.loc[exact, "duplicate_category"] = "exact_or_near_exact_duplicate"
    Path(duplicate_audit_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(aggregated_csv).parent.mkdir(parents=True, exist_ok=True)
    duplicate_df.to_csv(duplicate_audit_csv, index=False)
    ordered = [
        "protein_id",
        "mutation",
        "wt_aa",
        "position",
        "mut_aa",
        "median_destabilization_ddg",
        "mean_destabilization_ddg",
        "std_destabilization_ddg",
        "min_destabilization_ddg",
        "max_destabilization_ddg",
        "n_measurements",
        "pdb_id_values",
    ]
    aggregate_df[ordered].to_csv(aggregated_csv, index=False)

    category_counts = duplicate_df["duplicate_category"].value_counts().to_dict()
    summary = [
        "# Duplicate Measurement Summary",
        "",
        "Repeated protein/mutation measurements are preserved in the measurement-level dataset.",
        "The aggregated dataset uses the median `destabilization_ddg_kcal_mol`",
        "as the primary target.",
        "",
        f"- Measurement-level rows: {len(df):,}",
        f"- Aggregated mutation rows: {len(aggregate_df):,}",
        f"- Repeated protein/mutation groups: {len(duplicate_df):,}",
        f"- Repeated measurement rows: {int(df.duplicated(key, keep=False).sum()):,}",
        "",
        "## Duplicate Categories",
        "",
    ]
    for category, count in category_counts.items():
        summary.append(f"- {category}: {count:,}")
    summary.extend(
        [
            "",
            "Categories are heuristic audit labels, not automatic exclusion rules.",
            "No measurements are removed solely because repeated observations differ.",
        ]
    )
    Path(duplicate_summary_md).write_text("\n".join(summary) + "\n")
    return {
        "measurement_rows": len(df),
        "aggregated_rows": len(aggregate_df),
        "duplicate_groups": len(duplicate_df),
    }


def audit_ddg_distribution(
    processed_csv: str | Path = "data/processed/fireprotdb_mutations.csv",
    report_md: str | Path = "reports/ddg_distribution_audit.md",
    extreme_csv: str | Path = "reports/extreme_ddg_records.csv",
    figures_dir: str | Path = "results/data_audit",
) -> dict[str, float]:
    df = pd.read_csv(processed_csv, low_memory=False)
    target = df["destabilization_ddg_kcal_mol"].astype(float)
    percentiles = [0.1, 1, 5, 25, 50, 75, 95, 99, 99.9]
    stats = {
        "minimum": float(target.min()),
        "maximum": float(target.max()),
        "mean": float(target.mean()),
        "median": float(target.median()),
        "standard_deviation": float(target.std()),
    }
    stats.update({f"p{p}": float(target.quantile(p / 100.0)) for p in percentiles})

    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 4))
    plt.hist(target, bins=200, color="#4c78a8")
    plt.xlabel("Destabilization DDG (kcal/mol)")
    plt.ylabel("Count")
    plt.title("FireProtDB DDG Distribution")
    plt.tight_layout()
    plt.savefig(figures_dir / "ddg_histogram.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 2.5))
    plt.boxplot(target, vert=False, showfliers=True)
    plt.xlabel("Destabilization DDG (kcal/mol)")
    plt.title("FireProtDB DDG Boxplot")
    plt.tight_layout()
    plt.savefig(figures_dir / "ddg_boxplot.png", dpi=160)
    plt.close()

    lo, hi = stats["p1"], stats["p99"]
    zoom = target[(target >= lo) & (target <= hi)]
    plt.figure(figsize=(7, 4))
    plt.hist(zoom, bins=120, color="#59a14f")
    plt.xlabel("Destabilization DDG (kcal/mol)")
    plt.ylabel("Count")
    plt.title("FireProtDB DDG Distribution, 1st-99th Percentile")
    plt.tight_layout()
    plt.savefig(figures_dir / "ddg_log_frequency_or_zoomed_distribution.png", dpi=160)
    plt.close()

    extreme = pd.concat(
        [
            df.nsmallest(100, "destabilization_ddg_kcal_mol"),
            df.nlargest(100, "destabilization_ddg_kcal_mol"),
        ],
        ignore_index=True,
    )
    extreme.to_csv(extreme_csv, index=False)

    lines = [
        "# DDG Distribution Audit",
        "",
        "Outliers are audited but not removed from the main analysis.",
        "",
        "| Statistic | Value |",
        "|---|---:|",
    ]
    for name, value in stats.items():
        lines.append(f"| {name} | {value:.6g} |")
    lines.extend(
        [
            "",
            "Figures:",
            "",
            "- `results/data_audit/ddg_histogram.png`",
            "- `results/data_audit/ddg_boxplot.png`",
            "- `results/data_audit/ddg_log_frequency_or_zoomed_distribution.png`",
            "",
            f"Extreme records are saved in `{extreme_csv}`.",
        ]
    )
    Path(report_md).write_text("\n".join(lines) + "\n")
    return stats


def summarize_structure_candidates(
    raw_csv: str | Path = "data/raw/fireprotdb.csv",
    output_md: str | Path = "reports/structure_candidate_summary.md",
    output_csv: str | Path = "reports/structure_candidate_records.csv",
) -> dict[str, int]:
    raw = pd.read_csv(raw_csv, low_memory=False)
    simple_mut = raw["SUBSTITUTION"].astype(str).str.match(MUTATION_RE)
    valid_ddg = pd.to_numeric(raw["DDG"], errors="coerce").notna()
    valid = raw[simple_mut & valid_ddg].copy()
    records = []
    for idx, row in valid.iterrows():
        for pdb_id in split_pdb_ids(row["WWPDB"]):
            records.append(
                {
                    "source_row_index": int(idx),
                    "protein_id": row["PROTEIN"],
                    "pdb_id": pdb_id,
                    "uniprot_id": row["UNIPROTKB"],
                    "mutation": row["SUBSTITUTION"],
                    "destabilization_ddg_kcal_mol": pd.to_numeric(row["DDG"], errors="coerce"),
                }
            )
    candidate_records = pd.DataFrame(records)
    candidate_records.to_csv(output_csv, index=False)
    chain_col = "CHAIN" if "CHAIN" in raw.columns else None
    chain_count = (
        int(raw.loc[candidate_records["source_row_index"], chain_col].notna().sum())
        if chain_col and not candidate_records.empty
        else 0
    )

    summary = {
        "records_with_pdb": int(candidate_records["pdb_id"].notna().sum()),
        "records_without_pdb": int(valid["WWPDB"].isna().sum()),
        "unique_pdbs": int(candidate_records["pdb_id"].nunique(dropna=True)),
        "records_with_chain": chain_count,
        "records_without_chain": int(len(candidate_records) - chain_count),
        "unique_proteins_with_pdb": int(candidate_records["protein_id"].nunique(dropna=True)),
    }
    lines = [
        "# Structure Candidate Summary",
        "",
        "Initial candidates are valid single-substitution DDG rows with a PDB identifier.",
        "",
        f"- Records with PDB identifier: {summary['records_with_pdb']:,}",
        f"- Records without PDB identifier: {summary['records_without_pdb']:,}",
        f"- Unique PDB IDs: {summary['unique_pdbs']:,}",
        f"- Records with chain metadata: {summary['records_with_chain']:,}",
        f"- Records without chain metadata: {summary['records_without_chain']:,}",
        "- Unique proteins represented among PDB candidates: "
        f"{summary['unique_proteins_with_pdb']:,}",
        "",
        "The downloaded FireProtDB CSV does not provide reliable chain metadata, so chain",
        "selection must be performed by sequence alignment and ambiguous cases must fail.",
    ]
    Path(output_md).write_text("\n".join(lines) + "\n")
    return summary

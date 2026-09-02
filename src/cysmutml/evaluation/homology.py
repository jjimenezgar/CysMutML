"""Homology-aware grouping and split-comparison utilities."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from sklearn.model_selection import GroupKFold

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def _normalized_sequence(value: object) -> str:
    sequence = "".join(str(value).split()).upper()
    if not sequence or sequence in {"NAN", "NONE"}:
        raise ValueError("Sequence is missing")
    invalid = sorted(set(sequence) - VALID_AA)
    if invalid:
        raise ValueError(f"Sequence contains unsupported symbols: {''.join(invalid)}")
    return sequence


def unique_protein_sequences(table: pd.DataFrame) -> pd.DataFrame:
    """Return one verified canonical sequence per protein."""
    required = {"protein_id", "canonical_sequence"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns required for homology clustering: {sorted(missing)}")

    rows = []
    for protein_id, group in table.groupby("protein_id", dropna=False):
        sequences = {
            _normalized_sequence(value)
            for value in group["canonical_sequence"].dropna()
            if str(value).strip()
        }
        if not sequences:
            continue
        if len(sequences) != 1:
            raise ValueError(f"Protein {protein_id!r} has {len(sequences)} canonical sequences")
        rows.append({"protein_id": str(protein_id), "canonical_sequence": sequences.pop()})

    if len(rows) < 2:
        raise ValueError("At least two proteins with canonical sequences are required")
    return pd.DataFrame(rows).sort_values("protein_id").reset_index(drop=True)


def validate_cluster_mapping(mapping: pd.DataFrame) -> pd.DataFrame:
    """Validate a many-proteins-to-one-cluster mapping."""
    required = {"protein_id", "sequence_cluster"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"Cluster mapping is missing columns: {sorted(missing)}")

    clean = mapping.copy()
    clean["protein_id"] = clean["protein_id"].astype(str)
    clean["sequence_cluster"] = clean["sequence_cluster"].astype(str)
    if clean[["protein_id", "sequence_cluster"]].isna().any().any():
        raise ValueError("Cluster mapping contains missing identifiers")
    if clean["protein_id"].str.strip().eq("").any():
        raise ValueError("Cluster mapping contains empty protein identifiers")
    if clean["sequence_cluster"].str.strip().eq("").any():
        raise ValueError("Cluster mapping contains empty cluster identifiers")

    conflicts = clean.groupby("protein_id")["sequence_cluster"].nunique()
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        raise ValueError(
            "Proteins assigned to multiple clusters: "
            + ", ".join(conflicts.index.astype(str).tolist()[:10])
        )
    return clean.drop_duplicates("protein_id").reset_index(drop=True)


def attach_sequence_clusters(features: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """Attach clusters and fail rather than silently falling back to protein grouping."""
    clean_mapping = validate_cluster_mapping(mapping)
    attached = features.copy()
    attached["protein_id"] = attached["protein_id"].astype(str)
    attached = attached.merge(
        clean_mapping[["protein_id", "sequence_cluster"]],
        on="protein_id",
        how="left",
        validate="many_to_one",
    )
    missing = sorted(attached.loc[attached["sequence_cluster"].isna(), "protein_id"].unique())
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(
            f"Cluster mapping does not cover {len(missing)} proteins; first missing: {preview}"
        )
    return attached


def grouped_fold_assignments(
    table: pd.DataFrame,
    group_column: str,
    n_splits: int,
) -> pd.DataFrame:
    """Assign rows to folds and verify complete group isolation."""
    if group_column not in table:
        raise ValueError(f"Grouping column {group_column!r} is missing")
    groups = table[group_column].astype(str)
    unique_groups = groups.nunique()
    folds = min(int(n_splits), unique_groups)
    if folds < 2:
        raise ValueError("At least two unique groups are required")

    assignment = pd.Series(index=table.index, dtype="int64")
    splitter = GroupKFold(n_splits=folds)
    for fold, (_, test_idx) in enumerate(splitter.split(table, groups=groups), start=1):
        assignment.iloc[test_idx] = fold

    result = table[["protein_id", group_column]].copy()
    result["fold"] = assignment.astype(int)
    cluster_fold_counts = result.groupby(group_column)["fold"].nunique()
    if not cluster_fold_counts.eq(1).all():
        raise RuntimeError("A group was assigned to more than one test fold")
    return result


def build_mmseqs_cluster_map(
    input_csv: str | Path,
    output_csv: str | Path,
    min_sequence_identity: float = 0.30,
    coverage: float = 0.80,
    mmseqs_binary: str = "mmseqs",
    work_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Cluster canonical sequences with MMseqs2 easy-cluster."""
    if not 0 < min_sequence_identity <= 1:
        raise ValueError("min_sequence_identity must be in (0, 1]")
    if not 0 < coverage <= 1:
        raise ValueError("coverage must be in (0, 1]")

    executable = shutil.which(mmseqs_binary)
    if executable is None:
        raise FileNotFoundError(
            f"Could not find {mmseqs_binary!r}. Install MMseqs2 or pass --mmseqs-binary."
        )

    source = pd.read_csv(input_csv, low_memory=False)
    proteins = unique_protein_sequences(source)
    proteins["mmseqs_id"] = [f"seq_{index:06d}" for index in range(len(proteins))]
    token_to_protein = dict(zip(proteins["mmseqs_id"], proteins["protein_id"], strict=True))

    temporary = None
    if work_dir is None:
        temporary = tempfile.TemporaryDirectory(prefix="cysmutml_mmseqs_")
        workspace = Path(temporary.name)
    else:
        workspace = Path(work_dir)
        workspace.mkdir(parents=True, exist_ok=True)

    try:
        fasta_path = workspace / "proteins.fasta"
        fasta_text = "".join(
            f">{row.mmseqs_id}\n{row.canonical_sequence}\n"
            for row in proteins.itertuples(index=False)
        )
        fasta_path.write_text(fasta_text)

        result_prefix = workspace / "clusters"
        tmp_path = workspace / "tmp"
        command = [
            executable,
            "easy-cluster",
            str(fasta_path),
            str(result_prefix),
            str(tmp_path),
            "--min-seq-id",
            str(min_sequence_identity),
            "-c",
            str(coverage),
            "--cov-mode",
            "0",
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            raise RuntimeError(
                "MMseqs2 clustering failed:\n"
                + completed.stderr[-4000:]
                + completed.stdout[-1000:]
            )

        cluster_tsv = Path(f"{result_prefix}_cluster.tsv")
        pairs = pd.read_csv(
            cluster_tsv,
            sep="\t",
            header=None,
            names=["representative_token", "member_token"],
        )
        unknown = (set(pairs["representative_token"]) | set(pairs["member_token"])) - set(
            token_to_protein
        )
        if unknown:
            raise ValueError(f"MMseqs2 output contains unknown sequence IDs: {sorted(unknown)[:5]}")

        representative_tokens = sorted(pairs["representative_token"].unique())
        labels = {
            token: f"cluster_{index:05d}"
            for index, token in enumerate(representative_tokens, start=1)
        }
        mapping = pd.DataFrame(
            {
                "protein_id": pairs["member_token"].map(token_to_protein),
                "sequence_cluster": pairs["representative_token"].map(labels),
                "representative_protein_id": pairs["representative_token"].map(token_to_protein),
            }
        )
        mapping = validate_cluster_mapping(mapping)
        if len(mapping) != len(proteins):
            raise ValueError(
                f"MMseqs2 mapped {len(mapping)} of {len(proteins)} input protein sequences"
            )

        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        mapping.to_csv(output_csv, index=False)
        metadata = {
            "method": "MMseqs2 easy-cluster",
            "min_sequence_identity": min_sequence_identity,
            "coverage": coverage,
            "coverage_mode": 0,
            "proteins": int(len(mapping)),
            "clusters": int(mapping["sequence_cluster"].nunique()),
            "command": command,
        }
        output_csv.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2))
        return mapping
    finally:
        if temporary is not None:
            temporary.cleanup()


def compare_grouping_strategies(
    feature_csv: str | Path,
    cluster_mapping_csv: str | Path,
    results_dir: str | Path,
    model_names: list[str] | None = None,
    config_path: str | Path = "configs/default.yaml",
) -> pd.DataFrame:
    """Compare protein-grouped and homology-cluster-grouped CV."""
    from cysmutml.models.train import evaluate_models

    features = pd.read_csv(feature_csv, low_memory=False)
    mapping = pd.read_csv(cluster_mapping_csv)
    attached = attach_sequence_clusters(features, mapping)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    all_metrics = []
    strategies = [
        ("protein_grouped", "protein_id"),
        ("homology_clustered", "sequence_cluster"),
    ]
    with tempfile.TemporaryDirectory(prefix="cysmutml_homology_") as temporary:
        attached_csv = Path(temporary) / "features_with_clusters.csv"
        attached.to_csv(attached_csv, index=False)
        for strategy, group_column in strategies:
            metrics, _, _ = evaluate_models(
                attached_csv,
                results_dir / strategy,
                config_path=config_path,
                include_structural=False,
                model_names=model_names,
                group_column=group_column,
            )
            metrics.insert(0, "split_strategy", strategy)
            all_metrics.append(metrics)

    combined = pd.concat(all_metrics, ignore_index=True)
    combined.to_csv(results_dir / "split_comparison_fold_metrics.csv", index=False)
    summary = (
        combined.groupby(["split_strategy", "model"])[
            ["mae", "rmse", "r2", "pearson", "spearman"]
        ]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary.to_csv(results_dir / "split_comparison_summary.csv", index=False)
    cluster_sizes = mapping.groupby("sequence_cluster")["protein_id"].nunique()
    audit = {
        "proteins": int(mapping["protein_id"].nunique()),
        "sequence_clusters": int(mapping["sequence_cluster"].nunique()),
        "largest_cluster": int(cluster_sizes.max()),
        "median_cluster_size": float(cluster_sizes.median()),
        "feature_rows": int(len(attached)),
    }
    (results_dir / "cluster_audit.json").write_text(json.dumps(audit, indent=2))
    return summary

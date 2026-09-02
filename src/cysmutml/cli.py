"""Command-line interface for CysMutML."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cysmutml.data.audit import (
    audit_ddg_distribution,
    audit_duplicates_and_aggregate,
    summarize_structure_candidates,
)
from cysmutml.data.fireprotdb import download_fireprotdb_csv, prepare_data
from cysmutml.evaluation.ablation import run_ablation
from cysmutml.evaluation.error_analysis import write_error_analysis
from cysmutml.evaluation.homology import (
    build_mmseqs_cluster_map,
    compare_grouping_strategies,
)
from cysmutml.evaluation.structural_ablation import run_structural_ablation
from cysmutml.features.build import build_feature_table
from cysmutml.models.inference import predict_cys_mutations
from cysmutml.models.train import (
    evaluate_fast_baselines,
    evaluate_models,
    evaluate_physchem_model_comparison,
    train_final_model,
)
from cysmutml.ranking.engineering import rank_predictions
from cysmutml.structures.acquisition import (
    download_pdbs_for_candidates,
    download_uniprot_sequences_for_candidates,
)
from cysmutml.structures.structural_dataset import build_structural_mapping_and_features
from cysmutml.visualization.pymol import write_pymol_script, write_ranked_bfactor_pdb


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cysmutml")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prepare-data")
    p.add_argument("--raw", default="data/raw/fireprotdb.csv")
    p.add_argument("--output", default="data/processed/mutations.csv")
    p.add_argument("--download-fireprotdb", action="store_true")

    p = sub.add_parser("build-features")
    p.add_argument("--input", default="data/processed/mutations.csv")
    p.add_argument("--output", default="data/processed/features.csv")

    p = sub.add_parser("audit-data")
    p.add_argument("--processed", default="data/processed/fireprotdb_mutations.csv")
    p.add_argument("--raw", default="data/raw/fireprotdb.csv")

    p = sub.add_parser("download-structures")
    p.add_argument("--candidates", default="reports/structure_candidate_records.csv")
    p.add_argument("--structures-dir", default="data/structures")
    p.add_argument("--sequences-dir", default="data/sequences")

    p = sub.add_parser("build-structural-features")
    p.add_argument("--candidates", default="reports/structure_candidate_records.csv")
    p.add_argument("--min-identity", type=float, default=0.9)
    p.add_argument("--max-rows", type=int, default=None)

    p = sub.add_parser("train")
    p.add_argument("--features", default="data/processed/features.csv")
    p.add_argument("--model-output", default="models/cysmutml_model.joblib")
    p.add_argument("--metadata-output", default="models/model_metadata.json")
    p.add_argument("--model-name", default=None)

    p = sub.add_parser("evaluate")
    p.add_argument("--features", default="data/processed/features.csv")
    p.add_argument("--results-dir", default="results")
    p.add_argument(
        "--models",
        default=None,
        help="Comma-separated model names; default evaluates all configured models.",
    )

    p = sub.add_parser("evaluate-fast")
    p.add_argument("--features", default="data/processed/fireprotdb_features.csv")
    p.add_argument("--results-dir", default="results/fireprotdb_fast_baselines")
    p.add_argument("--save-oof", action="store_true")

    p = sub.add_parser("compare-physchem")
    p.add_argument("--features", default="data/processed/fireprotdb_aggregated_features.csv")
    p.add_argument("--results-dir", default="results/physchem_model_comparison")
    p.add_argument("--models", default="dummy_mean,ridge,hist_gradient_boosting")

    p = sub.add_parser("build-homology-clusters")
    p.add_argument("--input", default="data/processed/fireprotdb_mutations_aggregated.csv")
    p.add_argument("--output", default="data/processed/sequence_clusters.csv")
    p.add_argument("--min-sequence-identity", type=float, default=0.30)
    p.add_argument("--coverage", type=float, default=0.80)
    p.add_argument("--mmseqs-binary", default="mmseqs")
    p.add_argument("--work-dir", default=None)

    p = sub.add_parser("compare-grouping-strategies")
    p.add_argument("--features", default="data/processed/fireprotdb_aggregated_features.csv")
    p.add_argument("--clusters", default="data/processed/sequence_clusters.csv")
    p.add_argument("--results-dir", default="results/homology_validation")
    p.add_argument(
        "--models",
        default="dummy_mean,ridge,random_forest,hist_gradient_boosting",
    )
    p.add_argument("--target-proteins", type=int, default=150)
    p.add_argument("--random-seed", type=int, default=42)

    p = sub.add_parser("ablation")
    p.add_argument("--features", default="data/processed/features.csv")
    p.add_argument("--results-dir", default="results")

    p = sub.add_parser("error-analysis")
    p.add_argument(
        "--predictions",
        default="results/fireprotdb_fast_baselines/fast_baseline_out_of_fold_predictions.csv",
    )
    p.add_argument("--output-dir", default="results/fireprotdb_fast_baselines")
    p.add_argument("--model", default="ridge_fast")

    p = sub.add_parser("structural-ablation")
    p.add_argument(
        "--structural-features", default="data/processed/fireprotdb_structural_features.csv"
    )
    p.add_argument("--folds", default="data/processed/structural_cv_folds.csv")
    p.add_argument("--results-dir", default="results/structural_ablation")

    p = sub.add_parser("predict")
    p.add_argument("--pdb", required=True)
    p.add_argument("--chain", required=True)
    p.add_argument("--target", default="CYS", choices=["CYS"])
    p.add_argument("--model", default="models/cysmutml_model.joblib")
    p.add_argument("--output", default="results/prediction")
    p.add_argument("--protected-residues", default=None)
    p.add_argument("--config", default="configs/default.yaml")
    p.add_argument("--monocysteine-design", action="store_true")

    p = sub.add_parser("rank")
    p.add_argument("--predictions", required=True)
    p.add_argument("--output", required=True)

    p = sub.add_parser("visualize")
    p.add_argument("--ranking", required=True)
    p.add_argument("--pdb", required=True)
    p.add_argument("--output", required=True)

    args = parser.parse_args(argv)
    if args.command == "prepare-data":
        raw = Path(args.raw)
        if args.download_fireprotdb:
            download_fireprotdb_csv(raw)
        summary = prepare_data(raw, args.output)
        print(json.dumps(summary, indent=2))
    elif args.command == "build-features":
        df = build_feature_table(args.input, args.output)
        print(f"Wrote {len(df)} feature rows to {args.output}")
    elif args.command == "audit-data":
        duplicate_summary = audit_duplicates_and_aggregate(args.processed)
        ddg_summary = audit_ddg_distribution(args.processed)
        structure_summary = summarize_structure_candidates(args.raw)
        print(
            json.dumps(
                {
                    "duplicates": duplicate_summary,
                    "ddg": ddg_summary,
                    "structure_candidates": structure_summary,
                },
                indent=2,
            )
        )
    elif args.command == "download-structures":
        pdb_report = download_pdbs_for_candidates(
            args.candidates, args.structures_dir, "reports/structure_download_report.csv"
        )
        sequence_report = download_uniprot_sequences_for_candidates(
            args.candidates, args.sequences_dir, "reports/sequence_download_report.csv"
        )
        print(
            json.dumps(
                {
                    "pdbs": pdb_report["download_status"].value_counts().to_dict(),
                    "sequences": sequence_report["download_status"].value_counts().to_dict(),
                },
                indent=2,
            )
        )
    elif args.command == "build-structural-features":
        mapping, features = build_structural_mapping_and_features(
            args.candidates, min_identity=args.min_identity, max_rows=args.max_rows
        )
        print(
            json.dumps(
                {
                    "mapping_status": mapping["mapping_status"].value_counts().to_dict(),
                    "structural_feature_rows": int(len(features)),
                },
                indent=2,
            )
        )
    elif args.command == "evaluate":
        model_names = args.models.split(",") if args.models else None
        metrics, _, best = evaluate_models(args.features, args.results_dir, model_names=model_names)
        print(metrics.groupby("model")["mae"].mean().sort_values())
        print(f"Best by mean MAE: {best}")
    elif args.command == "evaluate-fast":
        metrics, cys = evaluate_fast_baselines(
            args.features, args.results_dir, save_oof=args.save_oof
        )
        print("Overall MAE by model:")
        print(metrics.groupby("model")["mae"].mean().sort_values())
        print("Cys-only MAE by model:")
        print(cys.groupby("model")["mae"].mean().sort_values())
    elif args.command == "compare-physchem":
        model_names = [name.strip() for name in args.models.split(",") if name.strip()]
        metrics, cys, _ = evaluate_physchem_model_comparison(
            args.features, args.results_dir, model_names=model_names
        )
        print("Overall MAE by model:")
        print(metrics.groupby("model")["mae"].mean().sort_values())
        print("Cys-only MAE by model:")
        print(cys.groupby("model")["mae"].mean().sort_values())
    elif args.command == "build-homology-clusters":
        mapping = build_mmseqs_cluster_map(
            args.input,
            args.output,
            min_sequence_identity=args.min_sequence_identity,
            coverage=args.coverage,
            mmseqs_binary=args.mmseqs_binary,
            work_dir=args.work_dir,
        )
        print(
            json.dumps(
                {
                    "proteins": int(mapping["protein_id"].nunique()),
                    "sequence_clusters": int(mapping["sequence_cluster"].nunique()),
                    "output": args.output,
                },
                indent=2,
            )
        )
    elif args.command == "compare-grouping-strategies":
        model_names = [name.strip() for name in args.models.split(",") if name.strip()]
        summary = compare_grouping_strategies(
            args.features,
            args.clusters,
            args.results_dir,
            model_names=model_names,
            target_proteins=args.target_proteins,
            random_seed=args.random_seed,
        )
        print(summary.to_string(index=False))
    elif args.command == "ablation":
        print(run_ablation(args.features, args.results_dir))
    elif args.command == "error-analysis":
        outputs = write_error_analysis(args.predictions, args.output_dir, args.model)
        print(json.dumps({key: str(value) for key, value in outputs.items()}, indent=2))
    elif args.command == "structural-ablation":
        outputs = run_structural_ablation(args.structural_features, args.folds, args.results_dir)
        print(
            json.dumps(
                {key: len(value) for key, value in outputs.items()},
                indent=2,
            )
        )
    elif args.command == "train":
        metadata = train_final_model(
            args.features, args.model_output, args.metadata_output, model_name=args.model_name
        )
        print(json.dumps(metadata, indent=2, default=str))
    elif args.command == "predict":
        _, warnings = predict_cys_mutations(
            args.pdb,
            args.chain,
            args.model,
            args.output,
            protected_residues=args.protected_residues,
            config_path=args.config,
            monocysteine_design=args.monocysteine_design,
        )
        ranking_path = Path(args.output) / "residue_ranking.csv"
        rank_predictions(
            Path(args.output) / "mutation_predictions.csv",
            ranking_path,
            config_path=args.config,
        )
        write_ranked_bfactor_pdb(
            ranking_path, args.pdb, Path(args.output) / "ranked_structure.pdb"
        )
        write_pymol_script(ranking_path, args.pdb, Path(args.output) / "visualize_rankings.pml")
        if warnings:
            print("Out-of-domain warnings:")
            print("\n".join(warnings))
    elif args.command == "rank":
        rank_predictions(args.predictions, args.output)
        print(f"Wrote ranking to {args.output}")
    elif args.command == "visualize":
        write_pymol_script(args.ranking, args.pdb, args.output)
        print(f"Wrote PyMOL script to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

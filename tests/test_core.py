from pathlib import Path

import joblib
import pandas as pd
from sklearn.model_selection import GroupKFold

from cysmutml.amino_acids import physicochemical_features
from cysmutml.data.audit import audit_duplicates_and_aggregate
from cysmutml.data.fireprotdb import normalize_fireprotdb_table
from cysmutml.evaluation.structural_ablation import make_structural_cv_folds
from cysmutml.features.build import build_feature_table
from cysmutml.models.inference import (
    _parse_protected_residues,
    generate_cys_feature_rows,
    predict_cys_mutations,
)
from cysmutml.models.train import train_final_model
from cysmutml.mutations import parse_mutation
from cysmutml.ranking.engineering import (
    accessibility_component_from_sasa,
    flexibility_component_from_proxy,
    lysine_environment_component_from_count,
    rank_predictions,
    rigidity_component_from_flexibility,
    stability_component_from_ddg,
)
from cysmutml.structures.features import chain_feature_rows
from cysmutml.structures.io import get_chain, parse_pdb
from cysmutml.structures.mapping import map_sequence_position_to_chain

ROOT = Path(__file__).resolve().parents[1]
PDB = ROOT / "examples" / "tiny_protein.pdb"


def test_mutation_parsing():
    mutation = parse_mutation("K12C")
    assert mutation.wt == "K"
    assert mutation.position == 12
    assert mutation.mut == "C"


def test_property_delta_direction():
    features = physicochemical_features("K", "C")
    assert features["delta_charge"] == -1
    assert features["mut_aa"] == "C"


def test_sign_convention_normalization():
    raw = pd.DataFrame({"protein": ["p"], "mutation": ["A1C"], "ddg": [1.5]})
    out, summary = normalize_fireprotdb_table(raw)
    assert summary["raw_records"] == 1
    assert out.loc[0, "destabilization_ddg_kcal_mol"] == 1.5


def test_structure_features():
    rows = chain_feature_rows(PDB, "A")
    assert len(rows) == 5
    assert rows[0]["ca_neighbors_6a"] > 0
    assert rows[0]["relative_sasa"] >= 0


def test_residue_mapping_and_wt_verification():
    chain = get_chain(parse_pdb(PDB), "A")
    mapped = map_sequence_position_to_chain("AKVEF", chain, 2, "K")
    assert mapped.status == "mapped"
    failed = map_sequence_position_to_chain("AKVEF", chain, 2, "A")
    assert failed.status == "failed"
    assert failed.failure_reason == "wt_mismatch"


def test_feature_schema_and_group_split(tmp_path):
    input_csv = ROOT / "data" / "raw" / "synthetic_fireprotdb_like.csv"
    features_csv = tmp_path / "features.csv"
    df = build_feature_table(input_csv, features_csv)
    assert {"delta_volume", "blosum62", "wt_hydrophobicity"}.issubset(df.columns)
    groups = df["protein_id"].astype(str)
    for train_idx, test_idx in GroupKFold(n_splits=3).split(
        df, df["destabilization_ddg_kcal_mol"], groups
    ):
        assert set(groups.iloc[train_idx]).isdisjoint(set(groups.iloc[test_idx]))


def test_deployed_model_feature_schema_if_present():
    model_path = ROOT / "models" / "cysmutml_model.joblib"
    if not model_path.exists():
        return
    artifact = joblib.load(model_path)
    assert artifact["categorical_features"] == ["wt_aa", "mut_aa"]
    assert artifact["numeric_features"] == [
        "wt_hydrophobicity",
        "wt_volume",
        "wt_mass",
        "wt_charge",
        "wt_polarity",
        "wt_aromatic",
        "mut_hydrophobicity",
        "mut_volume",
        "mut_mass",
        "mut_charge",
        "mut_polarity",
        "mut_aromatic",
        "delta_hydrophobicity",
        "delta_volume",
        "delta_mass",
        "delta_charge",
        "delta_polarity",
        "delta_aromatic",
        "blosum62",
    ]


def test_model_serialization(tmp_path):
    features_csv = tmp_path / "features.csv"
    build_feature_table(ROOT / "data" / "raw" / "synthetic_fireprotdb_like.csv", features_csv)
    model_path = tmp_path / "model.joblib"
    meta_path = tmp_path / "metadata.json"
    train_final_model(
        features_csv, model_path, meta_path, model_name="random_forest", include_structural=False
    )
    artifact = joblib.load(model_path)
    assert "model" in artifact
    assert meta_path.exists()


def test_inference_generation_and_ranking(tmp_path):
    features_csv = tmp_path / "features.csv"
    build_feature_table(ROOT / "data" / "raw" / "synthetic_fireprotdb_like.csv", features_csv)
    model_path = tmp_path / "model.joblib"
    train_final_model(
        features_csv,
        model_path,
        tmp_path / "metadata.json",
        model_name="random_forest",
        include_structural=False,
    )
    generated = generate_cys_feature_rows(PDB, "A")
    assert len(generated) == 5
    predictions, _ = predict_cys_mutations(PDB, "A", model_path, tmp_path)
    assert (tmp_path / "mutation_predictions.csv").exists()
    ranked = rank_predictions(
        tmp_path / "mutation_predictions.csv", tmp_path / "residue_ranking.csv"
    )
    assert ranked.iloc[0]["rank_engineering"] == 1
    assert "cys_suitability_score" in ranked
    assert "cys_site_suitability" in ranked
    assert "rigidification_potential" in ranked
    assert "final_engineering_score" in ranked
    assert "stability_component" in predictions
    assert "accessibility_component" in ranked
    assert "flexibility_component" in ranked
    assert "local_exposed_lys_count" in predictions
    assert "nearest_existing_cys_distance" in predictions


def test_duplicate_aggregation_uses_median(tmp_path):
    data = pd.DataFrame(
        {
            "protein_id": ["p1", "p1", "p2"],
            "mutation": ["A1C", "A1C", "K2C"],
            "wt_aa": ["A", "A", "K"],
            "position": [1, 1, 2],
            "mut_aa": ["C", "C", "C"],
            "destabilization_ddg_kcal_mol": [0.0, 2.0, 1.0],
            "exp_temperature": [25, 37, 25],
            "ph": [7.0, 8.0, 7.0],
            "method": ["Urea", "Urea", "CD"],
            "measure": ["CD", "CD", "CD"],
            "source_dataset": ["x", "x", "x"],
            "pdb_id": ["1abc", "1abc", "2abc"],
        }
    )
    input_csv = tmp_path / "processed.csv"
    data.to_csv(input_csv, index=False)
    audit_duplicates_and_aggregate(
        input_csv,
        tmp_path / "dup.csv",
        tmp_path / "dup.md",
        tmp_path / "agg.csv",
    )
    agg = pd.read_csv(tmp_path / "agg.csv")
    row = agg[agg["mutation"] == "A1C"].iloc[0]
    assert row["median_destabilization_ddg"] == 1.0
    assert row["n_measurements"] == 2


def test_structural_cv_folds_have_no_protein_overlap(tmp_path):
    df = pd.DataFrame(
        {
            "protein_id": ["p1", "p1", "p2", "p2", "p3", "p3"],
            "destabilization_ddg_kcal_mol": [0.1, 0.2, 1.0, 1.2, -0.1, -0.2],
        }
    )
    structural_csv = tmp_path / "struct.csv"
    folds_csv = tmp_path / "folds.csv"
    df.to_csv(structural_csv, index=False)
    folds = make_structural_cv_folds(structural_csv, folds_csv)
    for fold in folds["fold"].unique():
        train = set(folds.loc[folds["fold"] != fold, "protein_id"])
        test = set(folds.loc[folds["fold"] == fold, "protein_id"])
        assert train.isdisjoint(test)
    assert folds_csv.exists()


def test_hybrid_score_components_are_normalized():
    stability = stability_component_from_ddg(pd.Series([-1.0, 0.5, 2.0]))
    accessibility = accessibility_component_from_sasa(pd.Series([0.0, 0.5, 1.5]))
    rigidity = rigidity_component_from_flexibility(pd.Series([2.0, 1.0, 0.0]))
    flexibility = flexibility_component_from_proxy(pd.Series([2.0, 1.0, 0.0]))
    lys = lysine_environment_component_from_count(pd.Series([0, 3]))
    assert list(stability.round(3)) == [1.0, 0.5, 0.0]
    assert list(accessibility.round(3)) == [0.0, 0.5, 1.0]
    assert list(rigidity.round(3)) == [0.0, 0.5, 1.0]
    assert list(flexibility.round(3)) == [1.0, 0.5, 0.0]
    assert list(lys.round(3)) == [0.0, 0.632]


def test_rigidity_missing_values_use_neutral_fallback():
    neutral = rigidity_component_from_flexibility(pd.Series([float("nan")]))
    assert list(neutral) == [0.5]


def test_malformed_protected_residues_error():
    try:
        _parse_protected_residues("A45")
    except ValueError as exc:
        assert "CHAIN:RESIDUE" in str(exc)
    else:
        raise AssertionError("Malformed protected residue string should fail")


def test_ranking_score_reconstruction_and_penalties(tmp_path):
    predictions = pd.DataFrame(
        {
            "chain": ["A", "A"],
            "residue_number": ["1", "2"],
            "wt_aa": ["K", "L"],
            "mutation": ["K1C", "L2C"],
            "predicted_destabilization_ddg": [0.0, 0.0],
            "relative_sasa": [1.0, 0.0],
            "local_flexibility_proxy": [0.0, 1.0],
            "flexibility_method": ["BFACTOR", "BFACTOR"],
            "local_exposed_lys_count": [3, 0],
            "existing_cys_count_10A": [1, 0],
            "distance_to_nearest_protected": [1.0, 20.0],
            "distance_to_existing_cys": [1.0, 20.0],
        }
    )
    csv = tmp_path / "predictions.csv"
    out = tmp_path / "ranking.csv"
    predictions.to_csv(csv, index=False)
    ranked = rank_predictions(csv, out)
    first = ranked.iloc[0]
    site = (
        0.60 * first["stability_component"]
        + 0.35 * first["accessibility_component"]
        - 0.10 * first["existing_cys_penalty"]
        - 0.15 * first["protected_site_penalty"]
    )
    rigidification = (
        0.35 * first["flexibility_component"]
        + 0.40 * first["lysine_environment_component"]
        + 0.25 * first["accessibility_component"]
        - 0.05 * first["existing_cys_penalty"]
        - 0.10 * first["protected_site_penalty"]
    )
    final = 0.60 * site + 0.40 * rigidification
    assert abs(first["cys_site_suitability"] - site) < 1e-9
    assert abs(first["rigidification_potential"] - rigidification) < 1e-9
    assert abs(first["final_engineering_score"] - final) < 1e-9
    assert ranked["existing_cys_proximity_warning"].any()
    assert out.exists()


def test_local_lys_and_cys_environment_columns():
    generated = generate_cys_feature_rows(PDB, "A")
    assert "local_lys_count" in generated
    assert "local_exposed_lys_count" in generated
    assert "existing_cys_count_8A" in generated
    assert "existing_cys_count_10A" in generated
    assert (generated["lysine_radius_angstrom"] == 20.0).all()

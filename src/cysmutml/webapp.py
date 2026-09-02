"""Streamlit portfolio application for CysMutML."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio.PDB import PDBParser

from cysmutml.models.inference import predict_cys_mutations
from cysmutml.ranking.engineering import rank_predictions
from cysmutml.visualization.pymol import write_pymol_script, write_ranked_bfactor_pdb

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_PDB = ROOT / "examples" / "real_case" / "1csp.pdb"
MODEL_PATH = ROOT / "models" / "cysmutml_model.joblib"
CONFIG_PATH = ROOT / "configs" / "default.yaml"


@st.cache_data
def load_versioned_metrics() -> tuple[pd.DataFrame, pd.DataFrame]:
    result_dir = ROOT / "results" / "physchem_model_comparison"
    overall = pd.read_csv(result_dir / "regression_cv_metrics.csv")
    cys = pd.read_csv(result_dir / "cys_specific_metrics.csv")
    return overall, cys


def available_chains(pdb_path: str | Path) -> list[str]:
    structure = PDBParser(QUIET=True).get_structure("query", str(pdb_path))
    return sorted({chain.id for model in structure for chain in model})


def _metric_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    columns = ["mae", "rmse", "r2", "pearson", "spearman"]
    return metrics.groupby("model")[columns].mean().sort_values("mae").round(3)


def render_overview() -> None:
    st.subheader("A leakage-aware baseline for cysteine engineering")
    st.write(
        "CysMutML predicts mutation-associated destabilization from FireProtDB and "
        "combines that signal with transparent target-structure diagnostics."
    )
    first, second, third, fourth = st.columns(4)
    first.metric("Aggregated mutations", "352,005")
    second.metric("Protein groups", "542")
    third.metric("X→Cys records", "16,236")
    fourth.metric("Primary validation", "Group-aware")

    st.image(str(ROOT / "docs" / "figures" / "cysmutml_workflow.png"))
    st.info(
        "The learned stability model and the structural engineering heuristic are "
        "separate. The final ranking is not a calibrated probability."
    )


def render_benchmark() -> None:
    st.subheader("Model benchmark")
    overall, cys = load_versioned_metrics()
    st.caption(
        "Versioned full-data results use three-fold GroupKFold by protein. "
        "The homology-aware MVP is shown when its result artifacts are present."
    )

    left, right = st.columns(2)
    with left:
        st.markdown("**All mutations**")
        st.dataframe(_metric_summary(overall), use_container_width=True)
    with right:
        st.markdown("**X→Cys subset**")
        st.dataframe(_metric_summary(cys), use_container_width=True)

    mean_mae = overall.groupby("model")["mae"].mean().sort_values()
    st.bar_chart(mean_mae, horizontal=True, x_label="MAE (kcal/mol)")

    homology_folds = ROOT / "results" / "homology_validation" / (
        "split_comparison_fold_metrics.csv"
    )
    if homology_folds.exists():
        st.markdown("**150-protein homology-aware MVP**")
        folds = pd.read_csv(homology_folds)
        comparison = (
            folds.groupby(["split_strategy", "model"])[
                ["mae", "rmse", "r2", "spearman", "fit_seconds"]
            ]
            .mean()
            .round(3)
        )
        st.dataframe(comparison, use_container_width=True)
    else:
        st.warning(
            "The homology-split infrastructure is available, but numerical MVP "
            "artifacts have not yet been added to this checkout."
        )


def _run_prediction(pdb_path: Path, chain: str) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="cysmutml_app_") as temporary:
        output_dir = Path(temporary)
        _, warnings = predict_cys_mutations(
            pdb_path,
            chain,
            MODEL_PATH,
            output_dir,
            config_path=CONFIG_PATH,
        )
        ranking_path = output_dir / "residue_ranking.csv"
        ranking = rank_predictions(
            output_dir / "mutation_predictions.csv",
            ranking_path,
            config_path=CONFIG_PATH,
        )
        ranked_pdb = output_dir / "ranked_structure.pdb"
        pymol = output_dir / "visualize_rankings.pml"
        write_ranked_bfactor_pdb(ranking_path, pdb_path, ranked_pdb)
        write_pymol_script(ranking_path, pdb_path, pymol)
        return {
            "ranking": ranking,
            "ranking_csv": ranking_path.read_bytes(),
            "ranked_pdb": ranked_pdb.read_bytes(),
            "pymol": pymol.read_bytes(),
            "warnings": warnings,
        }


def render_prediction() -> None:
    st.subheader("Predict X→Cys candidates")
    source = st.radio("Structure", ["Built-in 1CSP example", "Upload PDB"], horizontal=True)
    uploaded = None
    if source == "Upload PDB":
        uploaded = st.file_uploader("PDB structure", type=["pdb"])

    temporary_path = None
    if uploaded is None and source == "Upload PDB":
        st.info("Upload a PDB file to continue.")
        return

    if uploaded is None:
        pdb_path = EXAMPLE_PDB
    else:
        suffix = Path(uploaded.name).suffix or ".pdb"
        handle = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        handle.write(uploaded.getvalue())
        handle.close()
        temporary_path = Path(handle.name)
        pdb_path = temporary_path

    try:
        chains = available_chains(pdb_path)
        if not chains:
            st.error("No protein chains were detected.")
            return
        chain = st.selectbox("Chain", chains)
        if st.button("Run CysMutML", type="primary"):
            with st.spinner("Predicting stability and calculating structural features..."):
                st.session_state["prediction"] = _run_prediction(pdb_path, chain)
    except Exception as error:
        st.error(f"Could not process this structure: {error}")
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    result = st.session_state.get("prediction")
    if not result:
        return

    ranking = result["ranking"]
    top_percent = st.slider("Candidates displayed (%)", 5, 100, 20, 5)
    count = max(1, round(len(ranking) * top_percent / 100))
    shown = ranking.nsmallest(count, "rank_engineering")
    display_columns = [
        "rank_engineering",
        "mutation",
        "predicted_destabilization_ddg",
        "relative_sasa",
        "cys_site_suitability",
        "rigidification_potential",
        "final_engineering_score",
    ]
    st.dataframe(shown[display_columns].round(3), use_container_width=True)
    chart = shown.set_index("mutation")[
        ["stability_component", "accessibility_component", "final_engineering_score"]
    ]
    st.bar_chart(chart)

    warnings = result["warnings"]
    if warnings:
        with st.expander("Out-of-domain warnings"):
            for warning in warnings:
                st.write(f"- {warning}")

    first, second, third = st.columns(3)
    first.download_button(
        "Download ranking CSV",
        result["ranking_csv"],
        file_name="cysmutml_ranking.csv",
        mime="text/csv",
    )
    second.download_button(
        "Download ranked PDB",
        result["ranked_pdb"],
        file_name="cysmutml_ranked.pdb",
        mime="chemical/x-pdb",
    )
    third.download_button(
        "Download PyMOL script",
        result["pymol"],
        file_name="visualize_cysmutml.pml",
        mime="text/plain",
    )


def render_methods() -> None:
    st.subheader("Methods and limitations")
    st.markdown(
        """
**Learned from FireProtDB**

- mutation-associated destabilization from amino-acid physicochemical descriptors;
- Ridge is deployed as the interpretable baseline.

**Calculated from the target PDB**

- relative SASA;
- B-factor-derived flexibility;
- local exposed Lys and native-Cys context;
- distances to user-protected residues.

**Not predicted**

CysMutML does not predict immobilization yield, retained activity, cysteine
reactivity, disulfide formation, or probability of experimental success.
"""
    )
    st.link_button(
        "Read the model card",
        "https://github.com/jjimenezgar/CysMutML/blob/main/MODEL_CARD.md",
    )


def main() -> None:
    st.set_page_config(page_title="CysMutML", page_icon="🧬", layout="wide")
    st.title("CysMutML")
    st.caption("Interpretable ML and structural bioinformatics for X→Cys prioritization")
    overview, benchmark, prediction, methods = st.tabs(
        ["Overview", "Model Benchmark", "Predict Cys Mutations", "Methods & Limitations"]
    )
    with overview:
        render_overview()
    with benchmark:
        render_benchmark()
    with prediction:
        render_prediction()
    with methods:
        render_methods()


if __name__ == "__main__":
    main()

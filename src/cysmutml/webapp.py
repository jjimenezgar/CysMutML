"""Streamlit portfolio application for CysMutML."""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from Bio.PDB import MMCIFParser, PDBParser

from cysmutml.models.inference import predict_cys_mutations
from cysmutml.ranking.engineering import rank_predictions
from cysmutml.visualization.pymol import write_pymol_script, write_ranked_bfactor_pdb

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_PDB = ROOT / "examples" / "real_case" / "1csp.pdb"
MODEL_PATH = ROOT / "models" / "cysmutml_model.joblib"
CONFIG_PATH = ROOT / "configs" / "default.yaml"
BRAND_IMAGE = ROOT / "docs" / "figures" / "cysmutml_github_cover.jpg"


def apply_style() -> None:
    st.markdown(
        """
        <style>
        .block-container { max-width: 1180px; padding-top: 2rem; }
        .hero { padding: 1.6rem 1.8rem; border-radius: 18px; background:
          linear-gradient(120deg, #071b3a 0%, #123d72 58%, #1d5da8 100%);
          color: white; margin-bottom: 1.2rem; }
        .hero h1 { margin: 0; font-size: 2.35rem; letter-spacing: -0.04em; }
        .hero p { margin: 0.45rem 0 0; color: #dbeafe; font-size: 1.05rem; }
        div[data-testid="stMetric"] { background: #f5f8fc; border: 1px solid #e3eaf3;
          padding: .7rem .9rem; border-radius: 12px; }
        .section-note { color: #526276; font-size: .92rem; }
        </style>
        """,
        unsafe_allow_html=True,
    )


@st.cache_data
def load_versioned_metrics() -> tuple[pd.DataFrame, pd.DataFrame]:
    result_dir = ROOT / "results" / "physchem_model_comparison"
    return (
        pd.read_csv(result_dir / "regression_cv_metrics.csv"),
        pd.read_csv(result_dir / "cys_specific_metrics.csv"),
    )


@st.cache_data(show_spinner=False)
def download_structure(source: str, identifier: str) -> tuple[bytes, str]:
    identifier = identifier.strip()
    if source == "RCSB / PDB ID":
        pdb_id = identifier.upper()
        if not re.fullmatch(r"[0-9A-Z]{4}", pdb_id):
            raise ValueError("Enter a valid four-character PDB code, for example 1CSP.")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        filename = f"{pdb_id}.pdb"
    else:
        accession = identifier.upper().split("-")[0]
        if not re.fullmatch(r"[A-Z0-9]{6,10}", accession):
            raise ValueError("Enter a valid UniProt accession, for example P0A7E1.")

        request = Request(
            f"https://rest.uniprot.org/uniprotkb/{accession}.txt",
            headers={"User-Agent": "CysMutML/1.2"},
        )
        try:
            with urlopen(request, timeout=30) as response:
                uniprot_text = response.read().decode("utf-8", errors="replace")
        except (HTTPError, URLError) as error:
            raise ValueError(f"UniProt accession not found: {accession}") from error

        pdb_ids = re.findall(r"^DR   PDB;\s*([0-9A-Za-z]{4});", uniprot_text, flags=re.MULTILINE)
        if pdb_ids:
            pdb_id = pdb_ids[0].upper()
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            filename = f"{accession}_{pdb_id}.pdb"
        else:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
            filename = f"AF-{accession}-F1-model_v4.pdb"

    try:
        with urlopen(Request(url, headers={"User-Agent": "CysMutML/1.2"}), timeout=60) as response:
            payload = response.read()
    except (HTTPError, URLError) as error:
        raise ValueError(f"Could not download structure from {url}") from error
    if not payload.startswith((b"HEADER", b"ATOM", b"MODEL", b"REMARK")):
        raise ValueError("The downloaded response is not a readable PDB structure.")
    return payload, filename


def available_chains(structure_path: str | Path) -> list[str]:
    path = Path(structure_path)
    if path.suffix.lower() in {".cif", ".mmcif"}:
        structure = MMCIFParser(QUIET=True).get_structure("query", str(path))
    else:
        structure = PDBParser(QUIET=True).get_structure("query", str(path))
    return sorted({chain.id for model in structure for chain in model})


def _metric_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    columns = ["mae", "rmse", "r2", "pearson", "spearman"]
    summary = metrics.groupby("model")[columns].mean().sort_values("mae").round(3)
    summary.index = summary.index.map(
        {
            "dummy_mean": "Baseline (mean)",
            "ridge": "Ridge",
            "random_forest": "Random forest",
            "hist_gradient_boosting": "Gradient boosting",
        }
    )
    summary.index.name = "Model"
    return summary.rename(
        columns={
            "mae": "MAE",
            "rmse": "RMSE",
            "r2": "R²",
            "pearson": "Pearson",
            "spearman": "Spearman",
        }
    )


def _humanize_ranking(ranking: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "rank_engineering",
        "mutation",
        "predicted_destabilization_ddg",
        "stability_score",
        "sasa_score",
        "flexibility_score",
        "lysine_boost",
        "existing_cys_penalty",
        "final_engineering_score",
    ]
    available = [column for column in columns if column in ranking.columns]
    shown = ranking[available].copy()
    shown = shown.rename(
        columns={
            "rank_engineering": "Priority",
            "mutation": "Candidate",
            "predicted_destabilization_ddg": "Predicted ΔΔG",
            "relative_sasa": "Relative exposure",
            "stability_score": "ML stability",
            "sasa_score": "Relative exposure",
            "flexibility_score": "Flexibility",
            "lysine_boost": "Nearby Lys boost",
            "existing_cys_penalty": "Nearby Cys penalty",
            "final_engineering_score": "Final priority",
        }
    )
    numeric = shown.select_dtypes(include="number").columns
    shown[numeric] = shown[numeric].round(3)
    return shown


def render_protein_viewer(
    structure_bytes: bytes,
    structure_format: str,
    chain: str,
    ranking: pd.DataFrame,
    top_n: int,
) -> None:
    """Render a lightweight 3Dmol.js view with the selected candidates highlighted."""
    residue_numbers = []
    for mutation in ranking.head(top_n)["mutation"].astype(str):
        match = re.search(r"(\d+)", mutation)
        if match:
            residue_numbers.append(int(match.group(1)))
    residue_numbers = sorted(set(residue_numbers))
    if not residue_numbers:
        st.info("No residue positions were available for the selected candidates.")
        return

    structure_text = structure_bytes.decode("utf-8", errors="replace")
    model_literal = json.dumps(structure_text)
    format_literal = json.dumps("cif" if structure_format in {"cif", "mmcif"} else "pdb")
    chain_literal = json.dumps(chain)
    residues_literal = json.dumps(residue_numbers)
    components.html(
        f"""
        <div id="cysmutml-viewer" style="width:100%;height:540px;position:relative;"></div>
        <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <script>
          const element = document.getElementById("cysmutml-viewer");
          const viewer = $3Dmol.createViewer(element, {{
            backgroundColor: "white",
            antialias: true
          }});
          viewer.addModel({model_literal}, {format_literal});
          viewer.setStyle({{}}, {{cartoon: {{color: "#cbd5e1"}}}});
          viewer.setStyle(
            {{chain: {chain_literal}, resi: {residues_literal}}},
            {{cartoon: {{color: "#f59e0b"}}, stick: {{radius: 0.22, color: "#dc2626"}}}}
          );
          viewer.zoomTo({{chain: {chain_literal}}});
          viewer.render();
          window.addEventListener("resize", () => viewer.resize(), false);
        </script>
        """,
        height=560,
        scrolling=False,
    )
    st.caption(
        f"Highlighted candidates: {', '.join(str(number) for number in residue_numbers)}. "
        "Gold residues are the selected positions; red sticks show their local environment."
    )


def render_overview() -> None:
    if BRAND_IMAGE.exists():
        st.image(str(BRAND_IMAGE), use_container_width=True)
    st.markdown(
        """
        <div class="hero">
          <p>Interpretable machine learning and structural analysis for X→Cys prioritisation.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.write(
        "CysMutML estimates mutation-associated destabilisation from FireProtDB and "
        "combines that signal with transparent, structure-derived ranking components."
    )
    first, second, third, fourth = st.columns(4)
    first.metric("Training rows", "352,005")
    second.metric("Protein groups", "542")
    third.metric("X→Cys rows", "16,236")
    fourth.metric("Validation", "Protein-aware")
    st.markdown(
        '<p class="section-note">The ML model and the structural ranking heuristic are deliberately separate. '
        "Neither is a calibrated probability of experimental success.</p>",
        unsafe_allow_html=True,
    )


def render_benchmark() -> None:
    st.subheader("Model benchmark")
    overall, cys = load_versioned_metrics()
    st.caption("Three-fold grouped validation on the full physicochemical dataset.")
    left, right = st.columns(2)
    with left:
        st.markdown("**All mutations**")
        st.dataframe(_metric_summary(overall), use_container_width=True)
    with right:
        st.markdown("**X→Cys subset**")
        st.dataframe(_metric_summary(cys), use_container_width=True)

    st.markdown("**Mean absolute error**")
    mean_mae = overall.groupby("model")["mae"].mean().sort_values()
    mean_mae.index = mean_mae.index.map(
        {
            "dummy_mean": "Baseline (mean)",
            "ridge": "Ridge",
            "random_forest": "Random forest",
            "hist_gradient_boosting": "Gradient boosting",
        }
    )
    st.bar_chart(mean_mae, horizontal=True, x_label="MAE (kcal/mol)")

    homology_folds = ROOT / "results" / "homology_validation" / "split_comparison_fold_metrics.csv"
    if homology_folds.exists():
        st.divider()
        st.subheader("Homology-aware MVP")
        st.caption("150 proteins, 5,634 rows, MMseqs2 at 30% identity / 80% coverage, seed 42.")
        folds = pd.read_csv(homology_folds)
        comparison = (
            folds.groupby(["split_strategy", "model"])[
                ["mae", "rmse", "r2", "spearman", "fit_seconds"]
            ]
            .mean()
            .round(3)
            .reset_index()
            .rename(
                columns={
                    "split_strategy": "Validation split",
                    "model": "Model",
                    "mae": "MAE",
                    "rmse": "RMSE",
                    "r2": "R²",
                    "spearman": "Spearman",
                    "fit_seconds": "Fit time (s)",
                }
            )
        )
        comparison["Validation split"] = comparison["Validation split"].map(
            {
                "protein_grouped": "Protein grouped",
                "homology_clustered": "Homology clustered",
            }
        )
        comparison["Model"] = comparison["Model"].map(
            {
                "dummy_mean": "Baseline (mean)",
                "ridge": "Ridge",
                "random_forest": "Random forest",
                "hist_gradient_boosting": "Gradient boosting",
            }
        )
        st.dataframe(comparison, use_container_width=True, hide_index=True)
        st.info("The homology-clustered split is intentionally stricter and exposes residual relatedness between proteins.")

        st.markdown("#### What the validation splits mean")
        st.caption(
            "**Protein grouped:** all mutations from one protein stay in the same fold, so the model "
            "is tested on proteins it did not see during training. **Homology clustered:** proteins "
            "with similar sequences are first grouped with MMseqs2 and the whole cluster stays in one "
            "fold. This is a stricter test of performance on less-related protein families."
        )

    st.markdown("#### How to read these metrics")
    st.caption(
        "MAE is the average absolute error in kcal/mol (lower is better). "
        "RMSE penalises larger errors more strongly. R² measures explained variance; "
        "Pearson and Spearman describe linear and rank correlation. Fit time is the "
        "average training time per fold."
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
            "structure_bytes": pdb_path.read_bytes(),
            "structure_format": pdb_path.suffix.lower().lstrip("."),
        }


def render_prediction() -> None:
    st.subheader("Predict X→Cys candidates")
    st.caption("Choose a structure source. Technical identifiers stay out of the result tables.")
    source = st.radio(
        "Structure source",
        ["Bundled 1CSP example", "Upload PDB/mmCIF", "RCSB / PDB ID", "UniProt accession"],
        horizontal=True,
    )

    temporary_path: Path | None = None
    source_label = "1CSP example"
    if source == "Upload PDB/mmCIF":
        uploaded = st.file_uploader("Upload a structure file", type=["pdb", "cif", "mmcif"])
        if uploaded is None:
            st.info("Upload a PDB or mmCIF file to continue.")
            return
        suffix = Path(uploaded.name).suffix or ".pdb"
        handle = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        handle.write(uploaded.getvalue())
        handle.close()
        temporary_path = Path(handle.name)
        pdb_path = temporary_path
        source_label = uploaded.name
    elif source in {"RCSB / PDB ID", "UniProt accession"}:
        placeholder = "1CSP" if source.startswith("RCSB") else "P0A7E1"
        identifier = st.text_input("Identifier", placeholder=placeholder).strip()
        if not identifier:
            st.info("Enter an identifier to continue.")
            return
        try:
            payload, filename = download_structure(source, identifier)
        except ValueError as error:
            st.error(str(error))
            return
        handle = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        handle.write(payload)
        handle.close()
        temporary_path = Path(handle.name)
        pdb_path = temporary_path
        source_label = filename
        st.caption(f"Loaded structure: {filename}")
    else:
        pdb_path = EXAMPLE_PDB

    try:
        chains = available_chains(pdb_path)
        if not chains:
            st.error("No protein chains were detected.")
            return
        chain = st.selectbox("Chain", chains, key="prediction_chain")
        if st.button("Run prediction", type="primary"):
            with st.spinner("Running the stability model and structural ranking..."):
                st.session_state["prediction"] = _run_prediction(pdb_path, chain)
                st.session_state["prediction_source"] = source_label
    except Exception as error:
        st.error(f"Could not process this structure: {error}")
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    result = st.session_state.get("prediction")
    if not result:
        return

    ranking = result["ranking"]
    st.success(f"Analysis complete: {st.session_state.get('prediction_source', source_label)}")
    top_percent = st.slider("Candidates to display", 5, 100, 20, 5)
    count = max(1, round(len(ranking) * top_percent / 100))
    shown = ranking.nsmallest(count, "rank_engineering")
    st.metric("Candidates shown", f"{count} of {len(ranking)}")
    st.dataframe(_humanize_ranking(shown), use_container_width=True, hide_index=True)

    st.markdown("#### What the columns mean")
    st.caption(
        "ML stability is the model prediction mapped to a 0–1 preference score. "
        "Relative exposure is the residue's relative SASA. Flexibility comes from the "
        "local B-factor signal. Nearby Lys boost rewards accessible lysines within the "
        "configured radius. Nearby Cys penalty discourages candidates close to existing "
        "cysteines. Final priority is the weighted sum of these signals, not a probability."
    )

    chart_columns = [
        column
        for column in [
            "stability_score",
            "sasa_score",
            "flexibility_score",
            "lysine_boost",
            "final_engineering_score",
        ]
        if column in shown
    ]
    if chart_columns:
        chart = shown.set_index("mutation")[chart_columns].rename(
            columns={
                "stability_score": "ML stability",
                "sasa_score": "Exposure",
                "flexibility_score": "Flexibility",
                "lysine_boost": "Lys boost",
                "final_engineering_score": "Final priority",
            }
        )
        st.bar_chart(chart)

    st.divider()
    st.subheader("3D structure view")
    viewer_top_n = st.slider(
        "Candidates highlighted",
        min_value=1,
        max_value=min(25, len(ranking)),
        value=min(10, len(ranking)),
        step=1,
        help="Highlights the highest-priority candidates on the selected chain.",
    )
    render_protein_viewer(
        result["structure_bytes"],
        result["structure_format"],
        chain,
        ranking,
        viewer_top_n,
    )

    warnings = result["warnings"]
    if warnings:
        with st.expander("Warnings and domain checks"):
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

        The deployed Ridge model uses physicochemical mutation descriptors, mutation
        deltas and BLOSUM62. It estimates mutation-associated destabilisation.

        **Calculated from the target structure**

        Relative exposure, B-factor-derived flexibility, local exposed-lysine context,
        existing-cysteine context. Protected residues, when supplied, are kept as an optional exclusion annotation and do not change the default MVP score.

        **Interpretation**

        These signals help prioritise candidates for inspection or experiment. They do
        not predict immobilisation yield, retained activity, cysteine reactivity,
        disulfide formation or probability of success.
        """
    )
    st.link_button(
        "Read the model card",
        "https://github.com/jjimenezgar/CysMutML/blob/portfolio-v1.1/MODEL_CARD.md",
    )


def main() -> None:
    st.set_page_config(page_title="CysMutML", page_icon="🧬", layout="wide")
    apply_style()
    overview, benchmark, prediction, methods = st.tabs(
        ["Overview", "Model benchmark", "Predict candidates", "Methods"]
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

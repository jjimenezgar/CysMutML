"""Build high-confidence structural mapping and feature datasets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.PDB import ShrakeRupley
from Bio.PDB.Chain import Chain

from cysmutml.amino_acids import MAX_ASA_TIEN_2013
from cysmutml.features.build import add_physicochemical_features
from cysmutml.mutations import parse_mutation
from cysmutml.structures.acquisition import read_fasta_sequence
from cysmutml.structures.features import (
    ca_coord,
    chain_residues,
    residue_one_letter,
)
from cysmutml.structures.io import parse_pdb


@dataclass(frozen=True)
class ChainMapping:
    chain_id: str
    mapped_residue_number: str | None
    insertion_code: str | None
    mapped_wt: str | None
    sequence_identity: float
    alignment_coverage: float
    status: str
    failure_reason: str | None


@dataclass(frozen=True)
class ChainAlignment:
    chain_id: str
    position_to_residue_index: dict[int, int]
    sequence_identity: float
    alignment_coverage: float


def _chain_sequence(chain: Chain) -> tuple[str, list]:
    residues = chain_residues(chain)
    return "".join(residue_one_letter(residue) for residue in residues), residues


def _align_chain(canonical_sequence: str, chain: Chain) -> ChainAlignment | None:
    chain_seq, residues = _chain_sequence(chain)
    if not chain_seq:
        return None
    start = canonical_sequence.find(chain_seq)
    if start >= 0:
        position_to_residue_index = {
            canonical_position: residue_index
            for residue_index, canonical_position in enumerate(
                range(start + 1, start + len(chain_seq) + 1)
            )
        }
        return ChainAlignment(
            chain.id,
            position_to_residue_index,
            1.0,
            len(chain_seq) / len(canonical_sequence) if canonical_sequence else 0.0,
        )
    return None
    alignment = pairwise2.align.globalms(
        canonical_sequence, chain_seq, 2, -1, -10, -0.5, one_alignment_only=True
    )[0]
    canonical_aln, chain_aln = alignment.seqA, alignment.seqB
    paired = [
        (a, b) for a, b in zip(canonical_aln, chain_aln, strict=True) if a != "-" and b != "-"
    ]
    matches = sum(1 for a, b in paired if a == b)
    identity = matches / len(paired) if paired else 0.0
    coverage = len(paired) / len(canonical_sequence) if canonical_sequence else 0.0
    canonical_pos = 0
    chain_pos = 0
    position_to_residue_index = {}
    for can_char, pdb_char in zip(canonical_aln, chain_aln, strict=True):
        if can_char != "-":
            canonical_pos += 1
        if pdb_char != "-":
            chain_pos += 1
        if can_char != "-" and pdb_char != "-":
            position_to_residue_index[canonical_pos] = chain_pos - 1
    return ChainAlignment(chain.id, position_to_residue_index, identity, coverage)


def _mapping_from_alignment(
    alignment: ChainAlignment,
    residues: list,
    dataset_position: int,
    dataset_wt: str,
) -> ChainMapping:
    if dataset_position not in alignment.position_to_residue_index:
        return ChainMapping(
            alignment.chain_id,
            None,
            None,
            None,
            alignment.sequence_identity,
            alignment.alignment_coverage,
            "residue_missing_in_structure",
            "gap_at_mutation_position",
        )
    residue = residues[alignment.position_to_residue_index[dataset_position]]
    mapped_wt = residue_one_letter(residue)
    insertion_code = residue.id[2].strip() or None
    status = "mapped_verified" if mapped_wt == dataset_wt else "wt_mismatch"
    reason = None if status == "mapped_verified" else "mapped_wt_differs_from_dataset_wt"
    return ChainMapping(
        alignment.chain_id,
        str(residue.id[1]),
        insertion_code,
        mapped_wt,
        alignment.sequence_identity,
        alignment.alignment_coverage,
        status,
        reason,
    )


def _select_mapping(
    canonical_sequence: str,
    pdb_path: Path,
    dataset_position: int,
    dataset_wt: str,
    dataset_chain: str | None = None,
    min_identity: float = 0.9,
) -> ChainMapping:
    structure = parse_pdb(pdb_path)
    model = next(structure.get_models())
    chains = [model[dataset_chain]] if dataset_chain and dataset_chain in model else list(model)
    mappings = []
    for chain in chains:
        alignment = _align_chain(canonical_sequence, chain)
        if alignment is None:
            continue
        _, residues = _chain_sequence(chain)
        mappings.append(_mapping_from_alignment(alignment, residues, dataset_position, dataset_wt))
    verified = [
        mapping
        for mapping in mappings
        if mapping.status == "mapped_verified" and mapping.sequence_identity >= min_identity
    ]
    if dataset_chain and dataset_chain not in model:
        return ChainMapping(
            dataset_chain, None, None, None, 0.0, 0.0, "chain_ambiguous", "chain_absent"
        )
    if len(verified) == 1:
        return verified[0]
    if len(verified) > 1:
        best = max(verified, key=lambda item: (item.sequence_identity, item.alignment_coverage))
        same_best = [
            item
            for item in verified
            if np.isclose(item.sequence_identity, best.sequence_identity)
            and np.isclose(item.alignment_coverage, best.alignment_coverage)
        ]
        if len(same_best) > 1 and dataset_chain is None:
            return ChainMapping(
                ",".join(item.chain_id for item in same_best),
                None,
                None,
                dataset_wt,
                best.sequence_identity,
                best.alignment_coverage,
                "chain_ambiguous",
                "multiple_verified_chains_with_equal_alignment",
            )
        return best
    if mappings:
        return max(mappings, key=lambda item: (item.sequence_identity, item.alignment_coverage))
    return ChainMapping(None, None, None, None, 0.0, 0.0, "other_failure", "no_chains")


def _select_from_aligned_chains(
    aligned_chains: list[tuple[ChainAlignment, list]],
    dataset_position: int,
    dataset_wt: str,
    min_identity: float,
) -> ChainMapping:
    mappings = [
        _mapping_from_alignment(alignment, residues, dataset_position, dataset_wt)
        for alignment, residues in aligned_chains
    ]
    verified = [
        mapping
        for mapping in mappings
        if mapping.status == "mapped_verified" and mapping.sequence_identity >= min_identity
    ]
    if len(verified) == 1:
        return verified[0]
    if len(verified) > 1:
        best = max(verified, key=lambda item: (item.sequence_identity, item.alignment_coverage))
        same_best = [
            item
            for item in verified
            if np.isclose(item.sequence_identity, best.sequence_identity)
            and np.isclose(item.alignment_coverage, best.alignment_coverage)
        ]
        if len(same_best) > 1:
            return ChainMapping(
                ",".join(item.chain_id for item in same_best),
                None,
                None,
                dataset_wt,
                best.sequence_identity,
                best.alignment_coverage,
                "chain_ambiguous",
                "multiple_verified_chains_with_equal_alignment",
            )
        return best
    if mappings:
        return max(mappings, key=lambda item: (item.sequence_identity, item.alignment_coverage))
    return ChainMapping(None, None, None, None, 0.0, 0.0, "other_failure", "no_chains")


def _chain_feature_map(pdb_path: Path, chain_id: str) -> dict[tuple[str, str | None], dict]:
    structure = parse_pdb(pdb_path)
    atom_count = sum(1 for _ in structure.get_atoms())
    if atom_count > 5_000:
        raise ValueError(f"structure_too_large_for_interactive_sasa:{atom_count}")
    ShrakeRupley().compute(structure, level="R")
    chain = next(structure.get_models())[chain_id]
    residues = chain_residues(chain)
    all_ca = np.asarray(
        [ca_coord(residue) for residue in residues if ca_coord(residue) is not None]
    )
    center = all_ca.mean(axis=0)
    max_center_distance = max(float(np.linalg.norm(point - center)) for point in all_ca) or 1.0
    chain_b = [float(atom.bfactor) for residue in residues for atom in residue.get_atoms()]
    chain_b_mean = float(np.mean(chain_b)) if chain_b else np.nan
    chain_b_std = float(np.std(chain_b)) if chain_b else np.nan
    out = {}
    for residue in residues:
        aa = residue_one_letter(residue)
        target_ca = ca_coord(residue)
        if target_ca is None:
            continue
        other_ca = np.asarray(
            [
                ca_coord(other)
                for other in residues
                if other.id != residue.id and ca_coord(other) is not None
            ]
        )
        distances = np.linalg.norm(other_ca - target_ca, axis=1) if len(other_ca) else np.array([])
        b_factors = [float(atom.bfactor) for atom in residue.get_atoms()]
        residue_b = float(np.mean(b_factors)) if b_factors else np.nan
        sasa = float(getattr(residue, "sasa", np.nan))
        insertion_code = residue.id[2].strip() or None
        out[(str(residue.id[1]), insertion_code)] = {
            "residue_sasa_abs": sasa,
            "residue_sasa_rel": min(sasa / MAX_ASA_TIEN_2013[aa], 1.5),
            "ca_contacts_6A": int(np.sum(distances <= 6.0)) if len(distances) else 0,
            "ca_contacts_8A": int(np.sum(distances <= 8.0)) if len(distances) else 0,
            "ca_contacts_10A": int(np.sum(distances <= 10.0)) if len(distances) else 0,
            "normalized_ca_distance_to_structure_center": float(
                np.linalg.norm(target_ca - center) / max_center_distance
            ),
            "mean_residue_b_factor": residue_b,
            "chain_normalized_b_factor": (
                (residue_b - chain_b_mean) / chain_b_std if chain_b_std > 0 else 0.0
            ),
            "secondary_structure": "unknown",
        }
    return out


def build_structural_mapping_and_features(
    candidate_csv: str | Path = "reports/structure_candidate_records.csv",
    structure_report_csv: str | Path = "reports/structure_download_report.csv",
    sequence_report_csv: str | Path = "reports/sequence_download_report.csv",
    mapping_csv: str | Path = "reports/residue_mapping_report.csv",
    failures_csv: str | Path = "reports/residue_mapping_failures.csv",
    features_csv: str | Path = "data/processed/fireprotdb_structural_features.csv",
    aggregated_features_csv: str
    | Path = "data/processed/fireprotdb_structural_features_aggregated.csv",
    min_identity: float = 0.9,
    max_rows: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(candidate_csv, low_memory=False)
    if max_rows is not None:
        candidates = candidates.head(max_rows).copy()
    structure_report = pd.read_csv(structure_report_csv, low_memory=False)
    sequence_report = pd.read_csv(sequence_report_csv, low_memory=False)
    pdb_paths = {
        row.pdb_id: Path(row.local_path)
        for row in structure_report.itertuples()
        if row.download_status in {"cached_valid", "downloaded_valid"} and pd.notna(row.local_path)
    }
    sequence_paths = {
        row.uniprot_id: Path(row.local_path)
        for row in sequence_report.itertuples()
        if row.download_status in {"cached", "downloaded"} and pd.notna(row.local_path)
    }
    sequence_cache = {
        accession: read_fasta_sequence(path) for accession, path in sequence_paths.items()
    }
    alignment_cache: dict[tuple[str, str], list[tuple[ChainAlignment, list]]] = {}
    chain_feature_cache: dict[tuple[str, str], dict[tuple[str, str | None], dict]] = {}
    mapping_rows = []
    feature_rows = []
    for idx, row in candidates.reset_index(drop=True).iterrows():
        try:
            mutation = parse_mutation(row["mutation"])
        except ValueError:
            mapping_rows.append({"row_id": idx, "mapping_status": "noncanonical_residue"})
            continue
        pdb_id = str(row["pdb_id"]).lower()
        uniprot_id = str(row["uniprot_id"])
        if pdb_id not in pdb_paths:
            status = "structure_missing"
            mapping = None
        elif uniprot_id not in sequence_cache:
            status = "sequence_unavailable"
            mapping = None
        else:
            cache_key = (uniprot_id, pdb_id)
            if cache_key not in alignment_cache:
                structure = parse_pdb(pdb_paths[pdb_id])
                aligned = []
                for chain in next(structure.get_models()):
                    alignment = _align_chain(sequence_cache[uniprot_id], chain)
                    if alignment is None:
                        continue
                    _, residues = _chain_sequence(chain)
                    aligned.append((alignment, residues))
                alignment_cache[cache_key] = aligned
            mapping = _select_from_aligned_chains(
                alignment_cache[cache_key], mutation.position, mutation.wt, min_identity
            )
            status = mapping.status
        mapping_row = {
            "row_id": idx,
            "execution_subset_max_rows": max_rows,
            "protein_id": row["protein_id"],
            "pdb_id": pdb_id,
            "dataset_chain": row.get("chain") if "chain" in row else None,
            "selected_chain": mapping.chain_id if mapping else None,
            "original_mutation": row["mutation"],
            "dataset_position": mutation.position,
            "dataset_wt": mutation.wt,
            "dataset_mut": mutation.mut,
            "mapped_pdb_resseq": mapping.mapped_residue_number if mapping else None,
            "insertion_code": mapping.insertion_code if mapping else None,
            "mapped_wt": mapping.mapped_wt if mapping else None,
            "sequence_identity": mapping.sequence_identity if mapping else None,
            "alignment_coverage": mapping.alignment_coverage if mapping else None,
            "mapping_status": status,
            "failure_reason": mapping.failure_reason if mapping else status,
        }
        mapping_rows.append(mapping_row)
        if status != "mapped_verified":
            continue
        chain_key = (pdb_id, mapping.chain_id)
        try:
            if chain_key not in chain_feature_cache:
                chain_feature_cache[chain_key] = _chain_feature_map(
                    pdb_paths[pdb_id], mapping.chain_id
                )
            structural = chain_feature_cache[chain_key][
                (mapping.mapped_residue_number, mapping.insertion_code)
            ]
        except Exception as exc:
            mapping_rows[-1]["mapping_status"] = "other_failure"
            mapping_rows[-1]["failure_reason"] = f"structural_feature_failure:{exc}"
            continue
        feature_rows.append(
            {
                **row.to_dict(),
                **mapping_row,
                **structural,
                "destabilization_ddg_kcal_mol": row["destabilization_ddg_kcal_mol"],
                "wt_aa": mutation.wt,
                "position": mutation.position,
                "mut_aa": mutation.mut,
            }
        )
    mapping_df = pd.DataFrame(mapping_rows)
    feature_df = (
        add_physicochemical_features(pd.DataFrame(feature_rows)) if feature_rows else pd.DataFrame()
    )
    Path(mapping_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(features_csv).parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(mapping_csv, index=False)
    mapping_df[mapping_df["mapping_status"] != "mapped_verified"].to_csv(failures_csv, index=False)
    feature_df.to_csv(features_csv, index=False)
    if not feature_df.empty:
        key = ["protein_id", "original_mutation", "pdb_id", "selected_chain", "mapped_pdb_resseq"]
        aggregated = (
            feature_df.groupby(key, dropna=False)
            .agg(
                median_destabilization_ddg=("destabilization_ddg_kcal_mol", "median"),
                mean_destabilization_ddg=("destabilization_ddg_kcal_mol", "mean"),
                n_measurements=("destabilization_ddg_kcal_mol", "size"),
            )
            .reset_index()
        )
        aggregated.to_csv(aggregated_features_csv, index=False)
    else:
        pd.DataFrame().to_csv(aggregated_features_csv, index=False)
    return mapping_df, feature_df

"""Sequence-to-structure residue mapping."""

from __future__ import annotations

from dataclasses import dataclass

from Bio import pairwise2
from Bio.PDB.Chain import Chain

from cysmutml.structures.features import chain_residues, residue_one_letter


@dataclass(frozen=True)
class MappingResult:
    dataset_position: int
    mapped_residue_number: str | None
    dataset_wt: str
    mapped_wt: str | None
    status: str
    failure_reason: str | None = None


def chain_sequence_and_residues(chain: Chain) -> tuple[str, list]:
    residues = chain_residues(chain)
    return "".join(residue_one_letter(res) for res in residues), residues


def map_sequence_position_to_chain(
    canonical_sequence: str,
    chain: Chain,
    position_1based: int,
    dataset_wt: str,
) -> MappingResult:
    chain_seq, residues = chain_sequence_and_residues(chain)
    if position_1based < 1 or position_1based > len(canonical_sequence):
        return MappingResult(
            position_1based, None, dataset_wt, None, "failed", "position_out_of_range"
        )

    alignment = pairwise2.align.globalms(
        canonical_sequence, chain_seq, 2, -1, -10, -0.5, one_alignment_only=True
    )[0]
    canonical_aln, chain_aln = alignment.seqA, alignment.seqB
    canonical_pos = 0
    chain_pos = 0
    for can_char, pdb_char in zip(canonical_aln, chain_aln, strict=True):
        if can_char != "-":
            canonical_pos += 1
        if pdb_char != "-":
            chain_pos += 1
        if canonical_pos == position_1based and can_char != "-":
            if pdb_char == "-":
                return MappingResult(
                    position_1based, None, dataset_wt, None, "failed", "gap_in_structure"
                )
            residue = residues[chain_pos - 1]
            mapped_wt = residue_one_letter(residue)
            residue_number = f"{residue.id[1]}{residue.id[2].strip()}"
            if mapped_wt != dataset_wt:
                return MappingResult(
                    position_1based, residue_number, dataset_wt, mapped_wt, "failed", "wt_mismatch"
                )
            return MappingResult(position_1based, residue_number, dataset_wt, mapped_wt, "mapped")
    return MappingResult(position_1based, None, dataset_wt, None, "failed", "alignment_failed")

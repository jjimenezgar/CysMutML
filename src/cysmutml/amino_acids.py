"""Canonical amino-acid properties used by CysMutML."""

from __future__ import annotations

from dataclasses import dataclass

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

THREE_TO_ONE = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

# Kyte-Doolittle hydrophobicity, approximate residue volume, molecular weight,
# and simple biochemical categories from standard amino-acid descriptor tables.
AA_PROPERTIES: dict[str, dict[str, float | int | str]] = {
    "A": {
        "hydrophobicity": 1.8,
        "volume": 88.6,
        "mass": 89.09,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "C": {
        "hydrophobicity": 2.5,
        "volume": 108.5,
        "mass": 121.16,
        "charge": 0,
        "polarity": 1,
        "aromatic": 0,
    },
    "D": {
        "hydrophobicity": -3.5,
        "volume": 111.1,
        "mass": 133.10,
        "charge": -1,
        "polarity": 1,
        "aromatic": 0,
    },
    "E": {
        "hydrophobicity": -3.5,
        "volume": 138.4,
        "mass": 147.13,
        "charge": -1,
        "polarity": 1,
        "aromatic": 0,
    },
    "F": {
        "hydrophobicity": 2.8,
        "volume": 189.9,
        "mass": 165.19,
        "charge": 0,
        "polarity": 0,
        "aromatic": 1,
    },
    "G": {
        "hydrophobicity": -0.4,
        "volume": 60.1,
        "mass": 75.07,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "H": {
        "hydrophobicity": -3.2,
        "volume": 153.2,
        "mass": 155.16,
        "charge": 0.1,
        "polarity": 1,
        "aromatic": 1,
    },
    "I": {
        "hydrophobicity": 4.5,
        "volume": 166.7,
        "mass": 131.17,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "K": {
        "hydrophobicity": -3.9,
        "volume": 168.6,
        "mass": 146.19,
        "charge": 1,
        "polarity": 1,
        "aromatic": 0,
    },
    "L": {
        "hydrophobicity": 3.8,
        "volume": 166.7,
        "mass": 131.17,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "M": {
        "hydrophobicity": 1.9,
        "volume": 162.9,
        "mass": 149.21,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "N": {
        "hydrophobicity": -3.5,
        "volume": 114.1,
        "mass": 132.12,
        "charge": 0,
        "polarity": 1,
        "aromatic": 0,
    },
    "P": {
        "hydrophobicity": -1.6,
        "volume": 112.7,
        "mass": 115.13,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "Q": {
        "hydrophobicity": -3.5,
        "volume": 143.8,
        "mass": 146.15,
        "charge": 0,
        "polarity": 1,
        "aromatic": 0,
    },
    "R": {
        "hydrophobicity": -4.5,
        "volume": 173.4,
        "mass": 174.20,
        "charge": 1,
        "polarity": 1,
        "aromatic": 0,
    },
    "S": {
        "hydrophobicity": -0.8,
        "volume": 89.0,
        "mass": 105.09,
        "charge": 0,
        "polarity": 1,
        "aromatic": 0,
    },
    "T": {
        "hydrophobicity": -0.7,
        "volume": 116.1,
        "mass": 119.12,
        "charge": 0,
        "polarity": 1,
        "aromatic": 0,
    },
    "V": {
        "hydrophobicity": 4.2,
        "volume": 140.0,
        "mass": 117.15,
        "charge": 0,
        "polarity": 0,
        "aromatic": 0,
    },
    "W": {
        "hydrophobicity": -0.9,
        "volume": 227.8,
        "mass": 204.23,
        "charge": 0,
        "polarity": 0,
        "aromatic": 1,
    },
    "Y": {
        "hydrophobicity": -1.3,
        "volume": 193.6,
        "mass": 181.19,
        "charge": 0,
        "polarity": 1,
        "aromatic": 1,
    },
}

# Maximum ASA values from Tien et al. 2013, used for relative SASA normalization.
MAX_ASA_TIEN_2013 = {
    "A": 129.0,
    "C": 167.0,
    "D": 193.0,
    "E": 223.0,
    "F": 240.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "K": 236.0,
    "L": 201.0,
    "M": 224.0,
    "N": 195.0,
    "P": 159.0,
    "Q": 225.0,
    "R": 274.0,
    "S": 155.0,
    "T": 172.0,
    "V": 174.0,
    "W": 285.0,
    "Y": 263.0,
}

BLOSUM62: dict[tuple[str, str], int] = {}
_BLOSUM_ROWS = {
    "A": "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0",
    "C": "-1 9 -3 -4 -2 -3 -3 -3 -1 -1 -1 -3 -1 -3 -3 -1 -1 -2 -2 -1",
    "D": "-2 -3 6 2 -3 -1 -1 -1 -3 -4 -3 1 -1 0 -2 0 -1 -4 -3 -3",
    "E": "-2 -4 2 5 -3 -2 0 -1 -3 -3 -2 0 -1 2 0 0 -1 -3 -2 -2",
    "F": "0 -2 -3 -3 6 -3 -1 -2 -1 0 0 -3 -2 -3 -4 -2 -2 1 3 -1",
    "G": "-1 -3 -1 -2 -3 6 -2 -4 -4 -4 -3 0 -2 -2 -2 0 -2 -2 -3 -3",
    "H": "-1 -3 -1 0 -1 -2 8 -3 -3 -3 -2 1 -2 0 0 -1 -2 -2 2 -3",
    "I": "0 -3 -1 -1 -2 -4 -3 4 -3 2 1 -3 -3 -3 -3 -2 -1 -3 -1 3",
    "K": "-2 -1 -3 -3 -1 -4 -3 2 4 -2 -2 -3 -3 -3 -3 -2 -1 -3 -1 1",
    "L": "-1 -1 -4 -3 0 -4 -3 2 -2 4 2 -3 -3 -2 -2 -2 -1 -2 -1 1",
    "M": "-1 -1 -3 -2 0 -3 -2 1 -2 2 5 -2 -2 0 -1 -1 -1 -1 -1 1",
    "N": "-2 -3 1 0 -3 0 1 -3 -3 -3 -2 6 -2 0 0 1 0 -4 -2 -3",
    "P": "-1 -1 -1 -1 -2 -2 -2 -3 -3 -3 -2 -2 7 -1 -2 -1 -1 -4 -3 -2",
    "Q": "-2 -3 0 2 -3 -2 0 -3 -3 -2 0 0 -1 5 1 0 -1 -2 -1 -2",
    "R": "-1 -3 -2 0 -4 -2 0 -3 -3 -2 -1 0 -2 1 5 -1 -1 -3 -2 -3",
    "S": "1 -1 0 0 -2 0 -1 -2 -2 -2 -1 1 -1 0 -1 4 1 -3 -2 -2",
    "T": "0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1 0 -1 -1 -1 1 5 -2 -2 0",
    "W": "-3 -2 -4 -3 1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 11 2 -3",
    "Y": "-2 -2 -3 -2 3 -3 2 -1 -1 -1 -1 -2 -3 -1 -2 -2 -2 2 7 -1",
    "V": "0 -1 -3 -2 -1 -3 -3 3 1 1 1 -3 -2 -2 -3 -2 0 -3 -1 4",
}
_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
for aa1, row in _BLOSUM_ROWS.items():
    for aa2, value in zip(_ORDER, row.split(), strict=True):
        BLOSUM62[(aa1, aa2)] = int(value)


@dataclass(frozen=True)
class Mutation:
    wt: str
    position: int
    mut: str

    @property
    def label(self) -> str:
        return f"{self.wt}{self.position}{self.mut}"


def is_canonical_aa(value: str) -> bool:
    return value in CANONICAL_AA


def physicochemical_features(wt: str, mut: str) -> dict[str, float | int | str]:
    if wt not in CANONICAL_AA or mut not in CANONICAL_AA:
        raise ValueError(f"Noncanonical amino acid in mutation {wt}->{mut}")
    features: dict[str, float | int | str] = {"wt_aa": wt, "mut_aa": mut}
    for prefix, aa in (("wt", wt), ("mut", mut)):
        for name, value in AA_PROPERTIES[aa].items():
            features[f"{prefix}_{name}"] = value
    for name in ("hydrophobicity", "volume", "mass", "charge", "polarity", "aromatic"):
        features[f"delta_{name}"] = float(AA_PROPERTIES[mut][name]) - float(AA_PROPERTIES[wt][name])
    features["blosum62"] = BLOSUM62[(wt, mut)]
    return features

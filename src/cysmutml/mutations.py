"""Mutation parsing and validation."""

from __future__ import annotations

import re

from cysmutml.amino_acids import CANONICAL_AA, Mutation

MUTATION_RE = re.compile(r"^\s*([A-Z])\s*([0-9]+)\s*([A-Z])\s*$")


def parse_mutation(value: str) -> Mutation:
    match = MUTATION_RE.match(str(value))
    if not match:
        raise ValueError(f"Not a single substitution: {value!r}")
    wt, position, mut = match.groups()
    if wt not in CANONICAL_AA or mut not in CANONICAL_AA:
        raise ValueError(f"Noncanonical amino acid in mutation: {value!r}")
    return Mutation(wt=wt, position=int(position), mut=mut)

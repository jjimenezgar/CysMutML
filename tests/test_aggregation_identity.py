import pandas as pd
import pytest

from cysmutml.data.audit import validate_aggregation_identity


def table(ids, proteins=None):
    return pd.DataFrame({
        "protein_id": proteins if proteins is not None else [None, None],
        "wt_aa": ["A", "A"], "position": [10, 10], "mut_aa": ["C", "C"],
        "uniprot_id": ids,
    })


KEY = ["protein_id", "wt_aa", "position", "mut_aa"]


def test_anonymous_same_identity_is_preserved():
    validate_aggregation_identity(table(["P1", "P1"]), KEY)


@pytest.mark.parametrize("proteins", [[None, None], ["same_name", "same_name"]])
def test_distinct_proteins_cannot_be_aggregated(proteins):
    with pytest.raises(ValueError, match="Conflicting"):
        validate_aggregation_identity(table(["P1", "P2"], proteins), KEY)


def test_unresolved_anonymous_duplicates_are_rejected():
    with pytest.raises(ValueError, match="lack shared identity"):
        validate_aggregation_identity(table(["P1", None]), KEY)

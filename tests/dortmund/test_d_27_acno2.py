import pytest

import ugropy as ug


# =============================================================================
# 27|ACNO2|[57]ACNO2
# =============================================================================

# Dortmund
trials = [
    # nitrobenzene
    ("C1=CC=C(C=C1)[N+](=O)[O-]", {"ACH": 5, "ACNO2": 1}, "smiles"),
    ("[O-][N+](=O)C1=CC=NC=C1", {"ACH": 2, "AC2H2N": 1, "ACNO2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_acno2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 25|ACCL|[53]ACCL
# =============================================================================

# Dortmund
trials = [
    # chlorobenzene
    ("C1=CC=C(C=C1)Cl", {"ACCL": 1, "ACH": 5}, "smiles"),
    # Impossible
    ("ClC1=NC=CC=C1", {}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_accl_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

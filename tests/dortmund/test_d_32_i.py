import pytest

import ugropy as ug


# =============================================================================
# 32|I|[63]I
# =============================================================================

# Dortmund
trials = [
    # 1-iodoethane
    ("CCI", {"CH3": 1, "CH2": 1, "I": 1}, "smiles"),
    # Iodobenzene
    ("C1=CC=C(C=C1)I", {"ACH": 5, "AC": 1, "I": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_i_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

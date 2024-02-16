import pytest

import ugropy as ug


# =============================================================================
# 33|BR|[64]BR
# =============================================================================

# Dortmund
trials = [
    # 1-bromoethane
    ("CCBr", {"CH3": 1, "CH2": 1, "BR": 1}, "smiles"),
    # Bromobenzene
    ("C1=CC=C(C=C1)Br", {"ACH": 5, "AC": 1, "BR": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_br_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

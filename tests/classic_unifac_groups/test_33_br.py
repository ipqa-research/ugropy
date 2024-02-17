import pytest

import ugropy as ug


# =============================================================================
# 33- Br Main group: Br
# =============================================================================

# UNIFAC
trials_unifac = [
    # 1-bromoethane
    ("CCBr", {"CH3": 1, "CH2": 1, "BR": 1}, "smiles"),
    # Bromobenzene
    ("C1=CC=C(C=C1)Br", {"ACH": 5, "AC": 1, "BR": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_br_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result

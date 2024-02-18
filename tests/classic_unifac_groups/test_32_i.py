import pytest

import ugropy as ug


# =============================================================================
# 32- I Main group: I
# =============================================================================

# UNIFAC
trials_unifac = [
    # 1-iodoethane
    ("CCI", {"CH3": 1, "CH2": 1, "I": 1}, "smiles"),
    # Iodobenzene
    ("C1=CC=C(C=C1)I", {"ACH": 5, "AC": 1, "I": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_i_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result

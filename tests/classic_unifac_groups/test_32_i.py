import pytest

from ugropy import get_groups, psrk, unifac


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
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_i_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_i_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

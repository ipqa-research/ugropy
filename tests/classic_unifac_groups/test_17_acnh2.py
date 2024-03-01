import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 17- ACNH2 Main group: ACNH2
# =============================================================================

# UNIFAC
trials_unifac = [
    # aniline
    ("C1=CC=C(C=C1)N", {"ACH": 5, "ACNH2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 28- CS2 Main group: CS2
# =============================================================================

# UNIFAC
trials_unifac = [
    # carbon disulfide
    ("C(=S)=S", {"CS2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cs2_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cs2_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

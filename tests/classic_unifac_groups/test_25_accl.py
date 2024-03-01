import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 25- ACCL Main group: ACCL
# =============================================================================

# UNIFAC
trials_unifac = [
    # chlorobenzene
    ("C1=CC=C(C=C1)Cl", {"ACCL": 1, "ACH": 5}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_accl_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_accl_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

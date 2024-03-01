import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 36- ACRY Main group: ACRY
# =============================================================================

# UNIFAC
trials_unifac = [
    # acrylonitrile
    ("C=CC#N", {"ACRY": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acry_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acry_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 19- CCN Main group: CH3CN, CH2CN
# =============================================================================

# UNIFAC
trials_unifac = [
    # acetonitrile
    ("CC#N", {"CH3CN": 1}, "smiles"),
    # propionitrile
    ("CCC#N", {"CH3": 1, "CH2CN": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccn_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccn_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

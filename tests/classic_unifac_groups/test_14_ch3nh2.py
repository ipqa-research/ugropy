import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 14- CH3NH2 Main group: CH3NH2, CH2NH2, CHNH2
# =============================================================================

# UNIFAC
trials_unifac = [
    # methylamine
    ("CN", {"CH3NH2": 1}, "smiles"),
    # isopropylamine
    ("CC(C)N", {"CH3": 2, "CHNH2": 1}, "smiles"),
    # propylamine
    ("CCCN", {"CH3": 1, "CH2": 1, "CH2NH2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 6- CH3OH Main group: CH3OH
# =============================================================================

# UNIFAC
# methanol
trials_unifac = [("CO", {"CH3OH": 1}, "smiles")]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_oh_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_oh_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 7- H2O Main group: H2O
# =============================================================================

# UNIFAC
trials_unifac = [("O", {"H2O": 1}, "smiles")]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_h2o_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_h2o_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

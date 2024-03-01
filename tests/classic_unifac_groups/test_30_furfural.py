import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 30- furfural Main group: furfural
# =============================================================================

# UNIFAC
trials_unifac = [
    # furfural
    ("C1=COC(=C1)C=O", {"FURFURAL": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_furfural_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_furfural_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

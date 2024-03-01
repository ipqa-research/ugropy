import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 35- DMSO Main group: DMSO
# =============================================================================

# UNIFAC
trials_unifac = [
    # dimethyl sulfoxide
    ("CS(=O)C", {"DMSO": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_alquine_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_alquine_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

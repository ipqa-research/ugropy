import pytest

import ugropy as ug


# =============================================================================
# 6- CH3OH Main group: CH3OH
# =============================================================================

# UNIFAC
# methanol
trials_unifac = [("CO", {"CH3OH": 1}, "smiles")]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_oh_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result

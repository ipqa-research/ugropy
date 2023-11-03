import pytest

import ugropy as ug


# =============================================================================
# 7- H2O Main group: H2O
# =============================================================================

# UNIFAC
trials_unifac = [("O", {"H2O": 1}, "smiles")]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_h2o_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

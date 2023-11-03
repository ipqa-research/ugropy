import pytest

import ugropy as ug


# =============================================================================
# 36- ACRY Main group: ACRY
# =============================================================================

# UNIFAC
trials_unifac = [
    # acrylonitrile
    ("C=CC#N", {"ACRY": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acry_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

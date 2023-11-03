import pytest

import ugropy as ug


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
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccn_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 47- OCCOH Main group: C2H2O2, C2H4O2
# =============================================================================

# UNIFAC
trials_unifac = [
    # 2-Ethoxyethanol
    ("CCOCCO", {"CH3": 1, "CH2": 1, "C2H5O2": 1}, "smiles"),
    # 2-Ethoxy-1-propanol
    ("CCOC(C)CO", {"CH3": 2, "CH2": 1, "C2H4O2": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_occoh_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

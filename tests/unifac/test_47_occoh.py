import ugropy as ug

import pytest


# =============================================================================
# 47- OCCOH Main group: C2H2O2, C2H4O2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("2-Ethoxyethanol", {"CH3": 1, "CH2": 1, "C2H5O2": 1}, "name"),
    ("2-Ethoxy-1-propanol", {"CH3": 2, "CH2": 1, "C2H4O2": 1}, "name"),
]

@pytest.mark.OCCOH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_occoh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
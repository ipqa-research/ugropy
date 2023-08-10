import ugropy as ug

import pytest


# =============================================================================
# 26- CNO2 Main group: CH3NO2, CH2NO2, CHNO2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("nitromethane", {"CH3NO2": 1}, "name"),
    ("1-nitropropane", {"CH3": 1, "CH2": 1, "CH2NO2": 1}, "name"),
    ("2-nitropropane", {"CH3": 2, "CHNO2": 1}, "name"),
]

@pytest.mark.CNO2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cno2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
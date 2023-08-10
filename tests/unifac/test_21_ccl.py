import ugropy as ug

import pytest


# =============================================================================
# 21- CCL Main group: CH2CL, CHCL, CCL
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1-chlorobutane", {"CH3": 1, "CH2": 2, "CH2CL": 1}, "name"),
    ("2-chloropropane", {"CH3": 2, "CHCL": 1}, "name"),
    ("2-chloro-2-methylpropane", {"CH3": 3, "CCL": 1}, "name"),
]

@pytest.mark.CCL
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
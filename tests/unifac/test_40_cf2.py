import ugropy as ug

import pytest


# =============================================================================
# 40- CF2 Main group: CF3, CF2, CF
# =============================================================================

# UNIFAC
trials_unifac = [
    (" Perfluorohexane", {"CF3": 2 , "CF2": 4}, "name"),
    ("Perfluoromethylcyclohexane", {"CF3": 1, "CF2": 5, "CF": 1}, "name"),
]

@pytest.mark.CF2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cf2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
import ugropy as ug

import pytest


# =============================================================================
# 12- HCOO Main group: HCOO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("phenyl formate", {"ACH": 5, "AC": 1, "HCOO": 1}),
    ("ethyl formate", {"HCOO": 1, "CH2": 1, "CH3": 1}),
    ("methyl formate", {"HCOO": 1, "CH3": 1})
]

@pytest.mark.HCOO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_cho_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
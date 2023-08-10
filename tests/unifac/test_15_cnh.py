import ugropy as ug

import pytest


# =============================================================================
# 15- CNH Main group: CH3NH, CH2NH, CHNH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("dimethylamine", {"CH3": 1, "CH3NH": 1}, "name"),
    ("diethylamine", {"CH3": 2, "CH2": 1, "CH2NH": 1}, "name"),
    ("diisopropylamine", {"CH3": 4, "CH": 1, "CHNH": 1}, "name"),
]

@pytest.mark.CNH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
import pytest

import ugropy as ug


# =============================================================================
# 48- CH2S Main group: CH3S, CH2S, CHS
# =============================================================================

# UNIFAC
trials_unifac = [
    ("Dimethylsulfide", {"CH3": 1, "CH3S": 1}, "name"),
    ("Diethylsulfide", {"CH3": 2, "CH2": 1, "CH2S": 1}, "name"),
    ("Isopropyl Sulfide", {"CH3": 4, "CH": 1, "CHS": 1}, "name"),
]


@pytest.mark.CH2S
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2s_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

import ugropy as ug

import pytest


# =============================================================================
# 29- CH3SH Main group: CH3SH, CH2SH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("methanethiol", {"CH3SH": 1}, "name"),
    ("ethanethiol", {"CH2SH": 1, "CH3": 1}, "name")
]

@pytest.mark.CH3SH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3sh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
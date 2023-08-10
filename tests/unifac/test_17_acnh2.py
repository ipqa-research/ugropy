import ugropy as ug

import pytest


# =============================================================================
# 17- ACNH2 Main group: ACNH2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("aniline", {"ACH": 5, "ACNH2": 1}, "name"),
]

@pytest.mark.ACNH2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
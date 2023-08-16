import ugropy as ug

import pytest


# =============================================================================
# 18- Pyridine Main group: C5H5N, C5H4N, C5H3N
# =============================================================================

# UNIFAC
trials_unifac = [
    #("pyridine", {"C5H5N": 1}, "name"),
    ("3-methylpyridine", {"C5H4N": 1, "CH3": 1}, "name"),
    #("2,3-methylpyridine", {"C5H3N": 1, "CH3": 2}, "name"),
]

@pytest.mark.pyridine
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
@pytest.mark.skip
def test_ch3nh2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
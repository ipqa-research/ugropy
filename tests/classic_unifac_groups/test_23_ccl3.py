import pytest

import ugropy as ug


# =============================================================================
# 23- CCL3 Main group: CHCL3, CCL3
# =============================================================================

# UNIFAC
trials_unifac = [
    ("chloroform", {"CHCL3": 1}, "name"),
    ("trichloroethane", {"CH3": 1, "CCL3": 1}, "name"),
]


@pytest.mark.CCL3
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl3_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

import pytest

import ugropy as ug


# =============================================================================
# 28- CS2 Main group: CS2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("carbon disulfide", {"CS2": 1}, "name"),
]


@pytest.mark.CS2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cs2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

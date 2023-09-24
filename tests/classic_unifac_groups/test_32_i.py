import pytest

import ugropy as ug


# =============================================================================
# 32- I Main group: I
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1-iodoethane", {"CH3": 1, "CH2": 1, "I": 1}, "name"),
    ("Iodobenzene", {"ACH": 5, "AC": 1, "I": 1}, "name"),
]


@pytest.mark.I
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_i_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

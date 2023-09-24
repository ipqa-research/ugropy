import pytest

import ugropy as ug


# =============================================================================
# 33- Br Main group: Br
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1-bromoethane", {"CH3": 1, "CH2": 1, "BR": 1}, "name"),
    ("Bromobenzene", {"ACH": 5, "AC": 1, "BR": 1}, "name"),
]


@pytest.mark.BR
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_br_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

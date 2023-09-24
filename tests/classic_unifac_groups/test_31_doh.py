import pytest

import ugropy as ug


# =============================================================================
# 31- DOH Main group: DOH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1,2-ethanediol", {"DOH": 1}, "name"),
]


@pytest.mark.DOH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_doh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

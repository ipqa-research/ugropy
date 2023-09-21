import pytest

import ugropy as ug


# =============================================================================
# 6- CH3OH Main group: CH3OH
# =============================================================================

# UNIFAC
trials_unifac = [("methanol", {"CH3OH": 1}, "name")]


@pytest.mark.CH3OH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_oh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

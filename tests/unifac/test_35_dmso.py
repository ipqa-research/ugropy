import pytest

import ugropy as ug


# =============================================================================
# 35- DMSO Main group: DMSO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("dimethyl sulfoxide", {"DMSO": 1}, "name"),
]


@pytest.mark.DMSO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_alquine_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

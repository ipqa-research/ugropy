import ugropy as ug

import pytest


# =============================================================================
# 30- furfural Main group: furfural
# =============================================================================

# UNIFAC
trials_unifac = [
    ("furfural", {"FURFURAL": 1}, "name"),
]

@pytest.mark.FURFURAL
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_furfural_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
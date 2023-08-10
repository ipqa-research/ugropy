import ugropy as ug

import pytest


# =============================================================================
# 7- H2O Main group: H2O
# =============================================================================

# UNIFAC
trials_unifac = [
    ("water", {"H2O": 1})
]

@pytest.mark.H2O
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_h2o_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
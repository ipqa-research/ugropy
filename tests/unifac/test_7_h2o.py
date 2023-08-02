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
def test_ch2_unifac(name, result):
    substance = ug.Substance(name)
    assert substance.unifac_groups == result
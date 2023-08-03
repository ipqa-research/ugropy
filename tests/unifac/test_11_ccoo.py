import ugropy as ug

import pytest


# =============================================================================
# 11- CCOO Main group (aldehyde): CH3COO, CH2COO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("Tert-butyl acetate", {"CH3COO": 1, "CH3": 3, "C": 1}),
    ("triacetin", {"CH3COO": 3, "CH2": 2, "CH": 1}),
    ("butyl propanoate", {"CH3": 2, "CH2": 3, "CH2COO": 1}),
    ("butyl acetate", {"CH3": 1, "CH2": 3, "CH3COO": 1})
]

@pytest.mark.CCOO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_cho_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
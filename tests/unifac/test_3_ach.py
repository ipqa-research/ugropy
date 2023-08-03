import ugropy as ug

import pytest


# =============================================================================
# 3- ACH Main group: ACH, AC
# =============================================================================

# UNIFAC
trials_unifac = [
    ("phenanthrene", {"ACH": 10, "AC": 4}),
    ("anthracene", {"ACH": 10, "AC": 4}),
    ("1,1-Diphenylethylene", {"ACH": 10, "AC": 2, "CH2=C": 1}),
    ("biphenyl", {"ACH": 10, "AC": 2}),
    ("Naphthalene", {"ACH": 8, "AC": 2}),
    ("benzene", {"ACH": 6}),
    ("styrene", {"AC": 1, "CH2=CH": 1, "ACH": 5})
]

@pytest.mark.ACH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_ach_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
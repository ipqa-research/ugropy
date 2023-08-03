import ugropy as ug

import pytest


# =============================================================================
# 9- CH2CO Main group: CH3CO, CH2CO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("cyclopropanone", {"CH2": 1, "CH2CO": 1}),
    ("(9Z)-Cycloheptadec-9-en-1-one", {"CH2": 13, "CH=CH": 1, "CH2CO": 1}),
    ("1-(4-tert-butyl-2,6-dimethylphenyl)ethanone", {"ACH": 2, "AC": 2, "ACCH3": 2, "CH3": 3, "C": 1, "CH3CO": 1}),
    ("acetophenone", {"ACH": 5, "AC": 1, "CH3CO": 1}),
    ("acetone", {"CH3CO": 1, "CH3": 1}),
    ("3-pentanone", {"CH3": 2, "CH2": 1, "CH2CO": 1}),
    ("2-butanone", {"CH3": 1, "CH2": 1, "CH3CO": 1})
]

@pytest.mark.CH2CO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_ch2co_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
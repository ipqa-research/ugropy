import ugropy as ug

import pytest


# =============================================================================
# 1- CH2 Main group: CH3, CH2, CH, C
# =============================================================================
ch2_trials = [
    # Gmehling L. 1982
    ("hexane", {"CH3": 2, "CH2": 4}),
    ("2-methylpropane", {"CH3": 3, "CH": 1}),
    ("2,2-dimethylpropane", {"CH3": 4, "C": 1}),
    # Horstmann et al. - 2005
    ("cyclohexane", {"CH2": 6}),
    # ugropy
    ("ethane", {"CH3": 2}),
    ("3-Ethyl-2,2-dimethylpentane", {"CH3": 5, "CH2": 2, "CH": 1, "C": 1}),
]

@pytest.mark.CH2
@pytest.mark.parametrize("name, result", ch2_trials)
def test_ch2(name, result):
    substance = ug.Substance(name)
    assert substance.unifac_groups == result
import ugropy as ug

import pytest


# =============================================================================
# 1- CH2 Main group: CH3, CH2, CH, C
# =============================================================================
# UNIFAC
trials_unifac = [
    ("ethane", {"CH3": 2}),
    ("3-Ethyl-2,2-dimethylpentane", {"CH3": 5, "CH2": 2, "CH": 1, "C": 1}),
    ("Decahydronaphthalene", {"CH2": 8, "CH": 2}),
    ("Dicyclohexylmethane", {"CH2": 11, "CH": 2}),
    ("hexane", {"CH3": 2, "CH2": 4}),
    ("2-methylpropane", {"CH3": 3, "CH": 1}),
    ("2,2-dimethylpropane", {"CH3": 4, "C": 1}),
    ("cyclohexane", {"CH2": 6}),
]

@pytest.mark.CH2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_ch2(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
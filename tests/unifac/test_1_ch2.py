import ugropy as ug

import pytest


# =============================================================================
# 1- CH2 Main group: CH3, CH2, CH, C
# =============================================================================
# UNIFAC
trials_unifac = [
    ("CCC(CC)C(C)(C)C", {"CH3": 5, "CH2": 2, "CH": 1, "C": 1}, "smiles"),
    ("C1CCC2CCCCC2C1", {"CH2": 8, "CH": 2}, "smiles"),
    ("C1CCC(CC1)CC2CCCCC2", {"CH2": 11, "CH": 2}, "smiles"),
    ("ethane", {"CH3": 2}, "name"),
    ("hexane", {"CH3": 2, "CH2": 4}, "name"),
    ("2-methylpropane", {"CH3": 3, "CH": 1}, "name"),
    ("2,2-dimethylpropane", {"CH3": 4, "C": 1}, "name"),
    ("cyclohexane", {"CH2": 6}, "name"),
]

@pytest.mark.CH2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
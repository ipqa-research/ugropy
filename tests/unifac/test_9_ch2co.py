import ugropy as ug

import pytest


# =============================================================================
# 9- CH2CO Main group: CH3CO, CH2CO
# =============================================================================

# UNIFAC
trials_unifac = [
    # Cyclopropanone
    ("C1CC1=O", {"CH2": 1, "CH2CO": 1}, "smiles"),
    # (9Z)-Cycloheptadec-9-en-1-one
    ("C1CCCC=CCCCCCCCC(=O)CCC1", {"CH2": 13, "CH=CH": 1, "CH2CO": 1}, "smiles"),
    # 1-(4-tert-butyl-2,6-dimethylphenyl)ethanone
    ("CC1=CC(=CC(=C1C(=O)C)C)C(C)(C)C", {"ACH": 2, "AC": 2, "ACCH3": 2, "CH3": 3, "C": 1, "CH3CO": 1}, "smiles"),
    # acetophenone
    ("CC(=O)C1=CC=CC=C1", {"ACH": 5, "AC": 1, "CH3CO": 1}, "smiles"),
    ("acetone", {"CH3CO": 1, "CH3": 1}, "name"),
    ("3-pentanone", {"CH3": 2, "CH2": 1, "CH2CO": 1}, "name"),
    ("2-butanone", {"CH3": 1, "CH2": 1, "CH3CO": 1}, "name")
]

@pytest.mark.CH2CO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2co_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
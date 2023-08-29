import ugropy as ug

import pytest


# =============================================================================
# 16- (C)3N Main group: CH3N, CH2N
# =============================================================================

# UNIFAC
trials_unifac = [
    ("CCN(CC(=O)CC)C1=CC=CC=C1", {"CH3": 2, "CH2N": 1, "AC": 1, "ACH": 5, "CH2CO": 1, "CH2": 1}, "smiles"),
    ("CCN(C(C)C)C(C)C", {"CH2N": 1, "CH3": 5, "CH": 2}, "smiles"),
    ("CCN(C)CC", {"CH3N": 1, "CH2": 2, "CH3": 2}, "smiles"),
    ("trimethylamine", {"CH3": 2, "CH3N": 1}, "name"),
    ("triethylamine", {"CH3": 3, "CH2": 2, "CH2N": 1}, "name"),
]

@pytest.mark.C3N
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
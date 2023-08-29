import ugropy as ug

import pytest


# =============================================================================
# 15- CNH Main group: CH3NH, CH2NH, CHNH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("CC(C)NCN", {"CHNH": 1, "CH3": 2, "CH2NH2": 1}, "smiles"),
    ("CCC(C)(C)NC(C)C", {"CH3": 5, "CHNH": 1, "CH2": 1, "C": 1}, "smiles"),
    ("CCNC(C)CC", {"CH3": 3, "CH2NH": 1, "CH": 1, "CH2": 1}, "smiles"),
    ("CCCNC", {"CH3NH": 1, "CH2": 2, "CH3": 1}, "smiles"),
    ("CN1CCCC1CC(=O)CC1CCCN1C", {"CH3N": 2, "CH2": 7, "CH": 2, "CH2CO": 1}, "smiles"),
    ("dimethylamine", {"CH3": 1, "CH3NH": 1}, "name"),
    ("diethylamine", {"CH3": 2, "CH2": 1, "CH2NH": 1}, "name"),
    ("diisopropylamine", {"CH3": 4, "CH": 1, "CHNH": 1}, "name"),
]

@pytest.mark.CNH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
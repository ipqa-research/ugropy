import pytest

import ugropy as ug


# =============================================================================
# 34- C=-C Main group: CH=-C, C=-C
# =============================================================================

# UNIFAC
trials_unifac = [
    ("CC#CC1=CC=CC=C1", {"ACH": 5, "AC": 1, "C=-C": 1, "CH3": 1}, "smiles"),
    # 1-hexyne
    ("CCCCC#C", {"CH3": 1, "CH2": 3, "CH=-C": 1}, "smiles"),
    # 2-hexyne
    ("CCCC#CC", {"CH3": 2, "CH2": 2, "C=-C": 1}, "smiles"),
]


@pytest.mark.alquine
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_alquine_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

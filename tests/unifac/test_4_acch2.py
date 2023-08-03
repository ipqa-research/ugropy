import ugropy as ug

import pytest


# =============================================================================
# 4- ACCH2 Main group: ACCH3, ACCH2, ACCH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1-Phenyl-2-methyl-1,3-butadiene", {"ACH": 5, "AC": 1, "CH=C": 1, "CH2=CH": 1, "CH3": 1}),
    ("9-(3-Butenyl)anthracene", {"ACH": 9, "ACCH2": 1, "AC": 4, "CH2": 1, "CH2=CH": 1}),
    ("9-Methylanthracene", {"ACH": 9, "ACCH3": 1, "AC": 4}),
    ("3-Methylbiphenyl", {"ACH": 9, "ACCH3": 1, "AC": 2}),
    ("1,2,4-Trimethyl-3-Ethylbenzene", {"ACH": 2, "ACCH3": 3, "ACCH2": 1, "CH3": 1}),
    ("1-t-Butyl-3-ethylbenzene", {"ACH": 4, "ACCH2": 1, "CH3": 4, "AC": 1, "C": 1}),
    ("1-Ethyl-2,3-dimethylbenzene", {"ACH": 3, "ACCH2": 1, "ACCH3": 2, "CH3": 1}),
    ("1-Ethyl-2-methylbenzene", {"ACH": 4, "ACCH3": 1, "ACCH2": 1, "CH3": 1}),
    ("Benzene, 1-ethyl-4-(1-methylethyl)-", {"ACH": 4, "ACCH": 1, "ACCH2": 1, "CH3": 3}),
    ("cumene", {"CH3": 2, "ACH": 5, "ACCH": 1}),
    ("ethylbenzene", {"CH3": 1, "ACH": 5, "ACCH2": 1}),
    ("toluene", {"ACH": 5, "ACCH3": 1}),
]

@pytest.mark.ACCH2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_cch2_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
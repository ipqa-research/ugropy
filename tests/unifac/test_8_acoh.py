import ugropy as ug

import pytest


# =============================================================================
# 8- ACOH Main group: ACOH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("Phenanthrene-3,4-diol", {"ACH": 8, "AC": 4, "ACOH": 2}),
    ("3-(tert-butyl)benzene-1,2-diol", {"ACH": 3, "AC": 1, "ACOH": 2, "CH3": 3, "C": 1}),
    ("[1,1'-Biphenyl]-2,3',4-triol", {"ACH": 7, "AC": 2, "ACOH": 3}),
    ("phenol", {"ACH": 5, "ACOH": 1}),
]

@pytest.mark.ACOH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_ch2_unifac(name, result):
    substance = ug.Substance(name)
    assert substance.unifac_groups == result
import ugropy as ug

import pytest


# =============================================================================
# 10- CHO Main group (aldehyde): CHO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("glyoxal", {"CHO": 2}),
    ("salicylaldehyde", {"ACH": 4, "ACOH": 1, "AC": 1, "CHO": 1}),
    ("2-Methyl-3-butenal", {"CH3": 1, "CH": 1, "CH2=CH": 1, "CHO": 1}),
    ("2-Hexyl-3-Phenyl-2-Propenal", {"ACH": 5, "AC": 1, "CH=C": 1, "CH2": 5, "CH3": 1, "CHO": 1}),
    ("cinnamaldehyde", {"ACH": 5, "AC": 1, "CH=CH": 1, "CHO": 1}),
    ("benzaldehyde", {"ACH": 5, "AC": 1, "CHO": 1}),
    ("cyclohexanecarbaldehyde", {"CH2": 5, "CH": 1, "CHO": 1,}),
    ("pentanal", {"CH3": 1, "CH2": 3, "CHO": 1}),
    ("3-methylbutanal", {"CH3": 2, "CH2": 1, "CH": 1, "CHO": 1}),
    ("acetaldehyde", {"CH3": 1, "CHO": 1}),
]

@pytest.mark.CHO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_cho_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
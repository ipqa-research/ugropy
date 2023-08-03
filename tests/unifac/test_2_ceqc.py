import ugropy as ug

import pytest


# =============================================================================
# 2- C=C Main group: CH2=CH, CH=CH, CH2=C, CH=C, C=C
# =============================================================================

# UNIFAC
trials_unifac = [
    ("alpha-pinene", {"CH3": 3, "CH2": 2, "CH": 2, "C": 1, "CH=C": 1}),
    ("d-limonene", {"CH3": 2, "CH2": 3, "CH": 1, "CH2=C": 1, "CH=C": 1}),
    ("2,3-Dimethyl-1,3-cyclohexadiene", {"CH2": 2, "CH=C": 2, "CH3": 2}),
    ("3,3'-(Pentane-1,3-diyl)dicyclohexene", {"CH3": 1, "CH2": 9, "CH":3, "CH=CH": 2}),
    ("cyclohexene", {"CH2": 4, "CH=CH": 1}),
    ("2,3-dimethylbutene-2", {"C=C": 1, "CH3": 4}),
    ("2-methyl-2-butene", {"CH=C": 1, "CH3": 3}),
    ("2-methyl-1-butene", {"CH2": 1, "CH2=C": 1, "CH3": 2}),
    ("2-hexene", {"CH2": 2, "CH=CH": 1, "CH3": 2}),
    ("1-hexene", {"CH2": 3, "CH2=CH": 1, "CH3": 1}),
]

@pytest.mark.CeqC
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_ceqc_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
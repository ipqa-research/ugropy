import ugropy as ug

import pytest


# =============================================================================
# 2- C=C Main group: CH2=CH, CH=CH, CH2=C, CH=C, C=C
# =============================================================================

# UNIFAC
trials_unifac = [
    # alpha-pinene
    ("CC1=CCC2CC1C2(C)C", {"CH3": 3, "CH2": 2, "CH": 2, "C": 1, "CH=C": 1}, "smiles"),
    # d-limonene
    ("CC1=CCC(CC1)C(=C)C", {"CH3": 2, "CH2": 3, "CH": 1, "CH2=C": 1, "CH=C": 1}, "smiles"),
    # 2,3-Dimethyl-1,3-cyclohexadiene
    ("CC1=CCCC=C1C", {"CH2": 2, "CH=C": 2, "CH3": 2}, "smiles"),
    # 3,3'-(Pentane-1,3-diyl)dicyclohexene
    ("CCC(CCC1CCCC=C1)C2CCCC=C2", {"CH3": 1, "CH2": 9, "CH":3, "CH=CH": 2}, "smiles"),
    # cyclohexene
    ("C1CCC=CC1", {"CH2": 4, "CH=CH": 1}, "smiles"),
    ("2,3-dimethylbutene-2", {"C=C": 1, "CH3": 4}, "name"),
    ("2-methyl-2-butene", {"CH=C": 1, "CH3": 3}, "name"),
    ("2-methyl-1-butene", {"CH2": 1, "CH2=C": 1, "CH3": 2}, "name"),
    ("2-hexene", {"CH2": 2, "CH=CH": 1, "CH3": 2}, "name"),
    ("1-hexene", {"CH2": 3, "CH2=CH": 1, "CH3": 1}, "name"),
]

@pytest.mark.CeqC
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ceqc_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
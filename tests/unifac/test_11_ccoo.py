import ugropy as ug

import pytest


# =============================================================================
# 11- CCOO Main group: CH3COO, CH2COO
# =============================================================================

def condicion_para_saltar(saltar):
    return saltar

# UNIFAC
trials_unifac = [
    # Tert-butyl acetate
    ("CC(=O)OC(C)(C)C", {"CH3COO": 1, "CH3": 3, "C": 1}, "smiles"),
    # triacetin
    ("CC(=O)OCC(COC(=O)C)OC(=O)C", {"CH3COO": 3, "CH2": 2, "CH": 1}, "smiles"),
    ("butyl propanoate", {"CH3": 2, "CH2": 3, "CH2COO": 1}, "name"),
    ("butyl acetate", {"CH3": 1, "CH2": 3, "CH3COO": 1}, "name")
]

@pytest.mark.CCOO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccoo_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
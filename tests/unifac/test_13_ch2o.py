import ugropy as ug

import pytest


# =============================================================================
# 13- CH2O Main group: CH2O, CH-O, THF
# =============================================================================

# UNIFAC
trials_unifac = [
    # Problematic ones
    #("Ethyl methyl carbonate", {"CH3": 1, "CH2": 1, "COO": 1, "CH3O": 1})
    ("Dimethyl carbonate", {"CH3": 1, "COO": 1, "CH3O": 1}),
]

@pytest.mark.CH2O
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_cho_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
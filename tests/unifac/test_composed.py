import ugropy as ug

import pytest


# =============================================================================
# Composed structures
# =============================================================================

# UNIFAC
trials_unifac = [
    # impossibles
    #("hydrogen peroxide", {}, "name"),
    #("methane", {}, "name"),
    #("C1(=CC=CC=C1)OC(C)(C)C", {}, "smiles"),
    # TODO: debug the alien
    # ("C13=C(C(=C(C(=C1C)C)CC2=C(C(=C(C(=C2C)CC)O[H])N([H])[H])C)C)CCCC3", {"CH3": 1, "CH2": 2, "AC": 1, "ACCH3": 5, "ACCH2": 4, "ACOH": 1, "ACNH2": 1}, "smiles"),
    # toluene - O - tertbutyl (eter between toluene and tertbutyl)
    ("C1(=CC=CC=C1)COC(C)(C)C", {"ACH": 5, "AC": 1, "CH2O": 1, "C": 1, "CH3": 3}, "smiles"),
    #diphenyl methane
    ("C1=CC=C(C=C1)CC2=CC=CC=C2", {"ACH": 10, "ACCH2": 1, "AC": 1}, "smiles"),
]

@pytest.mark.composed
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_problematics_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
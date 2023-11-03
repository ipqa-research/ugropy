import pytest

import ugropy as ug


# =============================================================================
# 40- CF2 Main group: CF3, CF2, CF
# =============================================================================

# UNIFAC
trials_unifac = [
    ("OC(F)(Br)I", {"CF": 1, "BR": 1, "I": 1, "OH": 1}, "smiles"),
    ("OC(O)(F)F", {"CF2": 1, "OH": 2}, "smiles"),
    ("OC(F)(F)F", {"CF3": 1, "OH": 1}, "smiles"),
    # Perfluorohexane
    (
        "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F",
        {"CF3": 2, "CF2": 4},
        "smiles",
    ),
    # Perfluoromethylcyclohexane
    (
        "C1(C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)(F)F)(C(F)(F)F)F",
        {"CF3": 1, "CF2": 5, "CF": 1},
        "smiles",
    ),
    # Impossibles
    ("FC(F)F", {}, "smiles"),
    ("FCF", {}, "smiles"),
    ("CF", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cf2_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

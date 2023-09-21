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
    (" Perfluorohexane", {"CF3": 2, "CF2": 4}, "name"),
    ("Perfluoromethylcyclohexane", {"CF3": 1, "CF2": 5, "CF": 1}, "name"),
    # Impossibles
    ("FC(F)F", {}, "smiles"),
    ("FCF", {}, "smiles"),
    ("CF", {}, "smiles"),
]


@pytest.mark.CF2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cf2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

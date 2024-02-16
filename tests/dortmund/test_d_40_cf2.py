import pytest

import ugropy as ug


# =============================================================================
# 40|CF2|[74]CF3 [75]CF2 [76]CF
# =============================================================================

# Dortmund
trials = [
    ("OC(F)(Br)I", {"CF": 1, "BR": 1, "I": 1, "OH (P)": 1}, "smiles"),
    ("OC(O)(F)F", {"CF2": 1, "OH (P)": 2}, "smiles"),
    ("OC(F)(F)F", {"CF3": 1, "OH (P)": 1}, "smiles"),
    # Perfluorohexane
    (
        "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F",
        {"CF3": 2, "CF2": 4},
        "smiles",
    ),
    # Impossibles
    ("FC(F)F", {}, "smiles"),
    ("FCF", {}, "smiles"),
    ("CF", {}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cf2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

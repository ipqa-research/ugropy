import pytest

import ugropy as ug


# =============================================================================
# 11|CCOO|[21]CH3COO [22]CH2COO
# =============================================================================

# Dortmund
trials = [
    # Aspirin
    (
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        {"CH3COO": 1, "AC": 2, "ACH": 4, "COOH": 1},
        "smiles",
    ),
    # Tert-butyl acetate
    ("CC(=O)OC(C)(C)C", {"CH3COO": 1, "CH3": 3, "C": 1}, "smiles"),
    # triacetin
    ("CC(=O)OCC(COC(=O)C)OC(=O)C", {"CH3COO": 3, "CH2": 2, "CH": 1}, "smiles"),
    # butyl propanoate
    ("CCCCOC(=O)CC", {"CH3": 2, "CH2": 3, "CH2COO": 1}, "smiles"),
    # butyl acetate
    ("CCCCOC(=O)C", {"CH3": 1, "CH2": 3, "CH3COO": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ccoo_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

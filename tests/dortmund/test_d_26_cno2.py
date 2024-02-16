import pytest

import ugropy as ug


# =============================================================================
# 26|CNO2|[54]CH3NO2 [55]CH2NO2 [56]CHNO2
# =============================================================================

# Dortmund
trials = [
    (
        "CC(O)C(O)[N+]([O-])=O",
        {"CH3": 1, "CH": 1, "OH (P)": 1, "OH (S)": 1, "CHNO2": 1},
        "smiles",
    ),
    (
        "CCC(O)[N+]([O-])=O",
        {"CH3": 1, "CH2": 1, "OH (P)": 1, "CHNO2": 1},
        "smiles",
    ),
    (
        "CCCC1=CC=C(C[N+]([O-])=O)C=C1",
        {"ACH": 4, "AC": 1, "CH2NO2": 1, "ACCH2": 1, "CH3": 1, "CH2": 1},
        "smiles",
    ),
    ("[O-][N+](=O)CC1=CC=CC=C1", {"ACH": 5, "AC": 1, "CH2NO2": 1}, "smiles"),
    # nitromethane
    ("C[N+](=O)[O-]", {"CH3NO2": 1}, "smiles"),
    # 1-nitropropane
    ("CCC[N+](=O)[O-]", {"CH3": 1, "CH2": 1, "CH2NO2": 1}, "smiles"),
    # 2-nitropropane
    ("CC(C)[N+](=O)[O-]", {"CH3": 2, "CHNO2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cno2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

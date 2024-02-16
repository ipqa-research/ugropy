import pytest

import ugropy as ug


# =============================================================================
# 20|COOH|[42]COOH
# =============================================================================

# Dortmund
trials = [
    ("OC(O)=O", {"COOH": 1, "OH (P)": 1}, "smiles"),
    ("CCOC(O)=O", {"COOH": 1, "CH2O": 1, "CH3": 1}, "smiles"),
    # 2,4-Diaminobutyric acid
    (
        "C(CN)C(C(=O)O)N",
        {"COOH": 1, "CHNH2": 1, "CH2": 1, "CH2NH2": 1},
        "smiles",
    ),
    # acetic acid
    ("CC(=O)O", {"CH3": 1, "COOH": 1}, "smiles"),
    # formic acid
    ("C(=O)O", {"HCOOH": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cooh_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 16|(C)3N|[34]CH3N [35]CH2N
# =============================================================================

# Dortmund
trials = [
    (
        "CCN(CC(=O)CC)C1=CC=CC=C1",
        {"CH3": 2, "CH2N": 1, "AC": 1, "ACH": 5, "CH2CO": 1, "CH2": 1},
        "smiles",
    ),
    ("CCN(C(C)C)C(C)C", {"CH2N": 1, "CH3": 5, "CH": 2}, "smiles"),
    ("CCN(C)CC", {"CH3N": 1, "CH2": 2, "CH3": 2}, "smiles"),
    # trimethylamine
    ("CN(C)C", {"CH3": 2, "CH3N": 1}, "smiles"),
    # triethylamine
    ("CCN(CC)CC", {"CH3": 3, "CH2": 2, "CH2N": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ch3nh2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 10|HCO|[20]HCO
# =============================================================================

# Dortmund
trials = [
    ("C(=O)C=O", {"HCO": 2}, "smiles"),
    # salicylaldehyde
    (
        "C1=CC=C(C(=C1)C=O)O",
        {"ACH": 4, "ACOH": 1, "AC": 1, "HCO": 1},
        "smiles",
    ),
    # 2-Methyl-3-butenal
    ("CC(C=C)C=O", {"CH3": 1, "CH": 1, "CH2=CH": 1, "HCO": 1}, "smiles"),
    # Cinnamaldehyde
    (
        "C1=CC=C(C=C1)C=CC=O",
        {"ACH": 5, "AC": 1, "CH=CH": 1, "HCO": 1},
        "smiles",
    ),
    # benzaldehyde
    ("C1=CC=C(C=C1)C=O", {"ACH": 5, "AC": 1, "HCO": 1}, "smiles"),
    # pentanal
    ("CCCCC=O", {"CH3": 1, "CH2": 3, "HCO": 1}, "smiles"),
    # 3-methylbutanal
    ("CC(C)CC=O", {"CH3": 2, "CH2": 1, "CH": 1, "HCO": 1}, "smiles"),
    # acetaldehyde
    ("CC=O", {"CH3": 1, "HCO": 1}, "smiles"),
    # 2-Hexyl-3-Phenyl-2-Propenal
    (
        r"CCCCCC\C(C=O)=C/C1=CC=CC=C1",
        {"ACH": 5, "AC": 1, "CH=C": 1, "CH2": 5, "CH3": 1, "HCO": 1},
        "smiles",
    ),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_hco_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 12|HCOO|[23]HCOO
# =============================================================================

# Dortmund
trials = [
    # phenyl formate
    ("C1=CC=C(C=C1)OC=O", {"ACH": 5, "AC": 1, "HCOO": 1}, "smiles"),
    # methyl formate
    ("COC=O", {"HCOO": 1, "CH3": 1}, "smiles"),
    # ethyl formate
    ("CCOC=O", {"HCOO": 1, "CH2": 1, "CH3": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cho_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 8- ACOH Main group: ACOH
# =============================================================================

# Dortmund
trials = [
    # Phenanthrene-3,4-diol
    (
        "C1=CC=C2C(=C1)C=CC3=C2C(=C(C=C3)O)O",
        {"ACH": 8, "AC": 4, "ACOH": 2},
        "smiles",
    ),
    # 3-(tert-butyl)benzene-1,2-diol
    (
        "CC(C)(C)C1=C(C(=CC=C1)O)O",
        {"ACH": 3, "AC": 1, "ACOH": 2, "CH3": 3, "C": 1},
        "smiles",
    ),
    # [1,1'-Biphenyl]-2,3',4-triol
    (
        "C1=CC(=CC(=C1)O)C2=C(C=C(C=C2)O)O",
        {"ACH": 7, "AC": 2, "ACOH": 3},
        "smiles",
    ),
    ("C1=CC=C(C=C1)O", {"ACH": 5, "ACOH": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_acoh_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
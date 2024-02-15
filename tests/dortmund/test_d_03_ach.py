import pytest

import ugropy as ug


# =============================================================================
# 3|ACH|[9]ACH [10]AC
# =============================================================================

# Dortmund
trials = [
    ("C1=CC2=CC=CC=CC2=C1", {"ACH": 8, "AC": 2}, "smiles"),
    ("C1=CC=CC=CC=CC=CC=CC=CC=CC=C1", {"ACH": 18}, "smiles"),
    ("C1=CC=CC=CC=CC=CC=CC=C1", {"ACH": 14}, "smiles"),
    # phenanthrene
    ("C1=CC=C2C(=C1)C=CC3=CC=CC=C32", {"ACH": 10, "AC": 4}, "smiles"),
    # Anthracene
    ("C1=CC=C2C=C3C=CC=CC3=CC2=C1", {"ACH": 10, "AC": 4}, "smiles"),
    # 1,1-Diphenylethylene
    (
        "C=C(C1=CC=CC=C1)C2=CC=CC=C2",
        {"ACH": 10, "AC": 2, "CH2=C": 1},
        "smiles",
    ),
    # biphenyl
    ("C1=CC=C(C=C1)C2=CC=CC=C2", {"ACH": 10, "AC": 2}, "smiles"),
    ("C1=CC=C2C=CC=CC2=C1", {"ACH": 8, "AC": 2}, "smiles"),  # Naphthalene
    ("C1=CC=CC=C1", {"ACH": 6}, "smiles"),  # benzene
    ("C=CC1=CC=CC=C1", {"AC": 1, "CH2=CH": 1, "ACH": 5}, "smiles"),  # styrene
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ach_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 41|COO|[77]COO
# =============================================================================

# Dortmund
trials = [
    # Procaine
    (
        "CCN(CC)CCOC(=O)C1=CC=C(N)C=C1",
        {
            "ACNH2": 1,
            "ACH": 4,
            "AC": 1,
            "COO": 1,
            "CH2": 3,
            "CH3": 2,
            "CH2N": 1,
        },
        "smiles",
    ),
    # Methyl acrylate
    ("COC(=O)C=C", {"CH3": 1, "CH2=CH": 1, "COO": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_coo_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

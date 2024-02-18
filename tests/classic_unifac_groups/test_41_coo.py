import pytest

import ugropy as ug


# =============================================================================
# 41- COO Main group: COO
# =============================================================================

# UNIFAC
trials_unifac = [
    # Ascorbic acid
    (
        "OCC(O)C1OC(=O)C(O)=C1O",
        {"COO": 1, "C=C": 1, "OH": 4, "CH": 2, "CH2": 1},
        "smiles",
    ),
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
    # Cocaine
    (
        "COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C",
        {
            "CH3": 1,
            "CH2": 3,
            "CH": 4,
            "CH3N": 1,
            "AC": 1,
            "ACH": 5,
            "COO": 2,
            "CH3N": 1,
        },
        "smiles",
    ),
    # Methyl acrylate
    ("COC(=O)C=C", {"CH3": 1, "CH2=CH": 1, "COO": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_coo_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result

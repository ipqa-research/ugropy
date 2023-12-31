import pytest

import ugropy as ug


# =============================================================================
# Composed structures
# =============================================================================

# UNIFAC
trials_unifac = [
    (
        "CC(C)CC1=CC=C(C=C1)C(C)OC(C)(C)C",
        {"CH3": 6, "CH": 1, "C": 1, "ACH": 4, "ACCH2": 1, "AC": 1, "CHO": 1},
        "smiles",
    ),
    # two solutions
    (
        "CCCC1=CC=C(CC(=O)OC)C=C1",
        [
            {"CH3": 2, "CH2": 1, "ACH": 4, "ACCH2": 1, "AC": 1, "CH2COO": 1},
            {"CH3": 2, "CH2": 1, "ACH": 4, "ACCH2": 2, "COO": 1},
        ],
        "smiles",
    ),
    (
        "C1=CC(=CC=C1COC(C)(C)C)CCC",
        {"ACH": 4, "ACCH2": 1, "AC": 1, "CH2O": 1, "CH3": 4, "CH2": 1, "C": 1},
        "smiles",
    ),
    (
        "C13=C(C=C(C=C1)CC2=CC=CC(=C2)CC)CCCC3",
        {"CH3": 1, "CH2": 2, "ACH": 7, "ACCH2": 4, "AC": 1},
        "smiles",
    ),
    (
        "C13=C(C(=C(C(=C1C)C)CC2=C(C(=C(C(=C2C)CC)O[H])N([H])[H])C)C)CCCC3",
        {
            "CH3": 1,
            "CH2": 2,
            "AC": 1,
            "ACCH3": 5,
            "ACCH2": 4,
            "ACOH": 1,
            "ACNH2": 1,
        },
        "smiles",
    ),
    # toluene - O - tertbutyl (eter between toluene and tertbutyl)
    (
        "C1(=CC=CC=C1)COC(C)(C)C",
        {"ACH": 5, "AC": 1, "CH2O": 1, "C": 1, "CH3": 3},
        "smiles",
    ),
    # diphenyl methane
    ("C1=CC=C(C=C1)CC2=CC=CC=C2", {"ACH": 10, "ACCH2": 1, "AC": 1}, "smiles"),
    (
        "C1(=CC=CC=C1)COC(C)(C)C",
        {"ACH": 5, "AC": 1, "CH2O": 1, "CH3": 3, "C": 1},
        "smiles",
    ),
    (
        "C1(=CC=CC=C1)C(OC(C)(C)C)C",
        {"ACH": 5, "AC": 1, "CHO": 1, "CH3": 4, "C": 1},
        "smiles",
    ),
    ("C12=CC=CC=C1COC2", {"ACH": 4, "AC": 1, "CH2O": 1, "ACCH2": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_composed_unifac(identifier, result, identifier_type):
    unifac_result = ug.get_unifac_groups(identifier, identifier_type)
    psrk_result = ug.get_psrk_groups(identifier, identifier_type)
    try:
        assert unifac_result == result
        assert psrk_result == result
    except ValueError:
        for uni, psrk, sol in zip(unifac_result, psrk_result, result):
            assert uni == sol
            assert psrk == sol

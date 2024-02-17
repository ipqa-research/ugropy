import pytest

import ugropy as ug


# =============================================================================
# 42|CY-CH2|[78]CY-CH2 [79]CY-CH [80]CY-C
# =============================================================================

# Dortmund
trials = [
    ("C1CC2CCCC3CCCC1C23", {"CY-CH2": 8, "CY-CH": 4}, "smiles"),
    ("C1C2CCCCC2C2CCCCC12", {"CY-CH2": 9, "CY-CH": 4}, "smiles"),
    ("C1C2CC1CCCC2", {"CY-CH2": 6, "CY-CH": 2}, "smiles"),
    ("C1CCCCCCCC1", {"CY-CH2": 9}, "smiles"),
    ("C1C2CC3CC1CC(C2)C3", {"CY-CH2": 6, "CY-CH": 4}, "smiles"),
    ("C12C3C1C1C2C31", {"CY-CH": 6}, "smiles"),
    ("C1CC2CC1CCC2", {"CY-CH2": 6, "CY-CH": 2}, "smiles"),
    ("C1CC2CC3CCC2CC13", {"CY-CH2": 6, "CY-CH": 4}, "smiles"),
    ("C12C3C4C1C1C2C3C41", {"CY-CH": 8}, "smiles"),
    ("C1CC1", {"CY-CH2": 3}, "smiles"),
    ("C1CCC1", {"CY-CH2": 4}, "smiles"),
    ("C1C2CC1CCCC2", {"CY-CH2": 6, "CY-CH": 2}, "smiles"),
    (
        "CC12C3CCC4CCC1C234",
        {"CH3": 1, "CY-CH2": 4, "CY-CH": 3, "CY-C": 2},
        "smiles",
    ),
    ("C1CCC2CCCCC2C1", {"CY-CH2": 8, "CY-CH": 2}, "smiles"),
    ("C1CCC(CC1)CC2CCCCC2", {"CY-CH2": 10, "CH2": 1, "CY-CH": 2}, "smiles"),
    ("C1CCCCC1", {"CY-CH2": 6}, "smiles"),  # cyclohexane
    # alpha-pinene
    (
        "CC1=CCC2CC1C2(C)C",
        {"CH3": 3, "CY-CH2": 2, "CY-CH": 2, "CY-C": 1, "CH=C": 1},
        "smiles",
    ),
    # d-limonene
    (
        "CC1=CCC(CC1)C(=C)C",
        {"CH3": 2, "CY-CH2": 3, "CY-CH": 1, "CH2=C": 1, "CH=C": 1},
        "smiles",
    ),
    # 2,3-Dimethyl-1,3-cyclohexadiene
    ("CC1=CCCC=C1C", {"CY-CH2": 2, "CH=C": 2, "CH3": 2}, "smiles"),
    # 3,3'-(Pentane-1,3-diyl)dicyclohexene
    (
        "CCC(CCC1CCCC=C1)C2CCCC=C2",
        {"CH3": 1, "CH2": 3, "CY-CH2": 6, "CH": 1, "CY-CH": 2, "CH=CH": 2},
        "smiles",
    ),
    ("C1CCC=CC1", {"CY-CH2": 4, "CH=CH": 1}, "smiles"),  # cyclohexene
    # 1,2-Cyclohexanediol, 4-tert-butyl-1-phenyl-, stereoisomer
    (
        "CC(C)(C)C1CCC(C(C1)O)(C2=CC=CC=C2)O",
        {
            "CH3": 3,
            "CY-CH2": 3,
            "CY-CH": 2,
            "C": 1,
            "CY-C": 1,
            "OH (S)": 1,
            "OH (T)": 1,
            "AC": 1,
            "ACH": 5,
        },
        "smiles",
    ),
    # cyclohexanecarbaldehyde
    (
        "C1CCC(CC1)C=O",
        {
            "CY-CH2": 5,
            "CY-CH": 1,
            "HCO": 1,
        },
        "smiles",
    ),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cych2_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

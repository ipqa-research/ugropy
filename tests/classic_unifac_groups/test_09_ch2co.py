import pytest

import ugropy as ug


# =============================================================================
# 9- CH2CO Main group: CH3CO, CH2CO
# =============================================================================

# UNIFAC
trials_unifac = [
    (
        "CC(C)CC(=O)OCC(=O)C(C)O",
        {"CH3": 3, "CH": 2, "CH2COO": 1, "CH2CO": 1, "OH": 1},
        "smiles",
    ),
    (
        "O=C(CC1=CC=CC=C1)CC1=CC=CC=C1",
        {"ACH": 10, "AC": 1, "ACCH2": 1, "CH2CO": 1},
        "smiles",
    ),
    ("CC(=O)CC1=CC=CC=C1", {"CH3CO": 1, "ACH": 5, "ACCH2": 1}, "smiles"),
    (
        "CC(C)(C)C(=O)CC1=CC=CC=C1",
        {"CH3": 3, "ACH": 5, "AC": 1, "CH2CO": 1, "C": 1},
        "smiles",
    ),
    # Cyclopropanone
    ("C1CC1=O", {"CH2": 1, "CH2CO": 1}, "smiles"),
    # (9Z)-Cycloheptadec-9-en-1-one
    (
        "C1CCCC=CCCCCCCCC(=O)CCC1",
        {"CH2": 13, "CH=CH": 1, "CH2CO": 1},
        "smiles",
    ),
    # 1-(4-tert-butyl-2,6-dimethylphenyl)ethanone
    (
        "CC1=CC(=CC(=C1C(=O)C)C)C(C)(C)C",
        {"ACH": 2, "AC": 2, "ACCH3": 2, "CH3": 3, "C": 1, "CH3CO": 1},
        "smiles",
    ),
    # acetophenone
    ("CC(=O)C1=CC=CC=C1", {"ACH": 5, "AC": 1, "CH3CO": 1}, "smiles"),
    # acetone
    ("CC(=O)C", {"CH3CO": 1, "CH3": 1}, "smiles"),
    # 3-pentanones
    ("CCC(=O)CC", {"CH3": 2, "CH2": 1, "CH2CO": 1}, "smiles"),
    # 2-butanone
    ("CCC(=O)C", {"CH3": 1, "CH2": 1, "CH3CO": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2co_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

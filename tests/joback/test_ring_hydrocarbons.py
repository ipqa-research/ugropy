import pytest

import ugropy as ug


# =============================================================================
# ring=CH-, ring>CH-, ring=C<, ring>C<, ring-CH2-
# =============================================================================
# Joback
trials = [
    # Atrolactic acid
    ("C1=CC2=CC=CC=CC2=C1", {"ring=CH-": 8, "ring=C<": 2}, "smiles"),
    ("C1CC2CCCC3CCCC1C23", {"ring-CH2-": 8, "ring>CH-": 4}, "smiles"),
    ("C1C2CCCCC2C2CCCCC12", {"ring-CH2-": 9, "ring>CH-": 4}, "smiles"),
    ("C1C2CC1CCCC2", {"ring-CH2-": 6, "ring>CH-": 2}, "smiles"),
    ("C1CCCCCCCC1", {"ring-CH2-": 9}, "smiles"),
    ("C1CCCCCCCC1", {"ring-CH2-": 9}, "smiles"),
    ("C1C2CC3CC1CC(C2)C3", {"ring-CH2-": 6, "ring>CH-": 4}, "smiles"),
    ("C12C3C1C1C2C31", {"ring>CH-": 6}, "smiles"),
    ("C1CC2CC1CCC2", {"ring-CH2-": 6, "ring>CH-": 2}, "smiles"),
    ("C1CC2CC3CCC2CC13", {"ring-CH2-": 6, "ring>CH-": 4}, "smiles"),
    ("C12C3C4C1C1C2C3C41", {"ring>CH-": 8}, "smiles"),
    ("C1CC1", {"ring-CH2-": 3}, "smiles"),
    ("C1CCC1", {"ring-CH2-": 4}, "smiles"),
    ("C1C2CC1CCCC2", {"ring-CH2-": 6, "ring>CH-": 2}, "smiles"),
    ("CC1(C)CCCCC1", {"ring-CH2-": 5, "ring>C<": 1, "-CH3": 2}, "smiles"),
    (
        "CC12C3CCC4CCC1C234",
        {"-CH3": 1, "ring-CH2-": 4, "ring>CH-": 3, "ring>C<": 2},
        "smiles",
    ),
    ("C1CCC2CCCCC2C1", {"ring-CH2-": 8, "ring>CH-": 2}, "smiles"),
    (
        "C1CCC(CC1)CC2CCCCC2",
        {"-CH2-": 1, "ring-CH2-": 10, "ring>CH-": 2},
        "smiles",
    ),
    ("C1CCCCC1", {"ring-CH2-": 6}, "smiles"),  # cyclohexane
    ("C1=CC2=CC=CC=CC2=C1", {"ring=CH-": 8, "ring=C<": 2}, "smiles"),
    ("C1=CC=CC=CC=CC=CC=CC=CC=CC=C1", {"ring=CH-": 18}, "smiles"),
    ("C1=CC=CC=CC=CC=CC=CC=C1", {"ring=CH-": 14}, "smiles"),
    # phenanthrene
    (
        "C1=CC=C2C(=C1)C=CC3=CC=CC=C32",
        {"ring=CH-": 10, "ring=C<": 4},
        "smiles",
    ),
    # Anthracene
    ("C1=CC=C2C=C3C=CC=CC3=CC2=C1", {"ring=CH-": 10, "ring=C<": 4}, "smiles"),
    # 1,1-Diphenylethylene
    (
        "C=C(C1=CC=CC=C1)C2=CC=CC=C2",
        {"ring=CH-": 10, "ring=C<": 2, "=CH2": 1, "=C<": 1},
        "smiles",
    ),
    # biphenyl
    ("C1=CC=C(C=C1)C2=CC=CC=C2", {"ring=CH-": 10, "ring=C<": 2}, "smiles"),
    (
        "C1=CC=C2C=CC=CC2=C1",
        {"ring=CH-": 8, "ring=C<": 2},
        "smiles",
    ),  # Naphthalene
    ("C1=CC=CC=C1", {"ring=CH-": 6}, "smiles"),  # benzene
    (
        "C=CC1=CC=CC=C1",
        {"ring=C<": 1, "=CH2": 1, "=CH-": 1, "ring=CH-": 5},
        "smiles",
    ),  # styrene
    # Atrolactic acid
    (
        "CC(C1=CC=CC=C1)(C(=O)O)O",
        {
            "ring=CH-": 5,
            "ring=C<": 1,
            "-OH (alcohol)": 1,
            "-CH3": 1,
            ">C<": 1,
            "-COOH (acid)": 1,
        },
        "smiles",
    ),
    (
        "C=CC=CC1=CC=CC=C1",
        {"ring=CH-": 5, "ring=C<": 1, "=CH-": 3, "=CH2": 1},
        "smiles",
    ),
    # 9-(3-Butenyl)anthracene
    (
        "C=CCCC1=C2C=CC=CC2=CC3=CC=CC=C31",
        {"ring=CH-": 9, "ring=C<": 5, "-CH2-": 2, "=CH2": 1, "=CH-": 1},
        "smiles",
    ),
    # 9-Methylanthracene
    (
        "CC1=C2C=CC=CC2=CC3=CC=CC=C13",
        {"ring=CH-": 9, "-CH3": 1, "ring=C<": 5},
        "smiles",
    ),
    # 3-Methylbiphenyl
    (
        "CC1=CC(=CC=C1)C2=CC=CC=C2",
        {"ring=CH-": 9, "-CH3": 1, "ring=C<": 3},
        "smiles",
    ),
    # 1,2,4-Trimethyl-3-Ethylbenzene
    (
        "CCC1=C(C=CC(=C1C)C)C",
        {"ring=CH-": 2, "ring=C<": 4, "-CH3": 4, "-CH2-": 1},
        "smiles",
    ),
    # 1-t-Butyl-3-ethylbenzene
    (
        "CCC1=CC(=CC=C1)C(C)(C)C",
        {"ring=CH-": 4, "ring=C<": 2, "-CH2-": 1, "-CH3": 4, ">C<": 1},
        "smiles",
    ),
    # 1-Ethyl-2,3-dimethylbenzene
    (
        "CCC1=CC=CC(=C1C)C",
        {"ring=CH-": 3, "ring=C<": 3, "-CH2-": 1, "-CH3": 3},
        "smiles",
    ),
    # 1-Ethyl-2-methylbenzene
    (
        "CCC1=CC=CC=C1C",
        {"ring=CH-": 4, "ring=C<": 2, "-CH3": 2, "-CH2-": 1},
        "smiles",
    ),
    # Benzene, 1-ethyl-4-(1-methylethyl)-
    (
        "CCC1=CC=C(C=C1)C(C)C",
        {"ring=CH-": 4, "ring=C<": 2, ">CH-": 1, "-CH2-": 1, "-CH3": 3},
        "smiles",
    ),
    # cumene
    (
        "CC(C)C1=CC=CC=C1",
        {"-CH3": 2, "ring=C<": 1, "ring=CH-": 5, ">CH-": 1},
        "smiles",
    ),
    # ethylbenzene
    (
        "CCC1=CC=CC=C1",
        {"-CH3": 1, "ring=C<": 1, "ring=CH-": 5, "-CH2-": 1},
        "smiles",
    ),
    # toluene
    ("CC1=CC=CC=C1", {"ring=CH-": 5, "ring=C<": 1, "-CH3": 1}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_cyclic_hydrocarbon(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result

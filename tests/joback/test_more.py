import pytest

import ugropy as ug


# Joback
trials = [
    (
        "C1(=CC=CC=C1)OC(C)(C)C",
        {
            "-CH3": 3,
            ">C<": 1,
            "-O- (non-ring)": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
        },
        "smiles",
    ),
    (
        "CC(C)CC1=CC=C(C=C1)C(C)OC(C)(C)C",
        {
            "-CH3": 6,
            ">C<": 1,
            ">CH-": 2,
            "-CH2-": 1,
            "-O- (non-ring)": 1,
            "ring=CH-": 4,
            "ring=C<": 2,
        },
        "smiles",
    ),
    (
        "CCCC1=CC=C(CC(=O)OC)C=C1",
        {
            "-CH3": 2,
            "-CH2-": 3,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-COO- (ester)": 1,
        },
        "smiles",
    ),
    (
        "C1=CC(=CC=C1COC(C)(C)C)CCC",
        {
            "-CH3": 4,
            ">C<": 1,
            "-CH2-": 3,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-O- (non-ring)": 1,
        },
        "smiles",
    ),
    (
        "C13=C(C=C(C=C1)CC2=CC=CC(=C2)CC)CCCC3",
        {"-CH3": 1, "-CH2-": 2, "ring=CH-": 7, "ring=C<": 5, "ring-CH2-": 4},
        "smiles",
    ),
    (
        "C12=CC=CC=C1COC2",
        {"ring=CH-": 4, "ring=C<": 2, "-O- (ring)": 1, "ring-CH2-": 2},
        "smiles",
    ),
    # Imidazol
    ("N1C=CN=C1", {"-N= (ring)": 1, ">NH (ring)": 1, "ring=CH-": 3}, "smiles"),
    # thiophene
    ("C1=CSC=C1", {"-S- (ring)": 1, "ring=CH-": 4}, "smiles"),
    # 2-methylthiophene
    (
        "CC1=CC=CS1",
        {"-S- (ring)": 1, "ring=CH-": 3, "ring=C<": 1, "-CH3": 1},
        "smiles",
    ),
    # 2,3-dimethylthiophene
    (
        "CC1=C(SC=C1)C",
        {"-S- (ring)": 1, "ring=CH-": 2, "ring=C<": 2, "-CH3": 2},
        "smiles",
    ),
    (
        "OC1=CSC=C1",
        {"-S- (ring)": 1, "ring=CH-": 3, "ring=C<": 1, "-OH (phenol)": 1},
        "smiles",
    ),
    (
        "OC1=CC=CS1",
        {"-S- (ring)": 1, "ring=CH-": 3, "ring=C<": 1, "-OH (phenol)": 1},
        "smiles",
    ),
    (
        "OC1=CSC=C1O",
        {"-S- (ring)": 1, "ring=CH-": 2, "ring=C<": 2, "-OH (phenol)": 2},
        "smiles",
    ),
    (
        "OC1=CC(O)=CS1",
        {"-S- (ring)": 1, "ring=CH-": 2, "ring=C<": 2, "-OH (phenol)": 2},
        "smiles",
    ),
    (
        "OC1=CC=C(O)S1",
        {"-S- (ring)": 1, "ring=CH-": 2, "ring=C<": 2, "-OH (phenol)": 2},
        "smiles",
    ),
    # acrylonitrile
    ("C=CC#N", {"-CN": 1, "=CH-": 1, "=CH2": 1}, "smiles"),
    # 1,2-ethanediol
    ("C(CO)O", {"-CH2-": 2, "-OH (alcohol)": 2}, "smiles"),
    # furfural
    (
        "C1=COC(=C1)C=O",
        {"ring=CH-": 3, "ring=C<": 1, "-O- (ring)": 1, "O=CH- (aldehyde)": 1},
        "smiles",
    ),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_more(identifier, result, identifier_type):
    assert ug.get_groups(ug.joback, identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# -COO- (ester)
# =============================================================================
# Joback
trials = [
    # Ascorbic acid
    (
        "OCC(O)C1OC(=O)C(O)=C1O",
        {
            "-COO- (ester)": 1,
            "ring=C<": 2,
            "-OH (alcohol)": 4,
            ">CH-": 1,
            "ring>CH-": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    # Procaine
    (
        "CCN(CC)CCOC(=O)C1=CC=C(N)C=C1",
        {
            "-NH2": 1,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-COO- (ester)": 1,
            "-CH2-": 4,
            "-CH3": 2,
            ">N- (non-ring)": 1,
        },
        "smiles",
    ),
    # Cocaine
    (
        "COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C",
        {},
        "smiles",
    ),
    # Methyl acrylate
    (
        "COC(=O)C=C",
        {"-CH3": 1, "=CH2": 1, "=CH-": 1, "-COO- (ester)": 1},
        "smiles",
    ),
    # Aspirin
    (
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        {
            "-CH3": 1,
            "-COO- (ester)": 1,
            "ring=C<": 2,
            "ring=CH-": 4,
            "-COOH (acid)": 1,
        },
        "smiles",
    ),
    # Tert-butyl acetate
    ("CC(=O)OC(C)(C)C", {"-COO- (ester)": 1, "-CH3": 4, ">C<": 1}, "smiles"),
    # triacetin
    (
        "CC(=O)OCC(COC(=O)C)OC(=O)C",
        {"-CH3": 3, "-COO- (ester)": 3, "-CH2-": 2, ">CH-": 1},
        "smiles",
    ),
    # butyl propanoate
    ("CCCCOC(=O)CC", {"-CH3": 2, "-CH2-": 4, "-COO- (ester)": 1}, "smiles"),
    # butyl acetate
    ("CCCCOC(=O)C", {"-CH3": 2, "-CH2-": 3, "-COO- (ester)": 1}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_esters(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# O=CH- (aldehyde)
# =============================================================================
# Joback
trials = [
    ("C(=O)C=O", {"O=CH- (aldehyde)": 2}, "smiles"),
    # salicylaldehyde
    (
        "C1=CC=C(C(=C1)C=O)O",
        {
            "ring=CH-": 4,
            "-OH (phenol)": 1,
            "ring=C<": 2,
            "O=CH- (aldehyde)": 1,
        },
        "smiles",
    ),
    # 2-Methyl-3-butenal
    (
        "CC(C=C)C=O",
        {"-CH3": 1, ">CH-": 1, "=CH2": 1, "=CH-": 1, "O=CH- (aldehyde)": 1},
        "smiles",
    ),
    # Cinnamaldehyde
    (
        "C1=CC=C(C=C1)C=CC=O",
        {"ring=CH-": 5, "ring=C<": 1, "=CH-": 2, "O=CH- (aldehyde)": 1},
        "smiles",
    ),
    # benzaldehyde
    (
        "C1=CC=C(C=C1)C=O",
        {"ring=CH-": 5, "ring=C<": 1, "O=CH- (aldehyde)": 1},
        "smiles",
    ),
    # cyclohexanecarbaldehyde
    (
        "C1CCC(CC1)C=O",
        {
            "ring-CH2-": 5,
            "ring>CH-": 1,
            "O=CH- (aldehyde)": 1,
        },
        "smiles",
    ),
    # pentanal
    ("CCCCC=O", {"-CH3": 1, "-CH2-": 3, "O=CH- (aldehyde)": 1}, "smiles"),
    # 3-methylbutanal
    (
        "CC(C)CC=O",
        {"-CH3": 2, "-CH2-": 1, ">CH-": 1, "O=CH- (aldehyde)": 1},
        "smiles",
    ),
    # acetaldehyde
    ("CC=O", {"-CH3": 1, "O=CH- (aldehyde)": 1}, "smiles"),
    (
        r"CCCCCC\C(C=O)=C/C1=CC=CC=C1",
        {
            "-CH3": 1,
            "O=CH- (aldehyde)": 1,
            "-CH2-": 5,
            "=C<": 1,
            "=CH-": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
        },
        "smiles",
    ),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_aldehydes(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result

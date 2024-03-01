import pytest

from ugropy import get_groups, joback


# =============================================================================
# >C=O (non-ring), >C=O (ring)
# =============================================================================
# Joback
trials = [
    # (9Z)-Cycloheptadec-9-en-1-one
    (
        "C1CCCC=CCCCCCCCC(=O)CCC1",
        {"ring-CH2-": 14, "ring=CH-": 2, ">C=O (ring)": 1},
        "smiles",
    ),
    (
        "CC(C)CC(=O)OCC(=O)C(C)O",
        {
            "-CH3": 3,
            ">CH-": 2,
            "-CH2-": 2,
            "-COO- (ester)": 1,
            ">C=O (non-ring)": 1,
            "-OH (alcohol)": 1,
        },
        "smiles",
    ),
    (
        "O=C(CC1=CC=CC=C1)CC1=CC=CC=C1",
        {"ring=CH-": 10, "ring=C<": 2, "-CH2-": 2, ">C=O (non-ring)": 1},
        "smiles",
    ),
    (
        "CC(=O)CC1=CC=CC=C1",
        {
            ">C=O (non-ring)": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-CH2-": 1,
            "-CH3": 1,
        },
        "smiles",
    ),
    (
        "CC(C)(C)C(=O)CC1=CC=CC=C1",
        {
            "-CH3": 3,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-CH2-": 1,
            ">C=O (non-ring)": 1,
            ">C<": 1,
        },
        "smiles",
    ),
    # Cyclopropanone
    ("C1CC1=O", {"ring-CH2-": 2, ">C=O (ring)": 1}, "smiles"),
    # 1-(4-tert-butyl-2,6-dimethylphenyl)ethanone
    (
        "CC1=CC(=CC(=C1C(=O)C)C)C(C)(C)C",
        {
            "ring=CH-": 2,
            "ring=C<": 4,
            "-CH3": 6,
            ">C<": 1,
            ">C=O (non-ring)": 1,
        },
        "smiles",
    ),
    # acetophenone
    (
        "CC(=O)C1=CC=CC=C1",
        {"ring=CH-": 5, "ring=C<": 1, "-CH3": 1, ">C=O (non-ring)": 1},
        "smiles",
    ),
    # acetone
    ("CC(=O)C", {"-CH3": 2, ">C=O (non-ring)": 1}, "smiles"),
    # 3-pentanones
    ("CCC(=O)CC", {"-CH3": 2, "-CH2-": 2, ">C=O (non-ring)": 1}, "smiles"),
    # 2-butanone
    ("CCC(=O)C", {"-CH3": 2, "-CH2-": 1, ">C=O (non-ring)": 1}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_ketones(identifier, result, identifier_type):
    assert get_groups(joback, identifier, identifier_type).subgroups == result

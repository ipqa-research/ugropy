import pytest

import ugropy as ug


# =============================================================================
# -SH, -S- (non-ring), -S- (ring)
# =============================================================================
# Joback
trials = [
    # Dimethylsulfide
    ("CSC", {"-CH3": 2, "-S- (non-ring)": 1}, "smiles"),
    # Diethylsulfide
    ("CCSCC", {"-CH3": 2, "-CH2-": 2, "-S- (non-ring)": 1}, "smiles"),
    # Isopropyl Sulfide
    ("CC(C)SC(C)C", {"-CH3": 4, ">CH-": 2, "-S- (non-ring)": 1}, "smiles"),
    # carbon disulfide
    ("C(=S)=S", {}, "smiles"),
    (
        "CCCC1=CC=C(C=C1)C(C)S",
        {
            "-CH3": 2,
            ">CH-": 1,
            "-CH2-": 2,
            "-SH": 1,
            "ring=CH-": 4,
            "ring=C<": 2,
        },
        "smiles",
    ),
    (
        "CCCC1=CC=C(CS)C=C1",
        {"-CH3": 1, "-CH2-": 3, "-SH": 1, "ring=CH-": 4, "ring=C<": 2},
        "smiles",
    ),
    (
        "CC1=CC=C(CS)C=C1",
        {"-CH3": 1, "-CH2-": 1, "-SH": 1, "ring=CH-": 4, "ring=C<": 2},
        "smiles",
    ),
    (
        "SCC1=CC=NC=C1",
        {"-CH2-": 1, "-SH": 1, "ring=CH-": 4, "ring=C<": 1, "-N= (ring)": 1},
        "smiles",
    ),
    # methanethiol
    ("CS", {"-CH3": 1, "-SH": 1}, "smiles"),
    # ethanethiol
    ("CCS", {"-CH2-": 1, "-SH": 1, "-CH3": 1}, "smiles"),
    (
        "CCCC1=CC=C(C=C1)C(C)S",
        {
            "-CH3": 2,
            ">CH-": 1,
            "-CH2-": 2,
            "-SH": 1,
            "ring=CH-": 4,
            "ring=C<": 2,
        },
        "smiles",
    ),
    # sulfolane
    ("C1CCS(=O)(=O)C1", {}, "smiles"),
    # 2,4-dimethylsulfolane
    (
        "CC1CC(S(=O)(=O)C1)C",
        {},
        "smiles",
    ),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_sulfur(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result

import pytest

from ugropy import get_groups, joback


# =============================================================================
# -O- (non-ring), -O- (ring)
# =============================================================================
# Joback
trials = [
    # 4-flavanol
    (
        "OC1CC(OC2=CC=CC=C12)C1=CC=CC=C1",
        {
            "ring=CH-": 9,
            "ring=C<": 3,
            "ring>CH-": 2,
            "-O- (ring)": 1,
            "ring-CH2-": 1,
            "-OH (alcohol)": 1,
        },
        "smiles",
    ),
    (
        "O[C@@H]1CO[C@H](O)[C@@H](O)[C@@H]1O",
        {"ring-CH2-": 1, "-O- (ring)": 1, "ring>CH-": 4, "-OH (alcohol)": 4},
        "smiles",
    ),
    ("C1COCCOCCOCCOC1", {"ring-CH2-": 9, "-O- (ring)": 4}, "smiles"),
    ("C1COCCO1", {"-O- (ring)": 2, "ring-CH2-": 4}, "smiles"),
    ("CCOCOCC", {"-CH3": 2, "-O- (non-ring)": 2, "-CH2-": 3}, "smiles"),
    ("C1COCO1", {"-O- (ring)": 2, "ring-CH2-": 3}, "smiles"),
    ("C1COCCOCCOCCOCCOCCO1", {"-O- (ring)": 6, "ring-CH2-": 12}, "smiles"),
    # tetrahydrofuran
    ("C1CCOC1", {"-O- (ring)": 1, "ring-CH2-": 4}, "smiles"),
    # diisopropyl ether
    ("CC(C)OC(C)C", {"-CH3": 4, ">CH-": 2, "-O- (non-ring)": 1}, "smiles"),
    # diethyl ether
    ("CCOCC", {"-CH3": 2, "-CH2-": 2, "-O- (non-ring)": 1}, "smiles"),
    # dimethyl ether
    ("COC", {"-CH3": 2, "-O- (non-ring)": 1}, "smiles"),
    # 2H-Pyran, 2-(cyclohexyloxy)tetrahydro-
    (
        "C1CCC(CC1)OC2CCCCO2",
        {"ring-CH2-": 9, "ring>CH-": 2, "-O- (ring)": 1, "-O- (non-ring)": 1},
        "smiles",
    ),
    # Problematic ones
    (
        "COC(=O)OC1=CC=CC=C1",
        {
            "-CH3": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
            "ring=C<": 1,
            "ring=CH-": 5,
        },
        "smiles",
    ),
    (
        "O=C1OCCCO1",
        {"-COO- (ester)": 1, "-O- (ring)": 1, "ring-CH2-": 3},
        "smiles",
    ),
    (
        "CCOC(=O)OC1=CC=CC=C1",
        {
            "-CH3": 1,
            "-CH2-": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
            "ring=C<": 1,
            "ring=CH-": 5,
        },
        "smiles",
    ),
    (
        "CC(C)OC(=O)OC1=CC=CC=C1",
        {
            "-CH3": 2,
            ">CH-": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
            "ring=C<": 1,
            "ring=CH-": 5,
        },
        "smiles",
    ),
    (
        "CC(C)(C)OC(=O)OC1=CC=CC=C1",
        {
            "-CH3": 3,
            ">C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
            "ring=C<": 1,
            "ring=CH-": 5,
        },
        "smiles",
    ),
    # Benzyl 2-hydroxyethyl carbonate
    (
        "C1=CC=C(C=C1)COC(=O)OCCO",
        {
            "ring=C<": 1,
            "ring=CH-": 5,
            "-COO- (ester)": 1,
            "-O- (non-ring)": 1,
            "-CH2-": 3,
            "-OH (alcohol)": 1,
        },
        "smiles",
    ),
    # tert-Butyl ethyl carbonate
    (
        "CCOC(=O)OC(C)(C)C",
        {
            "-CH3": 4,
            ">C<": 1,
            "-COO- (ester)": 1,
            "-CH2-": 1,
            "-O- (non-ring)": 1,
        },
        "smiles",
    ),
    # Ethyl phenyl carbonate
    (
        "CCOC(=O)OC1=CC=CC=C1",
        {
            "-CH3": 1,
            "ring=C<": 1,
            "ring=CH-": 5,
            "-COO- (ester)": 1,
            "-CH2-": 1,
            "-O- (non-ring)": 1,
        },
        "smiles",
    ),
    # Carbonic acid, ethyl 2,3,6-trimethylcyclohexyl ester
    (
        "CCOC(=O)OC1C(CCC(C1C)C)C",
        {
            "-CH3": 4,
            "-CH2-": 1,
            "ring-CH2-": 2,
            "ring>CH-": 4,
            "-COO- (ester)": 1,
            "-O- (non-ring)": 1,
        },
        "smiles",
    ),
    # Diethyl carbonate
    (
        "CCOC(=O)OCC",
        {"-CH3": 2, "-CH2-": 2, "-COO- (ester)": 1, "-O- (non-ring)": 1},
        "smiles",
    ),
    # Methyl phenyl carbonate
    (
        "COC(=O)OC1=CC=CC=C1",
        {
            "ring=C<": 1,
            "ring=CH-": 5,
            "-COO- (ester)": 1,
            "-CH3": 1,
            "-O- (non-ring)": 1,
        },
        "smiles",
    ),
    # tert-Butyl methyl carbonate
    (
        "CC(C)(C)OC(=O)OC",
        {"-CH3": 4, ">C<": 1, "-COO- (ester)": 1, "-O- (non-ring)": 1},
        "smiles",
    ),
    # Methyl isopropyl carbonate
    (
        "CC(C)OC(=O)OC",
        {"-CH3": 3, ">CH-": 1, "-COO- (ester)": 1, "-O- (non-ring)": 1},
        "smiles",
    ),
    # Ethyl methyl carbonate
    (
        "CCOC(=O)OC",
        {"-CH2-": 1, "-COO- (ester)": 1, "-CH3": 2, "-O- (non-ring)": 1},
        "smiles",
    ),
    # Dimethyl carbonate
    (
        "COC(=O)OC",
        {"-COO- (ester)": 1, "-CH3": 2, "-O- (non-ring)": 1},
        "smiles",
    ),
    (
        "COCOC(C)OCOC",
        {"-CH3": 3, "-O- (non-ring)": 4, "-CH2-": 2, ">CH-": 1},
        "smiles",
    ),
    (
        "CC(C)OCOC(C)OCOC(C)C",
        {"-CH3": 5, ">CH-": 3, "-CH2-": 2, "-O- (non-ring)": 4},
        "smiles",
    ),
    (
        "CC(C)OCOCC(OCOC(C)C)OCOC(C)C",
        {"-CH3": 6, ">CH-": 4, "-CH2-": 4, "-O- (non-ring)": 6},
        "smiles",
    ),
    (
        "CC(C)OCOC(OCOC(C)C)OCOC(C)C",
        {"-CH3": 6, ">CH-": 4, "-CH2-": 3, "-O- (non-ring)": 6},
        "smiles",
    ),
    (
        "CC(C)OCOC(C)C",
        {"-CH3": 4, ">CH-": 2, "-CH2-": 1, "-O- (non-ring)": 2},
        "smiles",
    ),
    ("CCOCOCC", {"-CH3": 2, "-O- (non-ring)": 2, "-CH2-": 3}, "smiles"),
    ("COCOC", {"-CH3": 2, "-CH2-": 1, "-O- (non-ring)": 2}, "smiles"),
    # Problematics with acids
    (
        "COC(O)=O",
        {"-COOH (acid)": 1, "-CH3": 1, "-O- (non-ring)": 1},
        "smiles",
    ),
    (
        "CCOC(O)=O",
        {"-COOH (acid)": 1, "-CH2-": 1, "-O- (non-ring)": 1, "-CH3": 1},
        "smiles",
    ),
    (
        "CC(C)OC(O)=O",
        {"-COOH (acid)": 1, ">CH-": 1, "-O- (non-ring)": 1, "-CH3": 2},
        "smiles",
    ),
    (
        "CC(C)(C)OC(O)=O",
        {"-COOH (acid)": 1, ">C<": 1, "-CH3": 3, "-O- (non-ring)": 1},
        "smiles",
    ),
    (
        "OC(=O)OC1=CC=CC=C1",
        {"-COOH (acid)": 1, "ring=C<": 1, "ring=CH-": 5, "-O- (non-ring)": 1},
        "smiles",
    ),
    ("C1COCON1", {">NH (ring)": 1, "-O- (ring)": 2, "ring-CH2-": 3}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_ethers(identifier, result, identifier_type):
    assert get_groups(joback, identifier, identifier_type).subgroups == result

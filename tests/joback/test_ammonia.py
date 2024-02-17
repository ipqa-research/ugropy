import pytest

import ugropy as ug


# =============================================================================
# -NH2, >NH (non-ring), >NH (ring), >N- (non-ring), -N= (non-ring), -N= (ring),
# =NH, -CN, -NO2
# =============================================================================
# Joback
trials = [
    ("CCC=N", {"-CH3": 1, "-CH2-": 1, "=CH-": 1, "=NH": 1}, "smiles"),
    # methylamine
    ("CN", {"-CH3": 1, "-NH2": 1}, "smiles"),
    # isopropylamine
    ("CC(C)N", {"-CH3": 2, ">CH-": 1, "-NH2": 1}, "smiles"),
    # propylamine
    ("CCCN", {"-CH3": 1, "-CH2-": 2, "-NH2": 1}, "smiles"),
    # (+)-Stepharine (ALMOST: 1 bond modified)
    (
        "COC1=C(C2=C3C(CC24C=CC(=O)C=C4)NCCC3=C1)OC",
        {
            "-CH3": 2,
            "ring-CH2-": 3,
            "ring>CH-": 1,
            "ring>C<": 1,
            "ring=CH-": 5,
            "ring=C<": 5,
            "-O- (non-ring)": 2,
            ">C=O (ring)": 1,
            ">NH (ring)": 1,
        },
        "smiles",
    ),
    (
        "CC(C)N(CN)C(C)C",
        {"-CH3": 4, "-CH2-": 1, ">CH-": 2, ">N- (non-ring)": 1, "-NH2": 1},
        "smiles",
    ),
    (
        "CC(C)N(C)CN",
        {"-CH3": 3, ">CH-": 1, "-CH2-": 1, ">N- (non-ring)": 1, "-NH2": 1},
        "smiles",
    ),
    (
        "CC(C)NC(C)NC(C)(C)C",
        {">CH-": 2, ">NH (non-ring)": 2, "-CH3": 6, ">C<": 1},
        "smiles",
    ),
    (
        "CC(C)NC(C)N",
        {"-NH2": 1, ">CH-": 2, ">NH (non-ring)": 1, "-CH3": 3},
        "smiles",
    ),
    (
        "CCC(C)(C)NC(C)C",
        {"-CH3": 5, ">CH-": 1, ">NH (non-ring)": 1, "-CH2-": 1, ">C<": 1},
        "smiles",
    ),
    (
        "CCNC(C)CC",
        {"-CH3": 3, ">NH (non-ring)": 1, ">CH-": 1, "-CH2-": 2},
        "smiles",
    ),
    ("CCCNC", {">NH (non-ring)": 1, "-CH2-": 2, "-CH3": 2}, "smiles"),
    (
        "CN1CCCC1CC(=O)CC1CCCN1C",
        {},
        "smiles",
    ),
    # dimethylamine
    ("CNC", {"-CH3": 2, ">NH (non-ring)": 1}, "smiles"),
    # diethylamine
    ("CCNCC", {"-CH3": 2, "-CH2-": 2, ">NH (non-ring)": 1}, "smiles"),
    # diisopropylamine
    ("CC(C)NC(C)C", {"-CH3": 4, ">CH-": 2, ">NH (non-ring)": 1}, "smiles"),
    # Problematics
    (
        "CC(C)NCN",
        {">CH-": 1, ">NH (non-ring)": 1, "-CH3": 2, "-CH2-": 1, "-NH2": 1},
        "smiles",
    ),
    # Quinuclidine
    ("C1CN2CCC1CC2", {}, "smiles"),
    (
        "CCN(CC(=O)CC)C1=CC=CC=C1",
        {
            "-CH3": 2,
            ">N- (non-ring)": 1,
            "ring=C<": 1,
            "ring=CH-": 5,
            ">C=O (non-ring)": 1,
            "-CH2-": 3,
        },
        "smiles",
    ),
    (
        "CCN(C(C)C)C(C)C",
        {"-CH2-": 1, ">N- (non-ring)": 1, "-CH3": 5, ">CH-": 2},
        "smiles",
    ),
    ("CCN(C)CC", {">N- (non-ring)": 1, "-CH2-": 2, "-CH3": 3}, "smiles"),
    # trimethylamine
    ("CN(C)C", {"-CH3": 3, ">N- (non-ring)": 1}, "smiles"),
    # triethylamine
    ("CCN(CC)CC", {"-CH3": 3, "-CH2-": 3, ">N- (non-ring)": 1}, "smiles"),
    # aniline
    ("C1=CC=C(C=C1)N", {"ring=CH-": 5, "ring=C<": 1, "-NH2": 1}, "smiles"),
    # pyridine sharing side with benzene ring
    (
        "C1=CC2=C(C=C1)C1=C(C=NC=C1)C1=C2C=CC=C1",
        {"-N= (ring)": 1, "ring=CH-": 11, "ring=C<": 6},
        "smiles",
    ),
    # cyclic ester with pyridine
    (
        "O=C1CCC2=C(O1)C=CN=C2",
        {
            "-N= (ring)": 1,
            "ring=CH-": 3,
            "ring=C<": 2,
            "ring-CH2-": 2,
            "-COO- (ester)": 1,
        },
        "smiles",
    ),
    # Pyrydine sharin CH2 with phenyl
    (
        "CCCC1=CC=C(CC2=CC=NC=C2)C=C1",
        {"-N= (ring)": 1, "ring=C<": 3, "ring=CH-": 8, "-CH2-": 3, "-CH3": 1},
        "smiles",
    ),
    # Nicotine
    (
        "CN1CCCC1C2=CN=CC=C2",
        {},
        "smiles",
    ),
    (
        "CC1=CC(C)=C(C)C=N1",
        {"-CH3": 3, "ring=C<": 3, "ring=CH-": 2, "-N= (ring)": 1},
        "smiles",
    ),
    (
        "C1=CC=C(C=C1)C1=CC=NC=C1",
        {"-N= (ring)": 1, "ring=C<": 2, "ring=CH-": 9},
        "smiles",
    ),
    (
        "CC(=C)C1=CC=NC=C1",
        {
            "-N= (ring)": 1,
            "ring=C<": 1,
            "ring=CH-": 4,
            "=CH2": 1,
            "=C<": 1,
            "-CH3": 1,
        },
        "smiles",
    ),
    (
        "CC1=NC=CC(O)=C1",
        {
            "-N= (ring)": 1,
            "ring=C<": 2,
            "ring=CH-": 3,
            "-CH3": 1,
            "-OH (phenol)": 1,
        },
        "smiles",
    ),
    # pyridine
    ("C1=CC=NC=C1", {"-N= (ring)": 1, "ring=CH-": 5}, "smiles"),
    # 3-methylpyridine
    (
        "CC1=CN=CC=C1",
        {"-N= (ring)": 1, "ring=C<": 1, "ring=CH-": 4, "-CH3": 1},
        "smiles",
    ),
    # 2,3-Dimethylpyridine
    (
        "CC1=C(N=CC=C1)C",
        {"-N= (ring)": 1, "ring=C<": 2, "ring=CH-": 3, "-CH3": 2},
        "smiles",
    ),
    # acetonitrile
    ("CC#N", {"-CH3": 1, "-CN": 1}, "smiles"),
    # propionitrile
    ("CCC#N", {"-CH3": 1, "-CH2-": 1, "-CN": 1}, "smiles"),
    (
        "CCCC1=CC=C(C[N+]([O-])=O)C=C1",
        {"ring=CH-": 4, "ring=C<": 2, "-NO2": 1, "-CH2-": 3, "-CH3": 1},
        "smiles",
    ),
    (
        "[O-][N+](=O)CC1=CC=CC=C1",
        {"ring=CH-": 5, "ring=C<": 1, "-NO2": 1, "-CH2-": 1},
        "smiles",
    ),
    # nitromethane
    ("C[N+](=O)[O-]", {"-CH3": 1, "-NO2": 1}, "smiles"),
    # 1-nitropropane
    ("CCC[N+](=O)[O-]", {"-CH3": 1, "-CH2-": 2, "-NO2": 1}, "smiles"),
    # 2-nitropropane
    ("CC(C)[N+](=O)[O-]", {"-CH3": 2, ">CH-": 1, "-NO2": 1}, "smiles"),
    # nitrobenzene
    (
        "C1=CC=C(C=C1)[N+](=O)[O-]",
        {"ring=CH-": 5, "ring=C<": 1, "-NO2": 1},
        "smiles",
    ),
    (
        "[O-][N+](=O)C1=CC=NC=C1",
        {"ring=CH-": 4, "ring=C<": 1, "-NO2": 1, "-N= (ring)": 1},
        "smiles",
    ),
    # N-Methyl-2-pyrrolidone
    ("CN1CCCC1=O", {}, "smiles"),
    # Isophorone diisocyanate
    (
        "CC1(CC(CC(C1)(C)CN=C=O)N=C=O)C",
        {
            "-CH3": 3,
            "-CH2-": 1,
            "ring-CH2-": 3,
            "ring>CH-": 1,
            "ring>C<": 2,
            "=C=": 2,
            "=O (other than above)": 2,
            "-N= (non-ring)": 2,
        },
        "smiles",
    ),
    # morpholine
    ("C1COCCN1", {"ring-CH2-": 4, "-O- (ring)": 1, ">NH (ring)": 1}, "smiles"),
    # acetamide
    ("CC(=O)N", {"-CH3": 1, "-NH2": 1, ">C=O (non-ring)": 1}, "smiles"),
    # N-Methylacetamide
    (
        "CC(=O)NC",
        {"-CH3": 2, ">NH (non-ring)": 1, ">C=O (non-ring)": 1},
        "smiles",
    ),
    # N-Ethylacetamide
    (
        "CCNC(=O)C",
        {"-CH3": 2, ">NH (non-ring)": 1, ">C=O (non-ring)": 1, "-CH2-": 1},
        "smiles",
    ),
    # N,N-Dimethylacetamide
    (
        "CC(=O)N(C)C",
        {"-CH3": 3, ">N- (non-ring)": 1, ">C=O (non-ring)": 1},
        "smiles",
    ),
    # N-ethyl-N-methylacetamide
    (
        "CCN(C)C(=O)C",
        {"-CH3": 3, ">N- (non-ring)": 1, ">C=O (non-ring)": 1, "-CH2-": 1},
        "smiles",
    ),
    # N,N-Diethylacetamide
    (
        "CCN(CC)C(=O)C",
        {"-CH3": 3, ">N- (non-ring)": 1, ">C=O (non-ring)": 1, "-CH2-": 2},
        "smiles",
    ),
    # di amide + amine
    (
        "CCN(C(C)C)C(=O)NC(C)C",
        {
            "-CH3": 5,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            "-CH2-": 1,
            ">CH-": 2,
        },
        "smiles",
    ),
    (
        "CC(C)NC(=O)N(C)C(C)C",
        {
            "-CH3": 5,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 2,
        },
        "smiles",
    ),
    (
        "CCN(CC)C(=O)NC(C)C",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
            "-CH2-": 2,
        },
        "smiles",
    ),
    (
        "CCN(C)C(=O)NC(C)C",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    (
        "CC(C)NC(=O)N(C)C",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
        },
        "smiles",
    ),
    (
        "CCNC(=O)N(CC)C(C)C",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
            "-CH2-": 2,
        },
        "smiles",
    ),
    (
        "CCNC(=O)N(C)C(C)C",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    (
        "CCNC(=O)N(CC)CC",
        {
            "-CH3": 3,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            "-CH2-": 3,
        },
        "smiles",
    ),
    (
        "CCNC(=O)N(C)CC",
        {
            "-CH3": 3,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            "-CH2-": 2,
        },
        "smiles",
    ),
    (
        "CCNC(=O)N(C)C",
        {
            "-CH3": 3,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    (
        "CCN(C(C)C)C(=O)NC",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    (
        "CNC(=O)N(C)C(C)C",
        {
            "-CH3": 4,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
        },
        "smiles",
    ),
    (
        "CCN(CC)C(=O)NC",
        {
            "-CH3": 3,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            "-CH2-": 2,
        },
        "smiles",
    ),
    (
        "CCN(C)C(=O)NC",
        {
            "-CH3": 3,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    (
        "CNC(=O)N(C)C",
        {
            "-CH3": 3,
            ">N- (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">C=O (non-ring)": 1,
        },
        "smiles",
    ),
    (
        "CC(C)NC(=O)NC(C)C",
        {"-CH3": 4, ">NH (non-ring)": 2, ">C=O (non-ring)": 1, ">CH-": 2},
        "smiles",
    ),
    (
        "CCNC(=O)NC(C)C",
        {
            "-CH3": 3,
            ">NH (non-ring)": 2,
            ">C=O (non-ring)": 1,
            ">CH-": 1,
            "-CH2-": 1,
        },
        "smiles",
    ),
    (
        "CCNC(=O)NCC",
        {"-CH3": 2, ">NH (non-ring)": 2, ">C=O (non-ring)": 1, "-CH2-": 2},
        "smiles",
    ),
    (
        "CNC(=O)NC(C)C",
        {"-CH3": 3, ">NH (non-ring)": 2, ">C=O (non-ring)": 1, ">CH-": 1},
        "smiles",
    ),
    (
        "CCNC(=O)NC",
        {"-CH3": 2, ">NH (non-ring)": 2, ">C=O (non-ring)": 1, "-CH2-": 1},
        "smiles",
    ),
    (
        "CNC(=O)NC",
        {"-CH3": 2, ">NH (non-ring)": 2, ">C=O (non-ring)": 1},
        "smiles",
    ),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_ammonia(identifier, result, identifier_type):
    assert ug.get_groups(ug.joback, identifier, identifier_type) == result

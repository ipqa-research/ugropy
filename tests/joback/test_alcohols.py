import pytest

import ugropy as ug


# =============================================================================
# -OH (alcohol), -OH (phenol)
# =============================================================================
# Joback
trials = [
    # 1,2-Cyclohexanediol, 4-tert-butyl-1-phenyl-, stereoisomer
    (
        "CC(C)(C)C1CCC(C(C1)O)(C2=CC=CC=C2)O",
        {
            "-CH3": 3,
            ">C<": 1,
            "ring-CH2-": 3,
            "ring>CH-": 2,
            "ring>C<": 1,
            "-OH (alcohol)": 2,
            "ring=C<": 1,
            "ring=CH-": 5,
        },
        "smiles",
    ),
    # (2S,3S)-2-Methyl-1,3-hexanediol
    (
        "CCCC(C(C)CO)O",
        {"-CH3": 2, "-CH2-": 3, ">CH-": 2, "-OH (alcohol)": 2},
        "smiles",
    ),
    # 2-propanol
    ("CC(C)O", {"-CH3": 2, ">CH-": 1, "-OH (alcohol)": 1}, "smiles"),
    ("OC(O)=O", {"-COOH (acid)": 1, "-OH (alcohol)": 1}, "smiles"),
    # Phenanthrene-3,4-diol
    (
        "C1=CC=C2C(=C1)C=CC3=C2C(=C(C=C3)O)O",
        {"ring=CH-": 8, "ring=C<": 6, "-OH (phenol)": 2},
        "smiles",
    ),
    # 3-(tert-butyl)benzene-1,2-diol
    (
        "CC(C)(C)C1=C(C(=CC=C1)O)O",
        {"ring=CH-": 3, "ring=C<": 3, "-OH (phenol)": 2, "-CH3": 3, ">C<": 1},
        "smiles",
    ),
    # [1,1'-Biphenyl]-2,3',4-triol
    (
        "C1=CC(=CC(=C1)O)C2=C(C=C(C=C2)O)O",
        {"ring=CH-": 7, "ring=C<": 5, "-OH (phenol)": 3},
        "smiles",
    ),
    # phenol
    (
        "C1=CC=C(C=C1)O",
        {"ring=CH-": 5, "ring=C<": 1, "-OH (phenol)": 1},
        "smiles",
    ),
    # methanol
    ("CO", {"-CH3": 1, "-OH (alcohol)": 1}, "smiles"),
    # Gastrodigenin
    (
        "C1=CC(=CC=C1CO)O",
        {
            "ring=CH-": 4,
            "ring=C<": 2,
            "-OH (phenol)": 1,
            "-CH2-": 1,
            "-OH (alcohol)": 1,
        },
        "smiles",
    ),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_alcohols(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result

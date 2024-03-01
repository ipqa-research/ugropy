import pytest

from ugropy import get_groups, joback


# =============================================================================
# -COOH (acid)
# =============================================================================
# Joback
trials = [
    # 2,4-Diaminobutyric acid
    (
        "C(CN)C(C(=O)O)N",
        {"-COOH (acid)": 1, ">CH-": 1, "-NH2": 2, "-CH2-": 2},
        "smiles",
    ),
    ("OC(O)=O", {"-COOH (acid)": 1, "-OH (alcohol)": 1}, "smiles"),
    (
        "CCOC(O)=O",
        {"-COOH (acid)": 1, "-CH2-": 1, "-O- (non-ring)": 1, "-CH3": 1},
        "smiles",
    ),
    # acetic acid
    ("CC(=O)O", {"-CH3": 1, "-COOH (acid)": 1}, "smiles"),
    # formic acid
    ("C(=O)O", {}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_acids(identifier, result, identifier_type):
    assert get_groups(joback, identifier, identifier_type).subgroups == result

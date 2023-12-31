import pytest

import ugropy as ug


# =============================================================================
# -COOH (acid)
# =============================================================================
# Joback
trials = [
    ("OC(O)=O", {"-COOH (acid)": 1, "-OH (alcohol)": 1}, "smiles"),
    (
        "CCOC(O)=O",
        {"-COOH (acid)": 1, "-CH2-": 1, "-O- (non-ring)": 1, "-CH3": 1},
        "smiles",
    ),
    # 2,4-Diaminobutyric acid
    (
        "C(CN)C(C(=O)O)N",
        {"-COOH (acid)": 1, ">CH-": 1, "-NH2": 2, "-CH2-": 2},
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
    assert ug.get_joback_groups(identifier, identifier_type) == result

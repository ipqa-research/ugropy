import pytest

import ugropy as ug


# =============================================================================
# 5|OH|[14]OH (P) [81]OH (S) [82]OH (T)
# =============================================================================

# Dortmund
trials = [
    ("CCC#CO", {"CH3": 1, "CH2": 1, "C=-C": 1, "OH (P)": 1}, "smiles"),
    ("CCC=C(C)O", {"CH3": 2, "CH2": 1, "CH=C": 1, "OH (S)": 1}, "smiles"),
    (r"CC\C=C\O", {"CH3": 1, "CH2": 1, "CH=CH": 1, "OH (P)": 1}, "smiles"),
    (
        "CC(O)(CO)C(O)C#CO",
        {
            "CH3": 1,
            "CH2": 1,
            "CH": 1,
            "C": 1,
            "OH (P)": 2,
            "OH (S)": 1,
            "OH (T)": 1,
            "C=-C": 1,
        },
        "smiles",
    ),
    # (2S,3S)-2-Methyl-1,3-hexanediol
    (
        "CCCC(C(C)CO)O",
        {"CH3": 2, "CH2": 3, "CH": 2, "OH (P)": 1, "OH (S)": 1},
        "smiles",
    ),
    # 2-propanol
    ("CC(C)O", {"CH3": 2, "CH": 1, "OH (S)": 1}, "smiles"),
    ("OC(O)=O", {"COOH": 1, "OH (P)": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_oh_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

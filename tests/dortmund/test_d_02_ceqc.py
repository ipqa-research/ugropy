import pytest

import ugropy as ug


# =============================================================================
# 2|C=C|[5]CH2=CH [6]CH=CH [7]CH2=C [8]CH=C [70]C=C
# =============================================================================

# Dortmund
trials = [
    (
        "CC=C(C)C1=C(C=CC=C1C=C)C(C)=C(C)C",
        {"CH3": 5, "CH2=CH": 1, "CH=C": 1, "C=C": 1, "ACH": 3, "AC": 3},
        "smiles",
    ),
    (
        "CC=CC(C)=C(C)C=C",
        {"CH2=CH": 1, "CH=CH": 1, "C=C": 1, "CH3": 3},
        "smiles",
    ),
    ("CC(=C(C)C)C", {"C=C": 1, "CH3": 4}, "smiles"),  # 2,3-dimethylbutene-2
    ("CC=C(C)C", {"CH=C": 1, "CH3": 3}, "smiles"),  # 2-methyl-2-butene
    # 2-methyl-1-butene
    ("CCC(=C)C", {"CH2": 1, "CH2=C": 1, "CH3": 2}, "smiles"),
    ("CCCC=CC", {"CH2": 2, "CH=CH": 1, "CH3": 2}, "smiles"),  # 2-hexene
    ("CCCCC=C", {"CH2": 3, "CH2=CH": 1, "CH3": 1}, "smiles"),  # 1-hexene
    # impossibles
    ("C=C=C", {}, "smiles"),
    ("CC=CC(C)C(C)=C=C", {}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ceqc_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

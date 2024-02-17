import pytest

import ugropy as ug


# =============================================================================
# 1- 1|CH2|[1]CH3 [2]CH2 [3]CH [4]C
# =============================================================================
# UNIFAC
trials = [
    ("CC", {"CH3": 2}, "smiles"),  # ethane
    ("CCCCCC", {"CH3": 2, "CH2": 4}, "smiles"),  # hexane
    ("CC(C)C", {"CH3": 3, "CH": 1}, "smiles"),  # 2-methylpropane
    ("CC(C)(C)C", {"CH3": 4, "C": 1}, "smiles"),  # 2,2-dimethylpropane
    ("CCC(CC)C(C)(C)C", {"CH3": 5, "CH2": 2, "CH": 1, "C": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_dortmund_ch2(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

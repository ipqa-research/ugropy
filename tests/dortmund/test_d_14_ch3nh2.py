import pytest

import ugropy as ug


# =============================================================================
# 14|CH2NH2|[28]CH3NH2 [29]CH2NH2 [30]CHNH2 [85]CNH2
# =============================================================================

# Dortmund
trials_unifac = [
    # methylamine
    ("CN", {"CH3NH2": 1}, "smiles"),
    # isopropylamine
    ("CC(C)N", {"CH3": 2, "CHNH2": 1}, "smiles"),
    # propylamine
    ("CCCN", {"CH3": 1, "CH2": 1, "CH2NH2": 1}, "smiles"),
    ("CC(C)(C)N", {"CH3": 3, "CNH2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result


import pytest

import ugropy as ug


# =============================================================================
# 34|C=-C|[65]CH=-C [66]C=-C
# =============================================================================

# Dortmund
trials = [
    ("CC#CC1=CC=CC=C1", {"ACH": 5, "AC": 1, "C=-C": 1, "CH3": 1}, "smiles"),
    # 1-hexyne
    ("CCCCC#C", {"CH3": 1, "CH2": 3, "CH=-C": 1}, "smiles"),
    # 2-hexyne
    ("CCCC#CC", {"CH3": 2, "CH2": 2, "C=-C": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_alquine_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

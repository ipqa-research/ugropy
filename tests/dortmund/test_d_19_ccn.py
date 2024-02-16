import pytest

import ugropy as ug


# =============================================================================
# 19|CH2CN|[40]CH3CN [41]CH2CN
# =============================================================================

# Dortmund
trials = [
    # acetonitrile
    ("CC#N", {"CH3CN": 1}, "smiles"),
    # propionitrile
    ("CCC#N", {"CH3": 1, "CH2CN": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ccn_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

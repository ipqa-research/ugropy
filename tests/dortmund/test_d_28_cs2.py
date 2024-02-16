import pytest

import ugropy as ug


# =============================================================================
# 28|CS2|[58]CS2
# =============================================================================

# Dortmund
trials = [
    # carbon disulfide
    ("C(=S)=S", {"CS2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cs2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

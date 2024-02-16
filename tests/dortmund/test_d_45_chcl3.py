import pytest

import ugropy as ug


# =============================================================================
# 45|CHCL3|[50]CHCL3
# =============================================================================

# Dortmund
trials = [
    # chloroform
    ("C(Cl)(Cl)Cl", {"CHCL3": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_chcl3_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

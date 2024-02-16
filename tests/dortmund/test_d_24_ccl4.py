import pytest

import ugropy as ug


# =============================================================================
# 24|CCL4|[52]CCL4
# =============================================================================

# Dortmund
trials = [
    # tetrachloromethane
    ("C(Cl)(Cl)(Cl)Cl", {"CCL4": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ccl4_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

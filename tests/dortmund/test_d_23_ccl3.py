import pytest

import ugropy as ug


# =============================================================================
# 23|CCL3|[51]CCL3
# =============================================================================

# Dortmund
trials = [
    # trichloroethane
    ("CC(Cl)(Cl)Cl", {"CH3": 1, "CCL3": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ccl3_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

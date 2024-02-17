import pytest

import ugropy as ug


# =============================================================================
# 36|ACRY|[68]ACRY
# =============================================================================

# Dortmund
trials = [
    # acrylonitrile
    ("C=CC#N", {"ACRY": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_acry_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

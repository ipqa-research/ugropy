import pytest

import ugropy as ug


# =============================================================================
# 38|ACF|[71]ACF
# =============================================================================

# Dortmund
trials = [
    # hexafluorobenzene
    ("C1(=C(C(=C(C(=C1F)F)F)F)F)F", {"ACF": 6}, "smiles"),
    ("FC1=CC=NC=C1", {"AC2H2N": 1, "ACH": 2, "ACF": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_acf_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

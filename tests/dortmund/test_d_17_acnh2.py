import pytest

import ugropy as ug


# =============================================================================
# 17|ACNH2|[36]ACNH2
# =============================================================================

# Dortmund
trials = [
    # aniline
    ("C1=CC=C(C=C1)N", {"ACH": 5, "ACNH2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ch3nh2_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

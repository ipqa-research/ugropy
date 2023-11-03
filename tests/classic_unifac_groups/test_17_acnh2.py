import pytest

import ugropy as ug


# =============================================================================
# 17- ACNH2 Main group: ACNH2
# =============================================================================

# UNIFAC
trials_unifac = [
    # aniline
    ("C1=CC=C(C=C1)N", {"ACH": 5, "ACNH2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 28- CS2 Main group: CS2
# =============================================================================

# UNIFAC
trials_unifac = [
    # carbon disulfide
    ("C(=S)=S", {"CS2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cs2_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

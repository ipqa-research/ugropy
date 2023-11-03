import pytest

import ugropy as ug


# =============================================================================
# 25- ACCL Main group: ACCL
# =============================================================================

# UNIFAC
trials_unifac = [
    # chlorobenzene
    ("C1=CC=C(C=C1)Cl", {"ACCL": 1, "ACH": 5}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_accl_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

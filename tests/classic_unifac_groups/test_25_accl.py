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


@pytest.mark.ACCL
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_accl_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result
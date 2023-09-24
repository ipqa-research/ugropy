import pytest

import ugropy as ug


# =============================================================================
# 27- ACNO2 Main group: ACNO2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("nitrobenzene", {"ACH": 5, "ACNO2": 1}, "name"),
    ("[O-][N+](=O)C1=CC=NC=C1", {}, "smiles"),
]


@pytest.mark.ACNO2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acno2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

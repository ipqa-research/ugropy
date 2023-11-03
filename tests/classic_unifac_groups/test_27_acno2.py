import pytest

import ugropy as ug


# =============================================================================
# 27- ACNO2 Main group: ACNO2
# =============================================================================

# UNIFAC
trials_unifac = [
    # nitrobenzene
    ("C1=CC=C(C=C1)[N+](=O)[O-]", {"ACH": 5, "ACNO2": 1}, "smiles"),
    ("[O-][N+](=O)C1=CC=NC=C1", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acno2_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

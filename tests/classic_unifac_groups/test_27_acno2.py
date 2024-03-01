import pytest

from ugropy import get_groups, psrk, unifac


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
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acno2_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acno2_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result

import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 84- IMIDAZOL Main group: IMIDAZOL
# =============================================================================

# UNIFAC
trials_unifac = [
    # Imidazol
    ("N1C=CN=C1", {"IMIDAZOL": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result

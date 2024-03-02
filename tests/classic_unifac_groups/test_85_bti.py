import pytest

from ugropy import get_groups, unifac


# =============================================================================
# 85- BTI Main group: BTI
# =============================================================================

# UNIFAC
trials_unifac = [
    # BTI
    ("FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F", {"BTI": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result

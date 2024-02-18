import pytest

import ugropy as ug


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
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result

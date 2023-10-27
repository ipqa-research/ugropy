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


@pytest.mark.IMIDAZOL
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

import ugropy as ug

import pytest


# =============================================================================
# Impossible molecules
# =============================================================================

# UNIFAC
trials_unifac = [
    ("hydrogen peroxide", {}, "name"),
    ("methane", {}, "name"),
    ("C1(=CC=CC=C1)OC(C)(C)C", {}, "smiles"),
]

@pytest.mark.impossibles
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_impossibles_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
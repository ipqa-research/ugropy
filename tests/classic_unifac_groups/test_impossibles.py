import pytest

import ugropy as ug


# =============================================================================
# Impossible molecules
# =============================================================================

# UNIFAC
trials_unifac = [
    ("hydrogen peroxide", {}, "name"),
    ("methane", {}, "name"),
    ("C1(=CC=CC=C1)OC(C)(C)C", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_impossibles_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result


# PSRK
trials_psrk = [
    ("hydrogen peroxide", {}, "name"),
    ("C1(=CC=CC=C1)OC(C)(C)C", {}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_impossibles_psrk(identifier, result, identifier_type):
    assert ug.get_psrk_groups(identifier, identifier_type) == result

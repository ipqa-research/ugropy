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


@pytest.mark.impossibles
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_impossibles_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result


# PSRK
trials_psrk = [
    ("hydrogen peroxide", {}, "name"),
    ("C1(=CC=CC=C1)OC(C)(C)C", {}, "smiles"),
]


@pytest.mark.impossibles
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_impossibles_psrk(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.psrk_groups == result

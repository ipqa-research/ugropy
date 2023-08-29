import ugropy as ug

import pytest


# =============================================================================
# Impossible molecules
# =============================================================================

# UNIFAC
trials_unifac = [
    # ch2o
    ("C1COCON1", {}, "smiles"),
    # ceqc
    # propa‐1,2‐diene
    ("C=C=C", {}, "smiles"),
    ("CC=CC(C)C(C)=C=C", {}, "smiles"),
    # others
    ("hydrogen peroxide", {}, "name"),
    ("methane", {}, "name"),
    ("C1(=CC=CC=C1)OC(C)(C)C", {}, "smiles"),
]

@pytest.mark.impossibles
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_problematics_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
import pytest

import ugropy as ug


# =============================================================================
# 23- CCL3 Main group: CHCL3, CCL3
# =============================================================================

# UNIFAC
trials_unifac = [
    # chloroform
    ("C(Cl)(Cl)Cl", {"CHCL3": 1}, "smiles"),
    # trichloroethane
    ("CC(Cl)(Cl)Cl", {"CH3": 1, "CCL3": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl3_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result

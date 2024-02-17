import pytest

import ugropy as ug


# =============================================================================
# 24- CCL4 Main group: CCL4
# =============================================================================

# UNIFAC
trials_unifac = [
    # tetrachloromethane
    ("C(Cl)(Cl)(Cl)Cl", {"CCL4": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl4_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result

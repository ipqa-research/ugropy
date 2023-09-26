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


@pytest.mark.CCL4
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl4_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

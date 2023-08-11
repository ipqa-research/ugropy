import ugropy as ug

import pytest


# =============================================================================
# 41- COO Main group: COO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("Methyl acrylate", {"CH3": 1, "CH2=CH": 1, "COO": 1}, "name"),
]

@pytest.mark.COO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_coo_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
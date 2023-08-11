import ugropy as ug

import pytest


# =============================================================================
# 37- CL-(C=C) Main group: CL-(C=C)
# =============================================================================

# UNIFAC
trials_unifac = [
    ("trichloroethylene ", {"CH=C": 1, "CL-(C=C)": 3}, "name"),
]

@pytest.mark.CLCeqC
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cleqc_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
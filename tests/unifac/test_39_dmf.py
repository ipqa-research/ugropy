import ugropy as ug

import pytest


# =============================================================================
# 39- DMF Main group: DMF, HCON(..
# =============================================================================

# UNIFAC
trials_unifac = [
    (" N,N-Dimethylformamide", {"DMF": 1}, "name"),
    ("N,N-Diethylformamide", {"CH3": 2, "HCON(..": 1}, "name"),
]

@pytest.mark.DMF
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_dmf_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
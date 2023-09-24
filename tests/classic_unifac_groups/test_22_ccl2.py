import pytest

import ugropy as ug


# =============================================================================
# 22- CCL2 Main group: CH2CL2, CHCL2, CCL2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("OC(Cl)Cl", {"CHCL2": 1, "OH": 1}, "smiles"),
    ("dichloro methane", {"CH2CL2": 1}, "name"),
    ("1,1-dichloroethane", {"CH3": 1, "CHCL2": 1}, "name"),
    ("2,2-dichloropropane", {"CH3": 2, "CCL2": 1}, "name"),
]


@pytest.mark.CCL2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

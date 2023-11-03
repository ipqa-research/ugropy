import pytest

import ugropy as ug


# =============================================================================
# 22- CCL2 Main group: CH2CL2, CHCL2, CCL2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("OC(Cl)Cl", {"CHCL2": 1, "OH": 1}, "smiles"),
    # dichloro methane
    ("C(Cl)Cl", {"CH2CL2": 1}, "smiles"),
    # 1,1-dichloroethane
    ("CC(Cl)Cl", {"CH3": 1, "CHCL2": 1}, "smiles"),
    # 2,2-dichloropropane
    ("CC(C)(Cl)Cl", {"CH3": 2, "CCL2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl2_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

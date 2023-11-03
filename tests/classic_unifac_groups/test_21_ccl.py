import pytest

import ugropy as ug


# =============================================================================
# 21- CCL Main group: CH2CL, CHCL, CCL
# =============================================================================

# UNIFAC
trials_unifac = [
    # 1-chlorobutane
    ("CCCCCl", {"CH3": 1, "CH2": 2, "CH2CL": 1}, "smiles"),
    # 2-chloropropane
    ("CC(C)Cl", {"CH3": 2, "CHCL": 1}, "smiles"),
    # 2-chloro-2-methylpropane
    ("CC(C)(C)Cl", {"CH3": 3, "CCL": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result

import pytest

import ugropy as ug


# =============================================================================
# 22|CCL2|[47]CH2CL2 [48]CHCL2 [49]CCL2
# =============================================================================

# Dortmund
trials = [
    ("OC(Cl)Cl", {"CHCL2": 1, "OH (P)": 1}, "smiles"),
    # dichloro methane
    ("C(Cl)Cl", {"CH2CL2": 1}, "smiles"),
    # 1,1-dichloroethane
    ("CC(Cl)Cl", {"CH3": 1, "CHCL2": 1}, "smiles"),
    # 2,2-dichloropropane
    ("CC(C)(Cl)Cl", {"CH3": 2, "CCL2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ccl2_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

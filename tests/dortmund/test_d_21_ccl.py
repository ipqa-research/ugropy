import pytest

import ugropy as ug


# =============================================================================
# 21|CCL|[44]CH2CL [45]CHCL [46]CCL
# =============================================================================

# Dortmund
trials = [
    # 1-chlorobutane
    ("CCCCCl", {"CH3": 1, "CH2": 2, "CH2CL": 1}, "smiles"),
    # 2-chloropropane
    ("CC(C)Cl", {"CH3": 2, "CHCL": 1}, "smiles"),
    # 2-chloro-2-methylpropane
    ("CC(C)(C)Cl", {"CH3": 3, "CCL": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ccl_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result

import ugropy as ug

import pytest


# =============================================================================
# 13- CH2O Main group: CH2O, CH-O, THF
# =============================================================================

# UNIFAC
trials_unifac = [
    ("2H-Pyran, 2-(cyclohexyloxy)tetrahydro-", {"CH2": 8, "CH": 1, "CH2O": 1, "CH-O": 1}),
    # Problematic ones
    ("Benzyl 2-hydroxyethyl carbonate", {"CH2": 1, "OH": 1, "ACCH2": 1, "ACH": 5, "COO": 1, "CH2O": 1}),
    ("tert-Butyl ethyl carbonate", {"CH3": 4, "C": 1, "COO": 1, "CH2O": 1}),
    ("Ethyl phenyl carbonate", {"CH3": 1, "AC": 1, "ACH": 5, "COO": 1, "CH2O": 1}),
    ("Carbonic acid, ethyl 2,3,6-trimethylcyclohexyl ester", {"CH3": 4, "CH2": 2, "CH": 4, "COO": 1, "CH2O": 1}),
    ("Diethyl carbonate", {"CH3": 2, "CH2": 1, "COO": 1, "CH2O": 1}),
    ("SCHEMBL3938062", {"CH2": 4, "CH": 2, "OH": 1, "COO": 1, "CH3O": 1}),
    ("Methyl phenyl carbonate", {"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1}),
    ("tert-Butyl methyl carbonate", {"CH3": 3, "C": 1, "COO": 1, "CH3O": 1}),
    ("Methyl isopropyl carbonate", {"CH3": 2, "CH": 1, "COO": 1, "CH3O": 1}),
    ("Ethyl methyl carbonate", {"CH3": 1, "CH2": 1, "COO": 1, "CH3O": 1}),
    ("Dimethyl carbonate", {"CH3": 1, "COO": 1, "CH3O": 1}),
]

@pytest.mark.CH2O
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_cho_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result
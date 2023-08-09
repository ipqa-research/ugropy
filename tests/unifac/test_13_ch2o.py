import ugropy as ug

import pytest


# =============================================================================
# 13- CH2O Main group: CH2O, CH-O, THF
# =============================================================================

# UNIFAC
trials_unifac = [
    ("diisopropyl ether", {"CH3": 4, "CH": 1, "CH-O": 1}, "name"),
    ("diethyl ether", {"CH3": 2, "CH2": 1, "CH2O": 1}, "name"),
    ("dimethyl ether", {"CH3": 1, "CH3O": 1}, "name"),
    ("2H-Pyran, 2-(cyclohexyloxy)tetrahydro-", {"CH2": 8, "CH": 1, "CH2O": 1, "CH-O": 1}, "name"),
    # Problematic ones
    ("Benzyl 2-hydroxyethyl carbonate", {"CH2": 1, "OH": 1, "ACCH2": 1, "ACH": 5, "COO": 1, "CH2O": 1}, "name"),
    ("tert-Butyl ethyl carbonate", {"CH3": 4, "C": 1, "COO": 1, "CH2O": 1}, "name"),
    ("Ethyl phenyl carbonate", {"CH3": 1, "AC": 1, "ACH": 5, "COO": 1, "CH2O": 1}, "name"),
    ("Carbonic acid, ethyl 2,3,6-trimethylcyclohexyl ester", {"CH3": 4, "CH2": 2, "CH": 4, "COO": 1, "CH2O": 1}, "name"),
    ("Diethyl carbonate", {"CH3": 2, "CH2": 1, "COO": 1, "CH2O": 1}, "name"),
    ("SCHEMBL3938062", {"CH2": 4, "CH": 2, "OH": 1, "COO": 1, "CH3O": 1}, "name"),
    ("Methyl phenyl carbonate", {"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1}, "name"),
    ("tert-Butyl methyl carbonate", {"CH3": 3, "C": 1, "COO": 1, "CH3O": 1}, "name"),
    ("Methyl isopropyl carbonate", {"CH3": 2, "CH": 1, "COO": 1, "CH3O": 1}, "name"),
    ("Ethyl methyl carbonate", {"CH3": 1, "CH2": 1, "COO": 1, "CH3O": 1}, "name"),
    ("Dimethyl carbonate", {"CH3": 1, "COO": 1, "CH3O": 1}, "name"),
]

@pytest.mark.CH2O
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2o_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
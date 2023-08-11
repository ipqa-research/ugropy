import ugropy as ug

import pytest


# =============================================================================
# 13- CH2O Main group: CH2O, CH-O, THF
# =============================================================================

# UNIFAC
trials_unifac = [
    ("tetrahydrofuran", {"THF": 1}, "name"),
    ("diisopropyl ether", {"CH3": 4, "CH": 1, "CH-O": 1}, "name"),
    ("diethyl ether", {"CH3": 2, "CH2": 1, "CH2O": 1}, "name"),
    ("dimethyl ether", {"CH3": 1, "CH3O": 1}, "name"),
    ("SCHEMBL3938062", {"CH2": 4, "CH": 2, "OH": 1, "COO": 1, "CH3O": 1}, "name"),
    ("2H-Pyran, 2-(cyclohexyloxy)tetrahydro-", {"CH2": 8, "CH": 1, "CH2O": 1, "CH-O": 1}, "name"),
    # Problematic ones
    # Benzyl 2-hydroxyethyl carbonate
    ("C1=CC=C(C=C1)COC(=O)OCCO", {"ACCH2": 1, "ACH": 5, "COO": 1, "C2H5O2": 1}, "smiles"),
    # tert-Butyl ethyl carbonate
    ("CCOC(=O)OC(C)(C)C", {"CH3": 4, "C": 1, "COO": 1, "CH2O": 1}, "smiles"),
    # Ethyl phenyl carbonate
    ("CCOC(=O)OC1=CC=CC=C1", {"CH3": 1, "AC": 1, "ACH": 5, "COO": 1, "CH2O": 1}, "smiles"),
    # Carbonic acid, ethyl 2,3,6-trimethylcyclohexyl ester
    ("CCOC(=O)OC1C(CCC(C1C)C)C", {"CH3": 4, "CH2": 2, "CH": 4, "COO": 1, "CH2O": 1}, "smiles"),
    # Diethyl carbonate
    ("CCOC(=O)OCC", {"CH3": 2, "CH2": 1, "COO": 1, "CH2O": 1}, "smiles"),
    # Methyl phenyl carbonate
    ("COC(=O)OC1=CC=CC=C1", {"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1}, "smiles"),
    # tert-Butyl methyl carbonate
    ("CC(C)(C)OC(=O)OC", {"CH3": 3, "C": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # Methyl isopropyl carbonate
    ("CC(C)OC(=O)OC", {"CH3": 2, "CH": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # Ethyl methyl carbonate
    ("CCOC(=O)OC", {"CH3": 1, "CH2": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # Dimethyl carbonate
    ("COC(=O)OC", {"CH3": 1, "COO": 1, "CH3O": 1}, "smiles"),
]

@pytest.mark.CH2O
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2o_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
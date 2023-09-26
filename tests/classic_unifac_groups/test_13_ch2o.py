import pytest

import ugropy as ug


# =============================================================================
# 13- CH2O Main group: CH3O, CH2O, CH-O, THF
# =============================================================================

# UNIFAC
trials_unifac = [
    # 4-flavanol
    (
        "OC1CC(OC2=CC=CC=C12)C1=CC=CC=C1",
        {"ACH": 9, "AC": 2, "ACCH": 1, "CH-O": 1, "CH2": 1, "OH": 1},
        "smiles",
    ),
    (
        "O[C@@H]1CO[C@H](O)[C@@H](O)[C@@H]1O",
        {"CH2O": 1, "CH": 4, "OH": 4},
        "smiles",
    ),
    ("C1COCCOCCOCCOC1", {"CH2O": 4, "CH2": 5}, "smiles"),
    ("C1COCCO1", {"CH2O": 2, "CH2": 2}, "smiles"),
    ("CCOCOCC", {"CH3": 2, "CH2O": 2, "CH2": 1}, "smiles"),
    ("C1COCO1", {"CH2O": 2, "CH2": 1}, "smiles"),
    ("C1COCCOCCOCCOCCOCCO1", {"CH2O": 6, "CH2": 6}, "smiles"),
    # tetrahydrofuran
    ("C1CCOC1", {"THF": 1}, "smiles"),
    # diisopropyl ether
    ("CC(C)OC(C)C", {"CH3": 4, "CH": 1, "CH-O": 1}, "smiles"),
    # diethyl ether
    ("CCOCC", {"CH3": 2, "CH2": 1, "CH2O": 1}, "smiles"),
    # dimethyl ether
    ("COC", {"CH3": 1, "CH3O": 1}, "smiles"),
    # 2H-Pyran, 2-(cyclohexyloxy)tetrahydro-
    (
        "C1CCC(CC1)OC2CCCCO2",
        {"CH2": 8, "CH": 1, "CH2O": 1, "CH-O": 1},
        "smiles",
    ),
    # Problematic ones
    # Benzyl 2-hydroxyethyl carbonate
    (
        "C1=CC=C(C=C1)COC(=O)OCCO",
        {"ACCH2": 1, "ACH": 5, "COO": 1, "C2H5O2": 1},
        "smiles",
    ),
    # tert-Butyl ethyl carbonate
    ("CCOC(=O)OC(C)(C)C", {"CH3": 4, "C": 1, "COO": 1, "CH2O": 1}, "smiles"),
    # Ethyl phenyl carbonate
    (
        "CCOC(=O)OC1=CC=CC=C1",
        {"CH3": 1, "AC": 1, "ACH": 5, "COO": 1, "CH2O": 1},
        "smiles",
    ),
    # Carbonic acid, ethyl 2,3,6-trimethylcyclohexyl ester
    (
        "CCOC(=O)OC1C(CCC(C1C)C)C",
        {"CH3": 4, "CH2": 2, "CH": 4, "COO": 1, "CH2O": 1},
        "smiles",
    ),
    # Diethyl carbonate
    ("CCOC(=O)OCC", {"CH3": 2, "CH2": 1, "COO": 1, "CH2O": 1}, "smiles"),
    # Methyl phenyl carbonate
    (
        "COC(=O)OC1=CC=CC=C1",
        {"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1},
        "smiles",
    ),
    # tert-Butyl methyl carbonate
    ("CC(C)(C)OC(=O)OC", {"CH3": 3, "C": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # Methyl isopropyl carbonate
    ("CC(C)OC(=O)OC", {"CH3": 2, "CH": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # Ethyl methyl carbonate
    ("CCOC(=O)OC", {"CH3": 1, "CH2": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # Dimethyl carbonate
    ("COC(=O)OC", {"CH3": 1, "COO": 1, "CH3O": 1}, "smiles"),
    # I hate ether group
    ("COCOC(C)OCOC", {"CH3O": 2, "CH2O": 2, "CH": 1, "CH3": 1}, "smiles"),
    (
        "CC(C)OCOC(C)OCOC(C)C",
        {"CH3": 5, "CH": 1, "CHO": 2, "CH2O": 2},
        "smiles",
    ),
    (
        "CC(C)OCOCC(OCOC(C)C)OCOC(C)C",
        {"CH3": 6, "CH": 2, "CH2O": 4, "CHO": 2},
        "smiles",
    ),
    (
        "CC(C)OCOC(OCOC(C)C)OCOC(C)C",
        {"CH3": 6, "CHO": 3, "CH2O": 3, "CH": 1},
        "smiles",
    ),
    ("CC(C)OCOC(C)C", {"CH3": 4, "CHO": 1, "CH2O": 1, "CH": 1}, "smiles"),
    ("CCOCOCC", {"CH3": 2, "CH2O": 2, "CH2": 1}, "smiles"),
    ("COCOC", {"CH3O": 2, "CH2": 1}, "smiles"),
    # Impossibles
    ("C1COCON1", {}, "smiles"),
]


@pytest.mark.CH2O
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2o_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

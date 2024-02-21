import pytest

import ugropy as ug


# =============================================================================
# 13|CH2O|[24]CH3O [25]CH2O [26]CHO
# =============================================================================

# Dortmund
trials = [
    ("CCOCOCC", {"CH3": 2, "CH2O": 2, "CH2": 1}, "smiles"),
    # diisopropyl ether
    ("CC(C)OC(C)C", {"CH3": 4, "CH": 1, "CHO": 1}, "smiles"),
    # diethyl ether
    ("CCOCC", {"CH3": 2, "CH2": 1, "CH2O": 1}, "smiles"),
    # dimethyl ether
    ("COC", {"CH3": 1, "CH3O": 1}, "smiles"),
    # Problematic ones
    (
        "COC(=O)OC1=CC=CC=C1",
        {"CH3O": 1, "COO": 1, "AC": 1, "ACH": 5},
        "smiles",
    ),
    (
        "CCOC(=O)OC1=CC=CC=C1",
        {"CH3": 1, "CH2O": 1, "COO": 1, "AC": 1, "ACH": 5},
        "smiles",
    ),
    (
        "CC(C)OC(=O)OC1=CC=CC=C1",
        {"CH3": 2, "CHO": 1, "COO": 1, "AC": 1, "ACH": 5},
        "smiles",
    ),
    ("CC(C)(C)OC(=O)OC1=CC=CC=C1", {}, "smiles"),
    # TODO Benzyl 2-hydroxyethyl carbonate
    (
        "C1=CC=C(C=C1)COC(=O)OCCO",
        {"CH2": 1, "ACH": 5, "ACCH2": 1, "OH (P)": 1, "CH2O": 1, "COO": 1},
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
    # Problematics with acids
    ("COC(O)=O", {"COOH": 1, "CH3O": 1}, "smiles"),
    ("CCOC(O)=O", {"COOH": 1, "CH2O": 1, "CH3": 1}, "smiles"),
    ("CC(C)OC(O)=O", {"COOH": 1, "CHO": 1, "CH3": 2}, "smiles"),
    ("CC(C)(C)OC(O)=O", {"OH (P)": 1, "COO": 1, "C": 1, "CH3": 3}, "smiles"),
    (
        "OC(=O)OC1=CC=CC=C1",
        {"OH (P)": 1, "COO": 1, "AC": 1, "ACH": 5},
        "smiles",
    ),
    # Impossibles
    ("C1COCON1", {}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ch2o_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

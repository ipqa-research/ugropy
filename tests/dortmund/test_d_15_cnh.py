import pytest

import ugropy as ug


# =============================================================================
# 15|CH2NH|[31]CH3NH [32]CH2NH [33]CHNH
# =============================================================================

# Dortmund
trials = [
    # Dihydro-beta-carboline
    ("C1C2=C(C=CN1)C3=CC=CC=C3N2", {}, "smiles"),
    ("CC(C)N(CN)C(C)C", {}, "smiles"),
    ("CC(C)N(C)CN", {"CH3": 2, "CH": 1, "CH2NH2": 1, "CH3N": 1}, "smiles"),
    ("CC(C)NC(C)NC(C)(C)C", {"CHNH": 2, "CH3": 6, "C": 1}, "smiles"),
    ("CC(C)NC(C)N", {"CHNH2": 1, "CHNH": 1, "CH3": 3}, "smiles"),
    ("CCC(C)(C)NC(C)C", {"CH3": 5, "CHNH": 1, "CH2": 1, "C": 1}, "smiles"),
    ("CCNC(C)CC", {"CH3": 3, "CH2NH": 1, "CH": 1, "CH2": 1}, "smiles"),
    ("CCCNC", {"CH3NH": 1, "CH2": 2, "CH3": 1}, "smiles"),
    # dimethylamine
    ("CNC", {"CH3": 1, "CH3NH": 1}, "smiles"),
    # diethylamine
    ("CCNCC", {"CH3": 2, "CH2": 1, "CH2NH": 1}, "smiles"),
    # diisopropylamine
    ("CC(C)NC(C)C", {"CH3": 4, "CH": 1, "CHNH": 1}, "smiles"),
    # Problematics
    ("CC(C)NCN", {"CHNH": 1, "CH3": 2, "CH2NH2": 1}, "smiles"),
    # concatenate amine
    # TODO
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ch3nh2_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result

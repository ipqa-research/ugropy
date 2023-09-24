import pytest

import ugropy as ug


# =============================================================================
# 46- CON Main group: AMH2, AMHCH3, AMHCH2, AM(CH3)2, AMCH3CH2, AM(CH2)2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("acetamide", {"CH3": 1, "AMH2": 1}, "name"),
    ("N-Methylacetamide", {"CH3": 1, "AMHCH3": 1}, "name"),
    ("N-Ethylacetamide", {"CH3": 2, "AMHCH2": 1}, "name"),
    ("N,N-Dimethylacetamide", {"CH3": 1, "AM(CH3)2": 1}, "name"),
    ("N-ethyl-N-methylacetamide", {"CH3": 2, "AMCH3CH2": 1}, "name"),
    ("N,N-Diethylacetamide", {"CH3": 3, "AM(CH2)2": 1}, "name"),
    # di amide + amine
    ("CCN(C(C)C)C(=O)NC(C)C", {}, "smiles"),
    ("CC(C)NC(=O)N(C)C(C)C", {}, "smiles"),
    ("CCN(CC)C(=O)NC(C)C", {"AM(CH2)2": 1, "CHNH": 1, "CH3": 4}, "smiles"),
    ("CCN(C)C(=O)NC(C)C", {"AMCH3CH2": 1, "CHNH": 1, "CH3": 3}, "smiles"),
    ("CC(C)NC(=O)N(C)C", {"AM(CH3)2": 1, "CHNH": 1, "CH3": 2}, "smiles"),
    (
        "CCNC(=O)N(CC)C(C)C",
        {"AMHCH2": 1, "CH2N": 1, "CH3": 4, "CH": 1},
        "smiles",
    ),
    (
        "CCNC(=O)N(C)C(C)C",
        {"AMHCH2": 1, "CH3N": 1, "CH3": 3, "CH": 1},
        "smiles",
    ),
    ("CCNC(=O)N(CC)CC", {"AM(CH2)2": 1, "CH2NH": 1, "CH3": 3}, "smiles"),
    ("CCNC(=O)N(C)CC", {"AMCH3CH2": 1, "CH2NH": 1, "CH3": 2}, "smiles"),
    ("CCNC(=O)N(C)C", {"AM(CH3)2": 1, "CH2NH": 1, "CH3": 1}, "smiles"),
    (
        "CCN(C(C)C)C(=O)NC",
        {"AMHCH3": 1, "CH2N": 1, "CH": 1, "CH3": 3},
        "smiles",
    ),
    (
        "CNC(=O)N(C)C(C)C",
        {"AMHCH3": 1, "CH3N": 1, "CH": 1, "CH3": 2},
        "smiles",
    ),
    ("CCN(CC)C(=O)NC", {"AM(CH2)2": 1, "CH3NH": 1, "CH3": 2}, "smiles"),
    ("CCN(C)C(=O)NC", {"AMCH3CH2": 1, "CH3NH": 1, "CH3": 1}, "smiles"),
    ("CNC(=O)N(C)C", {"AM(CH3)2": 1, "CH3NH": 1}, "smiles"),
    ("CC(C)NC(=O)NC(C)C", {}, "smiles"),
    ("CCNC(=O)NC(C)C", {"AMHCH2": 1, "CHNH": 1, "CH3": 3}, "smiles"),
    ("CCNC(=O)NCC", {"AMHCH2": 1, "CH2NH": 1, "CH3": 2}, "smiles"),
    ("CNC(=O)NC(C)C", {"AMHCH3": 1, "CHNH": 1, "CH3": 2}, "smiles"),
    ("CCNC(=O)NC", {"AMHCH3": 1, "CH2NH": 1, "CH3": 1}, "smiles"),
    ("CNC(=O)NC", {"AMHCH3": 1, "CH3NH": 1}, "smiles"),
]


@pytest.mark.CON
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_con_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result

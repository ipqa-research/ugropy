import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 15- CNH Main group: CH3NH, CH2NH, CHNH
# =============================================================================

# UNIFAC
trials_unifac = [
    # (+)-Stepharine (ALMOST: 1 bond modified)
    (
        "COC1=C(OC)C2=C3C(CC22CCC(=O)C=C2)NCCC3=C1",
        {
            "CH3O": 2,
            "ACH": 1,
            "ACCH2": 1,
            "ACCH": 1,
            "AC": 3,
            "CH2NH": 1,
            "CH2": 2,
            "C": 1,
            "CH2CO": 1,
            "CH=CH": 1,
        },
        "smiles",
    ),
    # Dihydro-beta-carboline
    ("C1C2=C(C=CN1)C3=CC=CC=C3N2", {}, "smiles"),
    ("CC(C)N(CN)C(C)C", {}, "smiles"),
    ("CC(C)N(C)CN", {"CH3": 2, "CH": 1, "CH2NH2": 1, "CH3N": 1}, "smiles"),
    ("CC(C)NC(C)NC(C)(C)C", {"CHNH": 2, "CH3": 6, "C": 1}, "smiles"),
    ("CC(C)NC(C)N", {"CHNH2": 1, "CHNH": 1, "CH3": 3}, "smiles"),
    ("CCC(C)(C)NC(C)C", {"CH3": 5, "CHNH": 1, "CH2": 1, "C": 1}, "smiles"),
    ("CCNC(C)CC", {"CH3": 3, "CH2NH": 1, "CH": 1, "CH2": 1}, "smiles"),
    ("CCCNC", {"CH3NH": 1, "CH2": 2, "CH3": 1}, "smiles"),
    (
        "CN1CCCC1CC(=O)CC1CCCN1C",
        {"CH3N": 2, "CH2": 7, "CH": 2, "CH2CO": 1},
        "smiles",
    ),
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


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}

import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 11- CCOO Main group: CH3COO, CH2COO
# =============================================================================

# UNIFAC
trials_unifac = [
    # Dilaurin
    (
        "CCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCCCCCC",
        {"CH3": 2, "CH2": 20, "CH2COO": 2, "OH": 1, "CH": 1},
        "smiles",
    ),
    # Aspirin
    (
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        {"CH3COO": 1, "AC": 2, "ACH": 4, "COOH": 1},
        "smiles",
    ),
    # Tert-butyl acetate
    ("CC(=O)OC(C)(C)C", {"CH3COO": 1, "CH3": 3, "C": 1}, "smiles"),
    # triacetin
    ("CC(=O)OCC(COC(=O)C)OC(=O)C", {"CH3COO": 3, "CH2": 2, "CH": 1}, "smiles"),
    # butyl propanoate
    ("CCCCOC(=O)CC", {"CH3": 2, "CH2": 3, "CH2COO": 1}, "smiles"),
    # butyl acetate
    ("CCCCOC(=O)C", {"CH3": 1, "CH2": 3, "CH3COO": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccoo_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccoo_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccoo_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    assert mol.subgroups == result
    assert (
        fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary)
        != {}
    )

import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 48- CH2S Main group: CH3S, CH2S, CHS
# =============================================================================

# UNIFAC
trials_unifac = [
    # Dimethylsulfide
    ("CSC", {"CH3": 1, "CH3S": 1}, "smiles"),
    # Diethylsulfide
    ("CCSCC", {"CH3": 2, "CH2": 1, "CH2S": 1}, "smiles"),
    # Isopropyl Sulfide
    ("CC(C)SC(C)C", {"CH3": 4, "CH": 1, "CHS": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2s_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2s_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2s_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary) != {}

import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 22- CCL2 Main group: CH2CL2, CHCL2, CCL2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("OC(Cl)Cl", {"CHCL2": 1, "OH": 1}, "smiles"),
    # dichloro methane
    ("C(Cl)Cl", {"CH2CL2": 1}, "smiles"),
    # 1,1-dichloroethane
    ("CC(Cl)Cl", {"CH3": 1, "CHCL2": 1}, "smiles"),
    # 2,2-dichloropropane
    ("CC(C)(Cl)Cl", {"CH3": 2, "CCL2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl2_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl2_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccl2_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    
    if identifier != "C(Cl)Cl":
        assert mol.subgroups == result
        assert fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary) != {}
    else:
        assert mol.subgroups == {}

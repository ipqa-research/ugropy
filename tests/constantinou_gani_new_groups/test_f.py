import pytest

from ugropy import constantinou_gani_primary, get_groups
from ugropy.core import fit_atoms


trials_cg = [
    ("FC1CC(F)C(F)C1", {"CH2": 2, "CH": 3, "F": 3}, "smiles"),
    ("OC(F)=O", {"COOH": 1, "F": 1}, "smiles"),
    ("FC(F)(F)F", {"C": 1, "F": 4}, "smiles"),
]


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_cg)
def test_f_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    assert mol.subgroups == result
    assert (
        fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary)
        != {}
    )

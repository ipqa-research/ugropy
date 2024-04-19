import pytest

from ugropy import constantinou_gani_primary, get_groups
from ugropy.core import fit_atoms


trials_cg = [
    ("CC=C=C", {"CH2=C=CH": 1, "CH3": 1}, "smiles"),
    ("C=C=CC1=CC=CC=C1", {"CH2=C=CH": 1, "AC": 1, "ACH": 5}, "smiles"),
    ("C=C=CC=C=C", {"CH2=C=CH": 2}, "smiles"),
    ("FC=C=C", {"CH2=C=CH": 1, "F": 1}, "smiles"),
]


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_cg)
def test_alkene_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary) != {}

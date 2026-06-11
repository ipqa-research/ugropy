import pytest

from rdkit import Chem

from ugropy import instantiate_mol_object


@pytest.mark.not_github
def test_making_it_explode():
    with pytest.raises(ValueError):
        instantiate_mol_object("snake", "mol")

    with pytest.raises(ValueError):
        instantiate_mol_object("acetone", "Argentina")


@pytest.mark.not_github
def test_instante_from_name():
    mol_name = instantiate_mol_object("ethanol", "name")
    mol_smiles = instantiate_mol_object("CCO", "smiles")

    assert Chem.MolToSmiles(mol_name) == Chem.MolToSmiles(mol_smiles)

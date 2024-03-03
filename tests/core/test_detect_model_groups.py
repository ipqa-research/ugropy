import pytest

from ugropy import unifac
from ugropy.core import instantiate_mol_object
from ugropy.core.detect_model_groups import group_matches


def test_making_it_explode():
    with pytest.raises(ValueError):
        mol = instantiate_mol_object("CCC", "smiles")
        group_matches(mol_object=mol, group="CH3", model=unifac, action="cook")

import pytest

from ugropy import instantiate_mol_object


def test_making_it_explode():
    with pytest.raises(ValueError):
        instantiate_mol_object("snake", "mol")

    with pytest.raises(ValueError):
        instantiate_mol_object("acetone", "Argentina")
